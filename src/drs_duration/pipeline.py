from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import hilbert

from .io import load_record_txt
from .spectrum import compute_drs
from .window import detect_window_by_derivative
from .fit import fit_sinusoidal_in_window

from .plots import (
    save_drs_full_png,
    save_drs_zoom_png,
    save_final_sine_fit_png,
)


@dataclass
class PipelineConfig:
    zeta: float = 0.05
    tmin: float | None = None
    tmax: float = 15.0
    dT: float = 0.01

    zoom_frac: float = 0.20
    zoom_min_points: int = 15
    zoom_fallback_frac: float = 0.35

    alpha_operational: float = 0.30

    alpha_search_min: float = 0.20
    alpha_search_max: float = 0.50
    alpha_search_step: float = 0.01

    min_R2: float = 0.95
    min_points: int = 10
    min_td: float = 1.0

    allow_envelope_fallback: bool = True
    envelope_rel: float = 0.85


def _auto_zoom(T_values, RESP, zoom_frac=0.20, zoom_min_points=15, zoom_fallback_frac=0.35):
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    idx_peak = int(np.nanargmax(RESP))
    Tp = float(T_values[idx_peak])
    Dmax = float(RESP[idx_peak])

    Tmin_full = float(T_values[0])
    Tmax_full = float(T_values[-1])

    Tmin_zoom = max(Tmin_full, Tp * (1 - zoom_frac))
    Tmax_zoom = min(Tmax_full, Tp * (1 + zoom_frac))

    mask_zoom = (T_values >= Tmin_zoom) & (T_values <= Tmax_zoom)
    T_fit = T_values[mask_zoom]
    RESP_fit = RESP[mask_zoom]

    zoom_used = float(zoom_frac)

    if len(T_fit) < zoom_min_points:
        zoom_used = float(zoom_fallback_frac)
        Tmin_zoom = max(Tmin_full, Tp * (1 - zoom_used))
        Tmax_zoom = min(Tmax_full, Tp * (1 + zoom_used))
        mask_zoom = (T_values >= Tmin_zoom) & (T_values <= Tmax_zoom)
        T_fit = T_values[mask_zoom]
        RESP_fit = RESP[mask_zoom]

    return Tp, Dmax, Tmin_zoom, Tmax_zoom, T_fit, RESP_fit, zoom_used


def _is_physically_valid(fit_result, cfg: PipelineConfig) -> bool:
    if fit_result is None:
        return False

    R2 = fit_result.get("R2", np.nan)
    n_points = fit_result.get("n_points", 0)
    td = fit_result.get("td", np.nan)

    if not np.isfinite(R2) or not np.isfinite(td):
        return False

    return (
        float(R2) >= cfg.min_R2
        and int(n_points) >= cfg.min_points
        and float(td) > cfg.min_td
    )


def _evaluate_alpha(T_fit, RESP_fit, alpha: float, cfg: PipelineConfig):
    det = detect_window_by_derivative(T_fit, RESP_fit, alpha=alpha, plot=False)

    fit_res = fit_sinusoidal_in_window(
        T_fit,
        RESP_fit,
        detect_res=det,
        use_detect_function=False,
        verbose=False,
    )

    fit_result = fit_res["result"]
    valid = _is_physically_valid(fit_result, cfg)

    return det, fit_res, fit_result, valid


def _restricted_alpha_search(T_fit, RESP_fit, cfg: PipelineConfig):
    alphas = np.arange(
        cfg.alpha_search_min,
        cfg.alpha_search_max + 0.5 * cfg.alpha_search_step,
        cfg.alpha_search_step,
    )

    rows = []

    for alpha in alphas:
        det, fit_res, fit_result, valid = _evaluate_alpha(T_fit, RESP_fit, float(alpha), cfg)

        row = {
            "alpha": float(alpha),
            "valid": bool(valid),
            "T_left": float(det["T_left"]) if det else np.nan,
            "T_right": float(det["T_right"]) if det else np.nan,
            "window_method": det.get("method", "") if det else "",
            "n_points_window": int(fit_result["n_points"]) if fit_result else np.nan,
            "td": float(fit_result["td"]) if fit_result else np.nan,
            "R2": float(fit_result["R2"]) if fit_result else np.nan,
            "cycles": float(fit_result["cycles"]) if fit_result else np.nan,
            "ssr": float(fit_result["ssr"]) if fit_result else np.nan,
            "best_seed": str(fit_result["best_seed"]) if fit_result else "",
        }
        rows.append(row)

        if valid:
            return float(alpha), det, fit_res, fit_result, pd.DataFrame(rows)

    return np.nan, None, None, None, pd.DataFrame(rows)


def _evaluate_envelope_fallback(T_fit, RESP_fit, cfg: PipelineConfig):
    T_fit = np.asarray(T_fit, dtype=float)
    RESP_fit = np.asarray(RESP_fit, dtype=float)

    idx_peak = int(np.nanargmax(RESP_fit))
    env = np.abs(hilbert(RESP_fit))
    peak_env = float(env[idx_peak])

    li = None
    for i in range(idx_peak, -1, -1):
        if env[i] < cfg.envelope_rel * peak_env:
            li = i + 1
            break

    if li is None:
        li = max(0, idx_peak - int(np.ceil(cfg.min_points / 2)))

    ri = None
    for i in range(idx_peak, len(T_fit)):
        if env[i] < cfg.envelope_rel * peak_env:
            ri = i - 1
            break

    if ri is None:
        ri = min(len(T_fit) - 1, idx_peak + int(np.ceil(cfg.min_points / 2)))

    if ri - li + 1 < cfg.min_points:
        half = max(cfg.min_points // 2, 1)
        li = max(0, idx_peak - half)
        ri = min(len(T_fit) - 1, idx_peak + half)

    det = {
        "T_left": float(T_fit[li]),
        "T_right": float(T_fit[ri]),
        "T_peak": float(T_fit[idx_peak]),
        "idx_peak": idx_peak,
        "method": "envelope_fallback_official",
        "n_points_orig": int(ri - li + 1),
    }

    fit_res = fit_sinusoidal_in_window(
        T_fit,
        RESP_fit,
        detect_res=det,
        use_detect_function=False,
        verbose=False,
    )

    fit_result = fit_res["result"]
    valid = _is_physically_valid(fit_result, cfg)

    return det, fit_res, fit_result, valid


def run_pipeline(
    input_path: str | Path,
    out_dir: str | Path,
    cfg: PipelineConfig = PipelineConfig(),
    make_plots: bool = True,
):
    input_path = Path(input_path)

    base_out = Path(out_dir)
    base_out.mkdir(parents=True, exist_ok=True)

    stem = input_path.stem
    safe = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in stem).strip("_")
    out_dir = base_out / safe
    out_dir.mkdir(parents=True, exist_ok=True)

    t, ug, dt = load_record_txt(input_path)

    Tmin = cfg.tmin if cfg.tmin is not None else max(10 * dt, 0.01)
    T_values, RESP = compute_drs(
        ug=ug,
        dt=dt,
        zeta=cfg.zeta,
        Tmin=Tmin,
        Tmax=cfg.tmax,
        dT=cfg.dT,
    )

    pd.DataFrame({"T": T_values, "DRS": RESP}).to_csv(out_dir / "drs_full.csv", index=False)

    Tp, Dmax, Tmin_zoom, Tmax_zoom, T_fit, RESP_fit, zoom_used = _auto_zoom(
        T_values,
        RESP,
        zoom_frac=cfg.zoom_frac,
        zoom_min_points=cfg.zoom_min_points,
        zoom_fallback_frac=cfg.zoom_fallback_frac,
    )

    selected_alpha = float(cfg.alpha_operational)
    det, fit_res, fit_result, valid = _evaluate_alpha(
        T_fit,
        RESP_fit,
        selected_alpha,
        cfg,
    )

    alpha_strategy = "operational_alpha"
    restricted_df = pd.DataFrame()

    if not valid:
        selected_alpha, det, fit_res, fit_result, restricted_df = _restricted_alpha_search(
            T_fit,
            RESP_fit,
            cfg,
        )
        alpha_strategy = "restricted_search"
        valid = fit_result is not None and _is_physically_valid(fit_result, cfg)

    if not valid and cfg.allow_envelope_fallback:
        det, fit_res, fit_result, valid = _evaluate_envelope_fallback(
            T_fit,
            RESP_fit,
            cfg,
        )
        selected_alpha = np.nan
        alpha_strategy = "official_envelope_fallback"

    if not restricted_df.empty:
        restricted_df.to_csv(out_dir / "restricted_alpha_search.csv", index=False)

    final = {
        "input": str(input_path),
        "dt": float(dt),
        "Tp": float(Tp),
        "Dmax": float(Dmax),
        "zoom_Tmin": float(Tmin_zoom),
        "zoom_Tmax": float(Tmax_zoom),
        "zoom_used": float(zoom_used),

        "alpha_operational": float(cfg.alpha_operational),
        "selected_alpha": float(selected_alpha) if np.isfinite(selected_alpha) else np.nan,
        "alpha_strategy": alpha_strategy,

        "alpha_search_min": float(cfg.alpha_search_min),
        "alpha_search_max": float(cfg.alpha_search_max),
        "alpha_search_step": float(cfg.alpha_search_step),

        "min_R2": float(cfg.min_R2),
        "min_points": int(cfg.min_points),
        "min_td": float(cfg.min_td),
        "allow_envelope_fallback": bool(cfg.allow_envelope_fallback),
        "envelope_rel": float(cfg.envelope_rel),

        "valid_solution": bool(valid),

        "T_left": float(det["T_left"]) if det else np.nan,
        "T_right": float(det["T_right"]) if det else np.nan,
        "window_method": det.get("method", "") if det else "",

        "n_points_window": int(fit_result["n_points"]) if fit_result else np.nan,
        "td": float(fit_result["td"]) if fit_result else np.nan,
        "R2": float(fit_result["R2"]) if fit_result else np.nan,
        "cycles": float(fit_result["cycles"]) if fit_result else np.nan,
        "ssr": float(fit_result["ssr"]) if fit_result else np.nan,
        "best_seed": str(fit_result["best_seed"]) if fit_result else "",
    }

    pd.DataFrame([final]).to_csv(out_dir / "results.csv", index=False)

    if make_plots:
        save_drs_full_png(T_values, RESP, out_dir / "drs_full.png", dpi=300)
        save_drs_zoom_png(T_fit, RESP_fit, Tp, Dmax, out_dir / "drs_zoom.png", dpi=300)

        if fit_result is not None and det is not None:
            mask_win = (T_fit >= det["T_left"]) & (T_fit <= det["T_right"])
            T_win = np.asarray(T_fit[mask_win], dtype=float)
            R_win = np.asarray(RESP_fit[mask_win], dtype=float)

            save_final_sine_fit_png(T_win, R_win, fit_result, out_dir / "sine_fit.png", dpi=300)

    return final

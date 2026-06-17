from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .io import load_record_txt
from .spectrum import compute_drs
from .window import detect_window
from .fit import fit_sinusoidal_in_window

from .plots import (
    save_drs_full_png,
    save_drs_window_png,
    save_second_derivative_window_png,
)


@dataclass
class PipelineConfig:
    """
    Main configuration for the SPECduration workflow.

    The displacement response spectrum, S(T), is computed first.
    The second spectral derivative, S''(T), is used only to detect
    the predominant spectral window. The harmonic fitting is then
    performed on the original DRS values inside that window.
    """

    # DRS computation
    zeta: float = 0.05
    tmin: float | None = None
    tmax: float = 15.0
    dT: float = 0.01

    # Second-derivative window detector
    min_points_window: int = 10
    zero_tol_rel: float | None = None
    extend_window_if_needed: bool = True

    # Harmonic fitting acceptance criteria
    min_points_fit: int = 10
    min_R2_fit: float = 0.98

    # Internal plotting controls
    zoom_frac: float = 0.25
    zoom_min_points: int = 15


def _auto_zoom_for_plot(
    T_values,
    RESP,
    T_peak,
    zoom_frac=0.25,
    zoom_min_points=15,
):
    """
    Build a local zoom around the peak period for diagnostic use.

    This helper is used only to store zoom limits and peak information
    in the diagnostic output. It does not affect the detected window
    or the harmonic fitting.
    """
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    if T_values.size == 0:
        raise RuntimeError("No valid data available to build the diagnostic zoom.")

    Tmin_full = float(T_values[0])
    Tmax_full = float(T_values[-1])

    Tmin_zoom = max(Tmin_full, float(T_peak) * (1.0 - zoom_frac))
    Tmax_zoom = min(Tmax_full, float(T_peak) * (1.0 + zoom_frac))

    mask_zoom = (T_values >= Tmin_zoom) & (T_values <= Tmax_zoom)

    if np.sum(mask_zoom) < zoom_min_points:
        half_points = max(zoom_min_points // 2, 1)
        idx_peak = int(np.argmin(np.abs(T_values - T_peak)))

        i0 = max(0, idx_peak - half_points)
        i1 = min(len(T_values) - 1, idx_peak + half_points)

        Tmin_zoom = float(T_values[i0])
        Tmax_zoom = float(T_values[i1])

        mask_zoom = (T_values >= Tmin_zoom) & (T_values <= Tmax_zoom)

    T_zoom = T_values[mask_zoom]
    S_zoom = RESP[mask_zoom]

    idx_peak_zoom = int(np.argmin(np.abs(T_values - T_peak)))
    Dmax = float(RESP[idx_peak_zoom])

    return Tmin_zoom, Tmax_zoom, T_zoom, S_zoom, Dmax


def _safe_float(value, default=np.nan):
    """
    Convert a value to float while protecting the diagnostics from
    missing or non-finite values.
    """
    try:
        value = float(value)
        return value if np.isfinite(value) else default
    except Exception:
        return default


def _safe_int_or_nan(value):
    """
    Convert a value to integer when possible; otherwise return NaN.
    """
    try:
        if value is None:
            return np.nan

        value = float(value)

        if not np.isfinite(value):
            return np.nan

        return int(value)

    except Exception:
        return np.nan


def _is_valid_solution(fit_result, detect_result, cfg: PipelineConfig) -> bool:
    """
    Check whether the detected window and the harmonic fitting satisfy
    the minimum acceptance criteria.
    """
    if fit_result is None:
        return False

    if detect_result is None:
        return False

    if not bool(detect_result.get("valid_window", False)):
        return False

    R2 = fit_result.get("R2", np.nan)
    td = fit_result.get("td", np.nan)
    n_points = fit_result.get("n_points", 0)

    if not np.isfinite(R2):
        return False

    if not np.isfinite(td):
        return False

    if int(n_points) < int(cfg.min_points_fit):
        return False

    if float(R2) < float(cfg.min_R2_fit):
        return False

    return True


def run_pipeline(
    input_path: str | Path,
    out_dir: str | Path,
    cfg: PipelineConfig = PipelineConfig(),
    make_plots: bool = True,
):
    """
    Run the full SPECduration workflow for a single seismic record.

    The workflow computes the DRS, detects the predominant spectral
    window using S''(T), fits the harmonic model inside that window,
    exports CSV results, and optionally generates the main figures.
    """
    input_path = Path(input_path)

    base_out = Path(out_dir)
    base_out.mkdir(parents=True, exist_ok=True)

    stem = input_path.stem
    safe = "".join(
        c if c.isalnum() or c in ("-", "_") else "_"
        for c in stem
    ).strip("_")

    out_record_dir = base_out / safe
    out_record_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Input ground-motion record
    # ------------------------------------------------------------------

    _, ug, dt = load_record_txt(input_path)

    # ------------------------------------------------------------------
    # 2. Displacement response spectrum, S(T)
    # ------------------------------------------------------------------

    Tmin = cfg.tmin if cfg.tmin is not None else max(10.0 * dt, 0.01)

    T_values, RESP = compute_drs(
        ug=ug,
        dt=dt,
        zeta=cfg.zeta,
        Tmin=Tmin,
        Tmax=cfg.tmax,
        dT=cfg.dT,
    )

    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    pd.DataFrame(
        {
            "T": T_values,
            "DRS": RESP,
        }
    ).to_csv(
        out_record_dir / "drs_full.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # 3. Predominant window detection from S''(T)
    # ------------------------------------------------------------------

    detect_result = detect_window(
        T_in=T_values,
        RESP_in=RESP,
        min_points_in_window=cfg.min_points_window,
        zero_tol_rel=cfg.zero_tol_rel,
        extend_if_needed=cfg.extend_window_if_needed,
        plot=False,
    )

    # ------------------------------------------------------------------
    # 4. Harmonic fitting inside the detected window
    # ------------------------------------------------------------------

    fit_output = fit_sinusoidal_in_window(
        T_values,
        RESP,
        detect_res=detect_result,
        use_detect_function=False,
        verbose=False,
    )

    fit_result = fit_output["result"]

    valid_solution = _is_valid_solution(
        fit_result=fit_result,
        detect_result=detect_result,
        cfg=cfg,
    )

    # ------------------------------------------------------------------
    # 5. Local diagnostic zoom around T_peak
    # ------------------------------------------------------------------

    T_peak = float(detect_result["T_peak"])

    (
        zoom_Tmin,
        zoom_Tmax,
        T_zoom,
        S_zoom,
        Dmax,
    ) = _auto_zoom_for_plot(
        T_values=T_values,
        RESP=RESP,
        T_peak=T_peak,
        zoom_frac=cfg.zoom_frac,
        zoom_min_points=cfg.zoom_min_points,
    )

    # ------------------------------------------------------------------
    # 6. Full diagnostic output
    # ------------------------------------------------------------------

    diagnostics = {
        "input": str(input_path),
        "record": safe,
        "dt_record": float(dt),

        "zeta": float(cfg.zeta),
        "Tmin_DRS": float(Tmin),
        "Tmax_DRS": float(cfg.tmax),
        "dT_DRS": float(cfg.dT),

        "T_peak": float(detect_result["T_peak"]),
        "T_left": float(detect_result["T_left"]),
        "T_right": float(detect_result["T_right"]),

        "idx_peak": int(detect_result["idx_peak"]),
        "idx_left": int(detect_result["idx_left"]),
        "idx_right": int(detect_result["idx_right"]),

        "left_zero_index": _safe_int_or_nan(
            detect_result.get("left_zero_index")
        ),
        "right_zero_index": _safe_int_or_nan(
            detect_result.get("right_zero_index")
        ),

        "S2_left": _safe_float(detect_result.get("S2_left", np.nan)),
        "S2_left_outer": _safe_float(detect_result.get("S2_left_outer", np.nan)),
        "S2_right": _safe_float(detect_result.get("S2_right", np.nan)),
        "S2_right_outer": _safe_float(detect_result.get("S2_right_outer", np.nan)),

        "T_left_outer": _safe_float(detect_result.get("T_left_outer", np.nan)),
        "T_right_outer": _safe_float(detect_result.get("T_right_outer", np.nan)),
        "idx_left_outer": _safe_int_or_nan(
            detect_result.get("idx_left_outer")
        ),
        "idx_right_outer": _safe_int_or_nan(
            detect_result.get("idx_right_outer")
        ),

        "n_points_window": int(detect_result["n_points_window"]),
        "n_points_window_raw": int(detect_result["n_points_window_raw"]),
        "min_points_window": int(detect_result["min_points_in_window"]),

        "window_method": str(detect_result["window_method"]),
        "valid_window": bool(detect_result["valid_window"]),
        "window_extended": bool(detect_result["window_extended"]),
        "extension_points_left": int(detect_result["extension_points_left"]),
        "extension_points_right": int(detect_result["extension_points_right"]),

        "zero_tol": (
            float(detect_result["zero_tol"])
            if np.isfinite(detect_result.get("zero_tol", np.nan))
            else np.nan
        ),
        "zero_tol_rel": (
            float(detect_result["zero_tol_rel"])
            if np.isfinite(detect_result.get("zero_tol_rel", np.nan))
            else np.nan
        ),

        "max_abs_d2R": (
            float(detect_result["max_abs_d2R"])
            if np.isfinite(detect_result.get("max_abs_d2R", np.nan))
            else np.nan
        ),
        "max_abs_d2S": (
            float(detect_result["max_abs_d2S"])
            if np.isfinite(detect_result.get("max_abs_d2S", np.nan))
            else np.nan
        ),

        "regridded": bool(detect_result["regridded"]),

        "min_points_fit": int(cfg.min_points_fit),
        "min_R2_fit": float(cfg.min_R2_fit),
        "valid_solution": bool(valid_solution),

        "td": float(fit_result["td"]) if fit_result else np.nan,
        "R2": float(fit_result["R2"]) if fit_result else np.nan,
        "cycles": float(fit_result["cycles"]) if fit_result else np.nan,
        "ssr": float(fit_result["ssr"]) if fit_result else np.nan,
        "best_seed": str(fit_result["best_seed"]) if fit_result else "",
        "fit_n_points": int(fit_result["n_points"]) if fit_result else np.nan,

        "seed_strategy": str(fit_result.get("seed_strategy", "")) if fit_result else "",
        "n_initial_guesses": int(fit_result.get("n_initial_guesses", 0)) if fit_result else 0,
        "n_td_candidates": int(fit_result.get("n_td_candidates", 0)) if fit_result else 0,
        "n_phi_candidates": int(fit_result.get("n_phi_candidates", 0)) if fit_result else 0,
        "A0_rule": str(fit_result.get("A0_rule", "")) if fit_result else "",
        "td0_rule": str(fit_result.get("td0_rule", "")) if fit_result else "",
        "phi0_rule": str(fit_result.get("phi0_rule", "")) if fit_result else "",
        "td0_candidates": str(fit_result.get("td0_candidates", "")) if fit_result else "",
        "phi0_candidates": str(fit_result.get("phi0_candidates", "")) if fit_result else "",

        "zoom_Tmin": float(zoom_Tmin),
        "zoom_Tmax": float(zoom_Tmax),
        "Dmax": float(Dmax),
    }

    # ------------------------------------------------------------------
    # 7. Compact result table
    # ------------------------------------------------------------------

    results_clean = {
        "record": diagnostics["record"],
        "dt_record": diagnostics["dt_record"],
        "zeta": diagnostics["zeta"],
        "dT_DRS": diagnostics["dT_DRS"],

        "T_peak": diagnostics["T_peak"],
        "T_left": diagnostics["T_left"],
        "T_right": diagnostics["T_right"],

        "n_points_window": diagnostics["n_points_window"],
        "window_method": diagnostics["window_method"],
        "window_extended": diagnostics["window_extended"],
        "valid_window": diagnostics["valid_window"],

        "td": diagnostics["td"],
        "R2": diagnostics["R2"],
        "valid_solution": diagnostics["valid_solution"],
    }

    pd.DataFrame([results_clean]).to_csv(
        out_record_dir / "results.csv",
        index=False,
    )

    pd.DataFrame([diagnostics]).to_csv(
        out_record_dir / "results_diagnostics.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # 8. Spectral derivative data
    # ------------------------------------------------------------------

    pd.DataFrame(
        {
            "T": detect_result["T_u"],
            "DRS": detect_result["R_u"],

            "dS_dT": detect_result.get("dS", detect_result["dR"]),
            "d2S_dT2": detect_result.get("d2S", detect_result["d2R"]),

            "dR_dT": detect_result["dR"],
            "d2R_dT2": detect_result["d2R"],
        }
    ).to_csv(
        out_record_dir / "second_derivative.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # 9. Main output figures
    # ------------------------------------------------------------------

    if make_plots:
        save_drs_full_png(
            T_values,
            RESP,
            out_record_dir / "drs_full.png",
            dpi=300,
        )

        save_drs_window_png(
            T_values=T_values,
            RESP=RESP,
            detect_result=detect_result,
            fit_result=fit_result,
            outpath=out_record_dir / "drs_window.png",
            dpi=300,
        )

        save_second_derivative_window_png(
            detect_result=detect_result,
            outpath=out_record_dir / "second_derivative_window.png",
            dpi=300,
        )

    return results_clean
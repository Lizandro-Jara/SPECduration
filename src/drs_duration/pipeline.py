from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .io import load_record_txt
from .spectrum import compute_drs
from .window import detect_window_by_derivative
from .fit import fit_sinusoidal_in_window
from .alpha_scan import alpha_scan

# plots PRO (los 4)
from .plots import (
    save_drs_full_png,
    save_drs_zoom_png,
    save_r2_vs_alpha_png,
    save_final_sine_fit_png,
)


@dataclass
class PipelineConfig:
    zeta: float = 0.05
    tmin: float | None = None     # si None: usa max(10*dt, 0.01)
    tmax: float = 15.0
    dT: float = 0.01

    zoom_frac: float = 0.20
    zoom_min_points: int = 15
    zoom_fallback_frac: float = 0.35

    alpha_min: float = 0.01
    alpha_max: float = 0.60
    alpha_step: float = 0.01


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


def run_pipeline(
    input_path: str | Path,
    out_dir: str | Path,
    cfg: PipelineConfig = PipelineConfig(),
    make_plots: bool = True,
):
    input_path = Path(input_path)

    # =========================================================
    # NUEVO: subcarpeta automática por cada registro
    # outputs/<nombre_archivo_sin_extension>/
    # =========================================================
    base_out = Path(out_dir)
    base_out.mkdir(parents=True, exist_ok=True)

    stem = input_path.stem  # nombre sin .txt
    safe = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in stem).strip("_")
    out_dir = base_out / safe
    out_dir.mkdir(parents=True, exist_ok=True)
    # =========================================================

    # 1) cargar
    t, ug, dt = load_record_txt(input_path)

    # 2) DRS completo
    Tmin = cfg.tmin if cfg.tmin is not None else max(10 * dt, 0.01)
    T_values, RESP = compute_drs(ug=ug, dt=dt, zeta=cfg.zeta, Tmin=Tmin, Tmax=cfg.tmax, dT=cfg.dT)

    # Guardar tabla DRS
    pd.DataFrame({"T": T_values, "DRS": RESP}).to_csv(out_dir / "drs_full.csv", index=False)

    # 3) zoom
    Tp, Dmax, Tmin_zoom, Tmax_zoom, T_fit, RESP_fit, zoom_used = _auto_zoom(
        T_values, RESP,
        zoom_frac=cfg.zoom_frac,
        zoom_min_points=cfg.zoom_min_points,
        zoom_fallback_frac=cfg.zoom_fallback_frac
    )

    # 4) alpha scan
    df_alpha, best_alpha = alpha_scan(T_fit, RESP_fit, cfg.alpha_min, cfg.alpha_max, cfg.alpha_step, verbose=False)
    df_alpha.to_csv(out_dir / "alpha_scan.csv", index=False)

    # 5) evaluar con alpha óptimo
    det = detect_window_by_derivative(T_fit, RESP_fit, alpha=best_alpha, plot=False)
    fit_res = fit_sinusoidal_in_window(T_fit, RESP_fit, detect_res=det, use_detect_function=False, verbose=False)

    noC = fit_res["result_noC"]
    conC = fit_res["result_C"]

    # 6) guardar resultados finales
    final = {
        "input": str(input_path),
        "dt": float(dt),
        "Tp": Tp,
        "Dmax": Dmax,
        "zoom_Tmin": float(Tmin_zoom),
        "zoom_Tmax": float(Tmax_zoom),
        "zoom_used": float(zoom_used),
        "best_alpha": float(best_alpha),
        "T_left": float(det["T_left"]),
        "T_right": float(det["T_right"]),
        "td_noC": (float(noC["td"]) if noC else np.nan),
        "R2_noC": (float(noC["R2"]) if noC else np.nan),
        "td_C": (float(conC["td"]) if conC else np.nan),
        "R2_C": (float(conC["R2"]) if conC else np.nan),
    }
    pd.DataFrame([final]).to_csv(out_dir / "results.csv", index=False)

    # 7) figuras (LOS 4 PRO)
    if make_plots:
        # fig 1: drs_full.png
        save_drs_full_png(T_values, RESP, out_dir / "drs_full.png", dpi=300)

        # fig 2: drs_zoom.png
        save_drs_zoom_png(T_fit, RESP_fit, Tp, Dmax, out_dir / "drs_zoom.png", dpi=300)

        # fig 3: r2_vs_alpha.png
        save_r2_vs_alpha_png(df_alpha, best_alpha, out_dir / "r2_vs_alpha.png", dpi=300)

        # fig 4: sine_fit.png
        mask_win = (T_fit >= det["T_left"]) & (T_fit <= det["T_right"])
        T_win = np.asarray(T_fit[mask_win], dtype=float)
        R_win = np.asarray(RESP_fit[mask_win], dtype=float)

        save_final_sine_fit_png(T_win, R_win, noC, conC, out_dir / "sine_fit.png", dpi=300)

    return final

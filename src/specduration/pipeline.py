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

    idx_peak_zoom = int(np.argmin(np.abs(T_values - T_peak)))
    Dmax = float(RESP[idx_peak_zoom])

    return Tmin_zoom, Tmax_zoom, Dmax


def _safe_float(value, default=np.nan):
    try:
        value = float(value)
        return value if np.isfinite(value) else default
    except Exception:
        return default


def _safe_int_or_nan(value):
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


def _format_excel_sheet(writer, sheet_name):
    worksheet = writer.sheets[sheet_name]
    worksheet.freeze_panes = "A2"

    for column_cells in worksheet.columns:
        max_length = 0
        column_letter = column_cells[0].column_letter

        for cell in column_cells:
            value = "" if cell.value is None else str(cell.value)
            max_length = max(max_length, len(value))

        worksheet.column_dimensions[column_letter].width = min(max_length + 2, 45)


def _write_excel_report(
    outpath,
    summary_df,
    drs_df,
    derivative_df,
    diagnostics_df,
):
    try:
        with pd.ExcelWriter(outpath, engine="openpyxl") as writer:
            summary_df.to_excel(writer, sheet_name="Summary", index=False)
            drs_df.to_excel(writer, sheet_name="DRS", index=False)
            derivative_df.to_excel(writer, sheet_name="SecondDerivative", index=False)
            diagnostics_df.to_excel(writer, sheet_name="Diagnostics", index=False)

            for sheet_name in writer.sheets:
                _format_excel_sheet(writer, sheet_name)

    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "openpyxl is required to export results.xlsx. "
            "Install it with: pip install openpyxl"
        ) from exc


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

    drs_df = pd.DataFrame(
        {
            "Period T [s]": T_values,
            "DRS S(T) [m]": RESP,
        }
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

    zoom_Tmin, zoom_Tmax, Dmax = _auto_zoom_for_plot(
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
        "Input file": str(input_path),
        "Record": safe,
        "Record time step [s]": float(dt),

        "Damping ratio": float(cfg.zeta),
        "Minimum DRS period [s]": float(Tmin),
        "Maximum DRS period [s]": float(cfg.tmax),
        "DRS period step [s]": float(cfg.dT),

        "Peak period T_peak [s]": float(detect_result["T_peak"]),
        "Left window limit T_left [s]": float(detect_result["T_left"]),
        "Right window limit T_right [s]": float(detect_result["T_right"]),

        "Peak index": int(detect_result["idx_peak"]),
        "Left window index": int(detect_result["idx_left"]),
        "Right window index": int(detect_result["idx_right"]),

        "Left zero-crossing index": _safe_int_or_nan(
            detect_result.get("left_zero_index")
        ),
        "Right zero-crossing index": _safe_int_or_nan(
            detect_result.get("right_zero_index")
        ),

        "S''(T_left)": _safe_float(detect_result.get("S2_left", np.nan)),
        "S'' left outer point": _safe_float(detect_result.get("S2_left_outer", np.nan)),
        "S''(T_right)": _safe_float(detect_result.get("S2_right", np.nan)),
        "S'' right outer point": _safe_float(detect_result.get("S2_right_outer", np.nan)),

        "Left outer period [s]": _safe_float(detect_result.get("T_left_outer", np.nan)),
        "Right outer period [s]": _safe_float(detect_result.get("T_right_outer", np.nan)),
        "Left outer index": _safe_int_or_nan(
            detect_result.get("idx_left_outer")
        ),
        "Right outer index": _safe_int_or_nan(
            detect_result.get("idx_right_outer")
        ),

        "Window points": int(detect_result["n_points_window"]),
        "Raw window points": int(detect_result["n_points_window_raw"]),
        "Minimum window points": int(detect_result["min_points_in_window"]),

        "Window detection method": str(detect_result["window_method"]),
        "Valid window": bool(detect_result["valid_window"]),
        "Window extended": bool(detect_result["window_extended"]),
        "Left extension points": int(detect_result["extension_points_left"]),
        "Right extension points": int(detect_result["extension_points_right"]),

        "Maximum |S''(T)|": (
            float(detect_result["max_abs_d2S"])
            if np.isfinite(detect_result.get("max_abs_d2S", np.nan))
            else np.nan
        ),

        "Regridded spectrum": bool(detect_result["regridded"]),

        "Minimum fit points": int(cfg.min_points_fit),
        "Minimum fit R²": float(cfg.min_R2_fit),
        "Valid solution": bool(valid_solution),

        "Spectral duration td [s]": float(fit_result["td"]) if fit_result else np.nan,
        "R²": float(fit_result["R2"]) if fit_result else np.nan,
        "Cycles inside window": float(fit_result["cycles"]) if fit_result else np.nan,
        "Sum of squared residuals": float(fit_result["ssr"]) if fit_result else np.nan,
        "Best initial seed": str(fit_result["best_seed"]) if fit_result else "",
        "Fit points": int(fit_result["n_points"]) if fit_result else np.nan,

        "Seed strategy": str(fit_result.get("seed_strategy", "")) if fit_result else "",
        "Initial guesses": int(fit_result.get("n_initial_guesses", 0)) if fit_result else 0,
        "td initial candidates": int(fit_result.get("n_td_candidates", 0)) if fit_result else 0,
        "phi initial candidates": int(fit_result.get("n_phi_candidates", 0)) if fit_result else 0,
        "A0 rule": str(fit_result.get("A0_rule", "")) if fit_result else "",
        "td0 rule": str(fit_result.get("td0_rule", "")) if fit_result else "",
        "phi0 rule": str(fit_result.get("phi0_rule", "")) if fit_result else "",
        "td0 candidates": str(fit_result.get("td0_candidates", "")) if fit_result else "",
        "phi0 candidates": str(fit_result.get("phi0_candidates", "")) if fit_result else "",

        "Zoom minimum period [s]": float(zoom_Tmin),
        "Zoom maximum period [s]": float(zoom_Tmax),
        "Peak DRS value [m]": float(Dmax),
    }

    # ------------------------------------------------------------------
    # 7. Compact result table
    # ------------------------------------------------------------------

    summary = {
        "Record": diagnostics["Record"],
        "Record time step [s]": diagnostics["Record time step [s]"],
        "Damping ratio": diagnostics["Damping ratio"],
        "DRS period step [s]": diagnostics["DRS period step [s]"],

        "Peak period T_peak [s]": diagnostics["Peak period T_peak [s]"],
        "Left window limit T_left [s]": diagnostics["Left window limit T_left [s]"],
        "Right window limit T_right [s]": diagnostics["Right window limit T_right [s]"],

        "Window points": diagnostics["Window points"],
        "Spectral duration td [s]": diagnostics["Spectral duration td [s]"],
        "R²": diagnostics["R²"],
        "Valid solution": diagnostics["Valid solution"],
    }

    summary_df = pd.DataFrame([summary])
    diagnostics_df = pd.DataFrame([diagnostics])

    # ------------------------------------------------------------------
    # 8. Spectral derivative data
    # ------------------------------------------------------------------

    derivative_df = pd.DataFrame(
        {
            "Period T [s]": detect_result["T_u"],
            "DRS S(T) [m]": detect_result["S_u"],
            "First derivative S'(T)": detect_result["dS"],
            "Second derivative S''(T)": detect_result["d2S"],
        }
    )

    # ------------------------------------------------------------------
    # 9. CSV and Excel output files
    # ------------------------------------------------------------------

    summary_df.to_csv(
        out_record_dir / "summary.csv",
        sep=";",
        index=False,
    )

    diagnostics_df.to_csv(
        out_record_dir / "results_diagnostics.csv",
        sep=";",
        index=False,
    )

    drs_df.to_csv(
        out_record_dir / "drs_full.csv",
        sep=";",
        index=False,
    )

    derivative_df.to_csv(
        out_record_dir / "second_derivative.csv",
        sep=";",
        index=False,
    )

    _write_excel_report(
        outpath=out_record_dir / "results.xlsx",
        summary_df=summary_df,
        drs_df=drs_df,
        derivative_df=derivative_df,
        diagnostics_df=diagnostics_df,
    )

    # ------------------------------------------------------------------
    # 10. Main output figures
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

    return {
        "record": safe,
        "dt_record": float(dt),
        "zeta": float(cfg.zeta),
        "dT_DRS": float(cfg.dT),
        "T_peak": float(detect_result["T_peak"]),
        "T_left": float(detect_result["T_left"]),
        "T_right": float(detect_result["T_right"]),
        "n_points_window": int(detect_result["n_points_window"]),
        "td": float(fit_result["td"]) if fit_result else np.nan,
        "R2": float(fit_result["R2"]) if fit_result else np.nan,
        "valid_solution": bool(valid_solution),
    }
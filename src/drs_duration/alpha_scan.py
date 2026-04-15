import numpy as np
import pandas as pd

from .window import detect_window_by_derivative
from .fit import fit_sinusoidal_in_window


def alpha_scan(
    T_fit,
    RESP_fit,
    alpha_min=0.01,
    alpha_max=0.60,
    alpha_step=0.01,
    verbose=False,
):
    alphas = np.arange(alpha_min, alpha_max + 1e-12, alpha_step)
    rows = []

    for a in alphas:
        try:
            det = detect_window_by_derivative(
                T_fit,
                RESP_fit,
                alpha=float(a),
                plot=False,
            )

            res = fit_sinusoidal_in_window(
                T_fit,
                RESP_fit,
                detect_res=det,
                use_detect_function=False,
                verbose=False,
            )

            fit_result = res["result"]

            rows.append({
                "alpha": float(a),
                "T_left": float(det["T_left"]),
                "T_right": float(det["T_right"]),
                "n_points": int(det["n_points_orig"]),
                "window_method": det.get("method", ""),
                "R2": (fit_result["R2"] if fit_result else np.nan),
                "td": (fit_result["td"] if fit_result else np.nan),
                "ssr": (fit_result["ssr"] if fit_result else np.nan),
            })

            if verbose:
                print(
                    f"alpha={a:.3f}  "
                    f"R2={rows[-1]['R2']:.6f}  "
                    f"td={rows[-1]['td']:.6f}  "
                    f"n={rows[-1]['n_points']}"
                )

        except Exception:
            rows.append({
                "alpha": float(a),
                "T_left": np.nan,
                "T_right": np.nan,
                "n_points": np.nan,
                "window_method": "",
                "R2": np.nan,
                "td": np.nan,
                "ssr": np.nan,
            })

    df = pd.DataFrame(rows)

    # seleccionar alpha óptimo por máximo R2
    if df["R2"].notna().any():
        idx_best = df["R2"].idxmax()
        best_alpha = float(df.loc[idx_best, "alpha"])
    else:
        best_alpha = np.nan

    return df, best_alpha





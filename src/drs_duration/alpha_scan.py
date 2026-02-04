import numpy as np
import pandas as pd

from .window import detect_window_by_derivative
from .fit import fit_sinusoidal_in_window


def alpha_scan(T_fit, RESP_fit, alpha_min=0.01, alpha_max=0.60, alpha_step=0.01, verbose=False):
    alphas = np.arange(alpha_min, alpha_max + 1e-12, alpha_step)
    rows = []

    for a in alphas:
        try:
            det = detect_window_by_derivative(T_fit, RESP_fit, alpha=float(a), plot=False)
            res = fit_sinusoidal_in_window(
                T_fit, RESP_fit,
                detect_res=det,
                use_detect_function=False,
                verbose=False
            )

            noC = res["result_noC"]
            conC = res["result_C"]

            rows.append({
                "alpha": float(a),
                "T_left": det["T_left"],
                "T_right": det["T_right"],
                "n_points": det["n_points_orig"],
                "R2_noC": (noC["R2"] if noC else np.nan),
                "td_noC": (noC["td"] if noC else np.nan),
                "R2_C": (conC["R2"] if conC else np.nan),
                "td_C": (conC["td"] if conC else np.nan),
            })

            if verbose:
                print(f"alpha={a:.3f}  R2_noC={rows[-1]['R2_noC']:.4f}  td_noC={rows[-1]['td_noC']:.4f}")

        except Exception:
            rows.append({
                "alpha": float(a),
                "T_left": np.nan,
                "T_right": np.nan,
                "n_points": np.nan,
                "R2_noC": np.nan,
                "td_noC": np.nan,
                "R2_C": np.nan,
                "td_C": np.nan,
            })

    df = pd.DataFrame(rows)

    # seleccionar alpha óptimo por máximo R2_noC
    if df["R2_noC"].notna().any():
        idx_best = df["R2_noC"].idxmax()
        best_alpha = float(df.loc[idx_best, "alpha"])
    else:
        best_alpha = np.nan

    return df, best_alpha

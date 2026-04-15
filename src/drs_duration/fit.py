import math
import numpy as np
from scipy.optimize import curve_fit

from .window import detect_window_by_derivative


def sine_fit(T, A, td, phi):
    T = np.asarray(T, dtype=float)
    safe_T = np.where(T == 0, np.finfo(float).eps, T)
    return A * np.abs(np.sin(np.pi * td / safe_T + phi))


def fit_sinusoidal_in_window(
    T_in,
    RESP_in,
    detect_res=None,
    use_detect_function=True,
    alpha_for_detect=0.05,
    verbose=False,
    td_candidates_override=None,
    phi_candidates_override=None,
):
    if detect_res is None and use_detect_function:
        detect_res = detect_window_by_derivative(
            T_in,
            RESP_in,
            alpha=alpha_for_detect,
            plot=False,
        )
    if detect_res is None:
        raise RuntimeError("No se proporcionó detect_res y use_detect_function es False")

    T_left = float(detect_res["T_left"])
    T_right = float(detect_res["T_right"])

    T_arr = np.asarray(T_in, dtype=float)
    R_arr = np.asarray(RESP_in, dtype=float)

    mask = (T_arr >= T_left) & (T_arr <= T_right)
    T_win = T_arr[mask]
    R_win = R_arr[mask]

    # fallback usando rejilla uniforme y suavizado si pocos puntos
    if T_win.size < 4 and "T_u" in detect_res and "R_s" in detect_res:
        T_u = np.asarray(detect_res["T_u"], dtype=float)
        R_s = np.asarray(detect_res["R_s"], dtype=float)
        mask_u = (T_u >= T_left) & (T_u <= T_right)
        T_win = T_u[mask_u]
        R_win = R_s[mask_u]

    if T_win.size < 4:
        raise RuntimeError(f"Pocos puntos en ventana para ajustar: {T_win.size}")

    # normalización
    scale_p = float(np.nanmax(np.abs(R_win)))
    if not np.isfinite(scale_p) or scale_p == 0:
        scale_p = float(np.nanstd(R_win))
    if not np.isfinite(scale_p) or scale_p == 0:
        scale_p = 1.0

    Rn = R_win / scale_p
    max_abs_rwin = float(np.nanmax(np.abs(R_win)))
    A0pn = max_abs_rwin / scale_p if max_abs_rwin > 0 else 1.0

    td_candidates_default = [
        float(np.mean(T_win)),
        max(0.5, float(np.mean(T_win)) * 0.5),
        min(120.0, float(np.mean(T_win)) * 2.0),
        1.0,
    ]
    phi_candidates_default = [0.0, math.pi / 4, -math.pi / 4, math.pi / 2]

    td_candidates = td_candidates_override if td_candidates_override is not None else td_candidates_default
    phi_candidates = phi_candidates_override if phi_candidates_override is not None else phi_candidates_default

    lower = [0.0, 0.001, -math.pi]
    upper = [np.inf, 120.0, math.pi]

    best_ssr = np.inf
    best_popt = None
    best_seed = None

    for td0 in td_candidates:
        for phi0 in phi_candidates:
            p0 = [A0pn, td0, phi0]
            try:
                popt_n, _ = curve_fit(
                    sine_fit,
                    T_win,
                    Rn,
                    p0=p0,
                    bounds=(lower, upper),
                    maxfev=200000,
                )
                resid = Rn - sine_fit(T_win, *popt_n)
                ssr = float(np.sum(resid ** 2))
                if ssr < best_ssr:
                    best_ssr = ssr
                    best_popt = popt_n.copy()
                    best_seed = {"td0": td0, "phi0": phi0}
            except Exception:
                continue

    if best_popt is None:
        result = None
    else:
        popt = best_popt.copy()
        popt[0] *= scale_p  # reescalar A

        R_pred = sine_fit(T_win, *popt)
        ss_tot = float(np.sum((R_win - np.mean(R_win)) ** 2))
        ss_res = float(np.sum((R_win - R_pred) ** 2))
        R2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

        td_val = float(popt[1])
        deltaT = float(T_right - T_left)
        cycles = float(deltaT / td_val) if td_val != 0 else np.nan

        result = dict(
            popt=popt,
            R2=R2,
            td=td_val,
            cycles=cycles,
            best_seed=best_seed,
            ssr=best_ssr,
            n_points=int(len(T_win)),
            T_left=T_left,
            T_right=T_right,
            window_method=detect_res.get("method", ""),
        )

    if verbose:
        print(f"Ventana: {T_left:.6f}–{T_right:.6f} (n={len(T_win)})")
        if result:
            print(
                f"R2={result['R2']:.6f}, td={result['td']:.6f}, "
                f"seed={result['best_seed']}, ssr={result['ssr']:.6e}"
            )

    return dict(
        result=result,
        detect_res=detect_res,
        T_win=T_win,
        R_win=R_win,
    )


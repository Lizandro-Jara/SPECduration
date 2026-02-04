import math
import numpy as np
from scipy.optimize import curve_fit

from .window import detect_window_by_derivative


def sine_fit(T, A, td, phi):
    T = np.asarray(T, dtype=float)
    safe_T = np.where(T == 0, np.finfo(float).eps, T)
    return A * np.abs(np.sin(np.pi * td / safe_T + phi))


def sine_fit_C(T, A, td, phi, C):
    T = np.asarray(T, dtype=float)
    safe_T = np.where(T == 0, np.finfo(float).eps, T)
    return A * np.abs(np.sin(np.pi * td / safe_T + phi)) + C


def fit_sinusoidal_in_window(
    T_in,
    RESP_in,
    detect_res=None,
    use_detect_function=True,
    alpha_for_detect=0.05,
    verbose=False,
):
    if detect_res is None and use_detect_function:
        detect_res = detect_window_by_derivative(T_in, RESP_in, alpha=alpha_for_detect, plot=False)
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
    A0pn = float(np.nanmax(np.abs(R_win)) / scale_p) if np.nanmax(np.abs(R_win)) > 0 else 1.0

    td_candidates = [
        float(np.mean(T_win)),
        max(0.5, float(np.mean(T_win)) * 0.5),
        min(120.0, float(np.mean(T_win)) * 2.0),
        1.0,
    ]
    phi_candidates = [0.0, math.pi / 4, -math.pi / 4, math.pi / 2]

    lower = [0.0, 0.001, -math.pi]
    upper = [np.inf, 120.0, math.pi]

    lowerC = [0.0, 0.001, -math.pi, -0.5]
    upperC = [np.inf, 120.0, math.pi, 0.5]

    # --- sin C ---
    best_ssr = np.inf
    best_popt = None
    for td0 in td_candidates:
        for phi0 in phi_candidates:
            p0 = [A0pn, td0, phi0]
            try:
                popt_n, _ = curve_fit(
                    sine_fit, T_win, Rn,
                    p0=p0, bounds=(lower, upper),
                    maxfev=200000
                )
                resid = Rn - sine_fit(T_win, *popt_n)
                ssr = float(np.sum(resid ** 2))
                if ssr < best_ssr:
                    best_ssr = ssr
                    best_popt = popt_n.copy()
            except Exception:
                continue

    if best_popt is None:
        result_noC = None
    else:
        popt_noC = best_popt.copy()
        popt_noC[0] *= scale_p
        R_pred = sine_fit(T_win, *popt_noC)
        ss_tot = float(np.sum((R_win - np.mean(R_win)) ** 2))
        ss_res = float(np.sum((R_win - R_pred) ** 2))
        R2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan
        td_val = float(popt_noC[1])
        deltaT = float(T_right - T_left)
        cycles = float(deltaT / td_val) if td_val != 0 else np.nan
        result_noC = dict(popt=popt_noC, R2=R2, td=td_val, cycles=cycles)

    # --- con C ---
    best_ssrC = np.inf
    best_poptC = None
    for td0 in td_candidates:
        for phi0 in phi_candidates:
            p0c = [A0pn, td0, phi0, 0.0]
            try:
                popt_nC, _ = curve_fit(
                    sine_fit_C, T_win, Rn,
                    p0=p0c, bounds=(lowerC, upperC),
                    maxfev=200000
                )
                residC = Rn - sine_fit_C(T_win, *popt_nC)
                ssrC = float(np.sum(residC ** 2))
                if ssrC < best_ssrC:
                    best_ssrC = ssrC
                    best_poptC = popt_nC.copy()
            except Exception:
                continue

    if best_poptC is None:
        result_C = None
    else:
        popt_C = best_poptC.copy()
        popt_C[0] *= scale_p
        popt_C[3] *= scale_p
        R_pred_C = sine_fit_C(T_win, *popt_C)
        ss_totC = float(np.sum((R_win - np.mean(R_win)) ** 2))
        ss_resC = float(np.sum((R_win - R_pred_C) ** 2))
        R2C = 1 - ss_resC / ss_totC if ss_totC != 0 else np.nan
        td_valC = float(popt_C[1])
        deltaT = float(T_right - T_left)
        cyclesC = float(deltaT / td_valC) if td_valC != 0 else np.nan
        result_C = dict(popt=popt_C, R2=R2C, td=td_valC, cycles=cyclesC)

    if verbose:
        print(f"Ventana: {T_left:.6f}–{T_right:.6f} (n={len(T_win)})")
        if result_noC:
            print(f"Sin C: R2={result_noC['R2']:.6f}, td={result_noC['td']:.6f}")
        if result_C:
            print(f"Con C: R2={result_C['R2']:.6f}, td={result_C['td']:.6f}")

    return dict(result_noC=result_noC, result_C=result_C, detect_res=detect_res, T_win=T_win, R_win=R_win)

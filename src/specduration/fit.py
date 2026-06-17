
import math
import numpy as np
from scipy.optimize import curve_fit

from .window import detect_window


# ----------------------------------------------------------------------
# Harmonic spectral model
# ----------------------------------------------------------------------

def spectral_sine_model(T, A, td, phi):
    """
    Harmonic model fitted to the displacement response spectrum, S(T),
    inside the detected spectral window.

        S(T) = A * |sin(pi * td / T + phi)|

    Parameters
    ----------
    T : array-like
        Period values.
    A : float
        Amplitude parameter.
    td : float
        Spectral duration parameter.
    phi : float
        Phase parameter.

    Returns
    -------
    ndarray
        Modeled DRS values.
    """
    T = np.asarray(T, dtype=float)
    safe_T = np.where(T == 0, np.finfo(float).eps, T)

    return A * np.abs(np.sin(np.pi * td / safe_T + phi))


def sine_fit(T, A, td, phi):
    """
    Backward-compatible alias for spectral_sine_model.
    """
    return spectral_sine_model(T, A, td, phi)


# ----------------------------------------------------------------------
# Deterministic initial guesses
# ----------------------------------------------------------------------

def _build_initial_guesses(T_win, S_scaled, detect_res):
    """
    Build deterministic initial guesses for the nonlinear harmonic fit.

    The initial guesses do not define the window and do not use S'(T) or
    S''(T). They are used only to improve numerical stability and reduce
    sensitivity to local minima. The same generation rules are applied to
    all records.
    """
    T_win = np.asarray(T_win, dtype=float)
    S_scaled = np.asarray(S_scaled, dtype=float)

    T_left = float(np.min(T_win))
    T_right = float(np.max(T_win))
    T_mean = float(np.mean(T_win))
    T_peak = float(detect_res.get("T_peak", T_mean))

    window_width = max(float(T_right - T_left), 1e-6)

    A0 = float(np.nanmax(np.abs(S_scaled)))
    if not np.isfinite(A0) or A0 <= 0:
        A0 = 1.0

    td_candidates = []

    # Candidate durations based on the spectral peak period.
    for factor in [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0]:
        td_candidates.append(factor * T_peak)

    # Candidate durations based on the mean period inside the window.
    for factor in [1.0, 1.5, 2.0, 2.5]:
        td_candidates.append(factor * T_mean)

    # Candidate durations based on the detected window width.
    for factor in [1.0, 2.0, 4.0, 8.0, 12.0]:
        td_candidates.append(factor * window_width)

    # Absolute duration candidates for short and long records.
    td_candidates.extend(
        [
            1.0,
            3.0,
            5.0,
            7.5,
            10.0,
            15.0,
            20.0,
            25.0,
            30.0,
            40.0,
        ]
    )

    td_candidates = np.asarray(td_candidates, dtype=float)
    td_candidates = td_candidates[np.isfinite(td_candidates)]
    td_candidates = td_candidates[
        (td_candidates >= 0.001) & (td_candidates <= 120.0)
    ]
    td_candidates = np.unique(np.round(td_candidates, 6))

    phi_candidates = np.asarray(
        [
            -math.pi,
            -3.0 * math.pi / 4.0,
            -math.pi / 2.0,
            -math.pi / 4.0,
            0.0,
            math.pi / 4.0,
            math.pi / 2.0,
            3.0 * math.pi / 4.0,
            math.pi,
        ],
        dtype=float,
    )

    guesses = []

    for td0 in td_candidates:
        for phi0 in phi_candidates:
            guesses.append([A0, float(td0), float(phi0)])

    seed_info = {
        "A0_rule": "A0 = max(|S_scaled|) inside the detected window",
        "td0_rule": (
            "Deterministic candidates based on T_peak, T_mean, "
            "window_width, and fixed absolute values"
        ),
        "phi0_rule": "Fixed phase candidates distributed from -pi to pi",
        "n_td_candidates": int(len(td_candidates)),
        "n_phi_candidates": int(len(phi_candidates)),
        "n_initial_guesses": int(len(guesses)),
        "td0_candidates": td_candidates.tolist(),
        "phi0_candidates": phi_candidates.tolist(),
    }

    return guesses, seed_info


# ----------------------------------------------------------------------
# Harmonic fitting inside the detected window
# ----------------------------------------------------------------------

def fit_sinusoidal_in_window(
    T_in,
    RESP_in,
    detect_res=None,
    use_detect_function=True,
    verbose=False,
):
    """
    Fit the harmonic spectral model to the original DRS values inside the
    detected spectral window.

    The detected window is defined before the fitting stage using S''(T).
    This function only performs the nonlinear harmonic fit on S(T) between
    T_left and T_right.
    """

    # ------------------------------------------------------------------
    # 1. Get or compute the detected window
    # ------------------------------------------------------------------

    if detect_res is None and use_detect_function:
        detect_res = detect_window(
            T_in=T_in,
            RESP_in=RESP_in,
            plot=False,
        )

    if detect_res is None:
        raise RuntimeError(
            "detect_res was not provided and use_detect_function is False."
        )

    T_left = float(detect_res["T_left"])
    T_right = float(detect_res["T_right"])

    # ------------------------------------------------------------------
    # 2. Extract original S(T) values inside the window
    # ------------------------------------------------------------------

    T = np.asarray(T_in, dtype=float)
    S = np.asarray(RESP_in, dtype=float)

    mask = (
        np.isfinite(T)
        & np.isfinite(S)
        & (T >= T_left)
        & (T <= T_right)
    )

    T_win = T[mask]
    S_win = S[mask]

    if T_win.size < 4:
        raise RuntimeError(
            f"Not enough points inside the detected window for fitting: {T_win.size}"
        )

    order = np.argsort(T_win)
    T_win = T_win[order]
    S_win = S_win[order]

    # ------------------------------------------------------------------
    # 3. Numerical scaling
    # ------------------------------------------------------------------

    scale = float(np.nanmax(np.abs(S_win)))

    if not np.isfinite(scale) or scale == 0:
        scale = 1.0

    S_scaled = S_win / scale

    # ------------------------------------------------------------------
    # 4. Deterministic initial guesses
    # ------------------------------------------------------------------

    initial_guesses, seed_info = _build_initial_guesses(
        T_win=T_win,
        S_scaled=S_scaled,
        detect_res=detect_res,
    )

    lower_bounds = [0.0, 0.001, -math.pi]
    upper_bounds = [np.inf, 120.0, math.pi]

    best_popt = None
    best_ssr = np.inf
    best_seed = None

    # ------------------------------------------------------------------
    # 5. Multi-start nonlinear fitting
    # ------------------------------------------------------------------

    for p0 in initial_guesses:
        try:
            popt_scaled, _ = curve_fit(
                spectral_sine_model,
                T_win,
                S_scaled,
                p0=p0,
                bounds=(lower_bounds, upper_bounds),
                maxfev=200000,
            )

            residual = S_scaled - spectral_sine_model(T_win, *popt_scaled)
            ssr = float(np.sum(residual ** 2))

            if ssr < best_ssr:
                best_ssr = ssr
                best_popt = popt_scaled.copy()
                best_seed = {
                    "A0": float(p0[0]),
                    "td0": float(p0[1]),
                    "phi0": float(p0[2]),
                }

        except Exception:
            continue

    if best_popt is None:
        return dict(
            result=None,
            detect_res=detect_res,
            T_win=T_win,
            S_win=S_win,

            # Backward-compatible alias for older scripts.
            R_win=S_win,
        )

    # ------------------------------------------------------------------
    # 6. Rescale amplitude and compute fit metrics
    # ------------------------------------------------------------------

    popt = best_popt.copy()
    popt[0] = popt[0] * scale

    S_pred = spectral_sine_model(T_win, *popt)

    ss_tot = float(np.sum((S_win - np.mean(S_win)) ** 2))
    ss_res = float(np.sum((S_win - S_pred) ** 2))

    R2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan

    A = float(popt[0])
    td = float(popt[1])
    phi = float(popt[2])

    window_width = float(T_right - T_left)
    cycles = float(window_width / td) if td != 0 else np.nan

    result = dict(
        popt=popt,
        A=A,
        td=td,
        phi=phi,
        R2=float(R2),
        cycles=cycles,
        ssr=float(best_ssr),

        best_seed=best_seed,

        seed_strategy=(
            "Deterministic initial guesses applied with the same rule "
            "to all records. They do not modify the detected window."
        ),
        n_initial_guesses=int(seed_info["n_initial_guesses"]),
        n_td_candidates=int(seed_info["n_td_candidates"]),
        n_phi_candidates=int(seed_info["n_phi_candidates"]),
        td0_candidates=seed_info["td0_candidates"],
        phi0_candidates=seed_info["phi0_candidates"],
        A0_rule=seed_info["A0_rule"],
        td0_rule=seed_info["td0_rule"],
        phi0_rule=seed_info["phi0_rule"],

        n_points=int(T_win.size),
        T_left=T_left,
        T_right=T_right,
        window_width=window_width,
        window_method=detect_res.get("window_method", detect_res.get("method", "")),
    )

    if verbose:
        print(f"Window: {T_left:.6f}–{T_right:.6f} s")
        print(f"Number of points: {T_win.size}")
        print(f"td = {td:.6f} s")
        print(f"R2 = {R2:.6f}")
        print(f"Cycles = {cycles:.6f}")
        print(f"Best seed = {best_seed}")
        print(f"Initial guesses = {seed_info['n_initial_guesses']}")

    return dict(
        result=result,
        detect_res=detect_res,
        T_win=T_win,
        S_win=S_win,

        # Backward-compatible alias for older scripts.
        R_win=S_win,
    )
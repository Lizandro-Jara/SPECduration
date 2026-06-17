import numpy as np
from scipy.interpolate import interp1d


DEFAULT_MIN_POINTS_IN_WINDOW = 10


# ----------------------------------------------------------------------
# Input spectrum preparation
# ----------------------------------------------------------------------

def preprocess_T_RESP(T_in, RESP_in):
    """
    Clean the input DRS data and sort it by period.

    Parameters
    ----------
    T_in : array-like
        Period values.
    RESP_in : array-like
        Displacement response spectrum values, S(T).

    Returns
    -------
    T : ndarray
        Clean and sorted period vector.
    S : ndarray
        Clean and sorted DRS vector.
    """
    T = np.asarray(T_in, dtype=float)
    S = np.asarray(RESP_in, dtype=float)

    mask = np.isfinite(T) & np.isfinite(S)
    T = T[mask]
    S = S[mask]

    if T.size == 0:
        raise RuntimeError("No valid period values remain after filtering NaN or Inf.")

    idx = np.argsort(T)

    return T[idx], S[idx]


def to_uniform_grid_if_needed(T, S, tol_rel=0.2):
    """
    Interpolate the DRS onto a uniform period grid when the input spacing
    is irregular.

    The second derivative is more stable when the spectrum is evaluated
    on a nearly uniform period grid.
    """
    T = np.asarray(T, dtype=float)
    S = np.asarray(S, dtype=float)

    if len(T) < 3:
        return T, S, False

    dT = np.diff(T)
    median_dT = np.median(dT)

    if median_dT <= 0:
        return T, S, False

    irregularity = np.max(np.abs(dT - median_dT)) / (median_dT + 1e-12)

    if irregularity < tol_rel:
        return T, S, False

    T_new = np.arange(T[0], T[-1] + 0.5 * median_dT, median_dT)

    interpolator = interp1d(
        T,
        S,
        kind="cubic",
        bounds_error=False,
        fill_value="extrapolate",
    )

    S_new = interpolator(T_new)

    return T_new, S_new, True


# ----------------------------------------------------------------------
# Spectral derivatives
# ----------------------------------------------------------------------

def compute_second_derivative(T, S):
    """
    Compute the first and second numerical derivatives of the DRS.

    Notation
    --------
    S(T)   : displacement response spectrum
    S'(T)  : first spectral derivative
    S''(T) : second spectral derivative

    The derivatives are computed directly from the DRS values using
    numerical gradients. No smoothing, envelope, alpha threshold, or
    Savitzky-Golay filter is applied.
    """
    T = np.asarray(T, dtype=float)
    S = np.asarray(S, dtype=float)

    if len(T) < 5:
        raise RuntimeError("At least five DRS points are required to compute S''(T).")

    dS = np.gradient(S, T)
    d2S = np.gradient(dS, T)

    return dS, d2S


def _has_real_sign_crossing(a, b):
    """
    Check whether two consecutive S''(T) values define a real sign crossing.
    """
    if not np.isfinite(a) or not np.isfinite(b):
        return False

    if a == 0.0 or b == 0.0:
        return True

    return a * b < 0.0


# ----------------------------------------------------------------------
# Zero-crossing candidates around the DRS peak
# ----------------------------------------------------------------------

def _left_crossing_candidates(d2S, idx_peak):
    """
    Find sign-crossing candidates to the left of the peak period.

    Candidates are returned from the closest to the peak to the farthest.
    For the left side, the window boundary is the first internal point
    after the external crossing point.

    external_index = i
    internal_index = i + 1
    """
    candidates = []

    for i in range(idx_peak - 1, 0, -1):
        if _has_real_sign_crossing(d2S[i], d2S[i + 1]):
            internal_index = i + 1
            external_index = i
            candidates.append((internal_index, external_index))

    return candidates


def _right_crossing_candidates(d2S, idx_peak):
    """
    Find sign-crossing candidates to the right of the peak period.

    Candidates are returned from the closest to the peak to the farthest.
    For the right side, the window boundary is the last internal point
    before the external crossing point.

    internal_index = i - 1
    external_index = i
    """
    candidates = []
    n = len(d2S)

    for i in range(idx_peak + 1, n - 1):
        if _has_real_sign_crossing(d2S[i - 1], d2S[i]):
            internal_index = i - 1
            external_index = i
            candidates.append((internal_index, external_index))

    return candidates


def _select_first_valid_crossing_pair(
    left_candidates,
    right_candidates,
    idx_peak,
    min_points,
):
    """
    Select the narrowest valid pair of real S''(T) sign crossings around
    the peak period.

    If the closest crossings produce a very small window, the detector
    checks farther crossings until the minimum number of points is met.
    """
    valid_pairs = []

    for i_left, left_cross in left_candidates:
        for i_right, right_cross in right_candidates:

            if i_left >= idx_peak:
                continue

            if i_right <= idx_peak:
                continue

            n_points = i_right - i_left + 1

            if n_points >= min_points:
                valid_pairs.append(
                    {
                        "i_left": int(i_left),
                        "i_right": int(i_right),
                        "left_zero_index": int(left_cross),
                        "right_zero_index": int(right_cross),
                        "n_points": int(n_points),
                        "width_index": int(i_right - i_left),
                    }
                )

    if not valid_pairs:
        return None

    valid_pairs = sorted(valid_pairs, key=lambda x: x["width_index"])

    return valid_pairs[0]


# ----------------------------------------------------------------------
# Fallback window strategies
# ----------------------------------------------------------------------

def _low_curvature_fallback(abs_d2S, idx_peak, n, min_points):
    """
    Fallback strategy used when no useful real sign-crossing pair is found.

    It searches for low-curvature points on both sides of the peak, while
    keeping a minimum distance from the peak to avoid an excessively small
    window.
    """
    min_side_distance = max(min_points // 2, 3)

    left_end = idx_peak - min_side_distance
    right_start = idx_peak + min_side_distance

    if left_end <= 0:
        i_left = max(0, idx_peak - min_side_distance)
        left_zero_index = None
    else:
        left_range = np.arange(0, left_end + 1)
        local = abs_d2S[left_range]
        i_left = int(left_range[np.nanargmin(local)])
        left_zero_index = i_left

    if right_start >= n - 1:
        i_right = min(n - 1, idx_peak + min_side_distance)
        right_zero_index = None
    else:
        right_range = np.arange(right_start, n)
        local = abs_d2S[right_range]
        i_right = int(right_range[np.nanargmin(local)])
        right_zero_index = i_right

    if i_left >= idx_peak:
        i_left = max(0, idx_peak - min_side_distance)

    if i_right <= idx_peak:
        i_right = min(n - 1, idx_peak + min_side_distance)

    return {
        "i_left": int(i_left),
        "i_right": int(i_right),
        "left_zero_index": left_zero_index,
        "right_zero_index": right_zero_index,
    }


def _extend_window_to_min_points(i_left, i_right, idx_peak, n, min_points):
    """
    Final fallback: extend the detected window around the peak until the
    minimum number of points is reached.
    """
    i_left = int(i_left)
    i_right = int(i_right)

    if i_left > idx_peak:
        i_left = idx_peak

    if i_right < idx_peak:
        i_right = idx_peak

    initial_left = i_left
    initial_right = i_right

    while (i_right - i_left + 1) < min_points:
        can_left = i_left > 0
        can_right = i_right < n - 1

        if not can_left and not can_right:
            break

        left_distance = idx_peak - i_left
        right_distance = i_right - idx_peak

        if can_left and can_right:
            if left_distance <= right_distance:
                i_left -= 1
            else:
                i_right += 1
        elif can_left:
            i_left -= 1
        elif can_right:
            i_right += 1

    extension_points_left = initial_left - i_left
    extension_points_right = i_right - initial_right

    return (
        int(i_left),
        int(i_right),
        int(extension_points_left),
        int(extension_points_right),
    )


# ----------------------------------------------------------------------
# Boundary diagnostics
# ----------------------------------------------------------------------

def _extract_S2_boundary_diagnostics(
    T_u,
    d2S,
    i_left,
    i_right,
    left_zero_index,
    right_zero_index,
):
    """
    Extract S''(T) values at the detected window boundaries.

    These diagnostics verify whether T_left and T_right correspond to
    internal points immediately adjacent to a sign change of S''(T).

    Reported values
    ---------------
    S2_left
        S''(T) at T_left.
    S2_left_outer
        S''(T) at the closest external point on the left side.
    S2_right
        S''(T) at T_right.
    S2_right_outer
        S''(T) at the closest external point on the right side.
    """
    n = len(d2S)

    S2_left = float(d2S[i_left]) if 0 <= i_left < n else np.nan
    S2_right = float(d2S[i_right]) if 0 <= i_right < n else np.nan

    if left_zero_index is not None and 0 <= int(left_zero_index) < n:
        i_left_outer = int(left_zero_index)
    else:
        i_left_outer = i_left - 1 if i_left - 1 >= 0 else None

    if right_zero_index is not None and 0 <= int(right_zero_index) < n:
        i_right_outer = int(right_zero_index)
    else:
        i_right_outer = i_right + 1 if i_right + 1 < n else None

    S2_left_outer = (
        float(d2S[i_left_outer])
        if i_left_outer is not None and 0 <= i_left_outer < n
        else np.nan
    )

    S2_right_outer = (
        float(d2S[i_right_outer])
        if i_right_outer is not None and 0 <= i_right_outer < n
        else np.nan
    )

    T_left_outer = (
        float(T_u[i_left_outer])
        if i_left_outer is not None and 0 <= i_left_outer < n
        else np.nan
    )

    T_right_outer = (
        float(T_u[i_right_outer])
        if i_right_outer is not None and 0 <= i_right_outer < n
        else np.nan
    )

    return {
        "S2_left": S2_left,
        "S2_left_outer": S2_left_outer,
        "S2_right": S2_right,
        "S2_right_outer": S2_right_outer,
        "T_left_outer": T_left_outer,
        "T_right_outer": T_right_outer,
        "idx_left_outer": int(i_left_outer) if i_left_outer is not None else None,
        "idx_right_outer": int(i_right_outer) if i_right_outer is not None else None,
    }


# ----------------------------------------------------------------------
# Main second-derivative window detector
# ----------------------------------------------------------------------

def detect_window_by_second_derivative(
    T_in,
    RESP_in,
    min_points_in_window=DEFAULT_MIN_POINTS_IN_WINDOW,
    zero_tol_rel=None,
    extend_if_needed=True,
    plot=False,
):
    """
    Detect the predominant spectral window using S''(T).

    Main strategy
    -------------
    The detector identifies the DRS peak period, T_peak, and searches for
    real sign crossings of S''(T) on both sides of the peak. The selected
    window boundaries are the internal points adjacent to those crossings.

    Fallback strategy
    -----------------
    If no useful sign-crossing pair is found, low-curvature points away
    from the peak are used. If the resulting window still has too few
    points, it is extended around T_peak until the minimum number of
    points is reached.

    Notes
    -----
    S''(T) is used only for window detection. The harmonic fitting is
    performed later on the original DRS values, S(T), inside the selected
    window.

    The parameters zero_tol_rel and plot are kept for API compatibility.
    They are not used by the current real sign-crossing detector.
    """
    _ = zero_tol_rel
    _ = plot

    T0, S0 = preprocess_T_RESP(T_in, RESP_in)

    idx_peak_original = int(np.nanargmax(S0))
    T_peak_original = float(T0[idx_peak_original])

    T_u, S_u, regridded = to_uniform_grid_if_needed(T0, S0)

    T_u = np.asarray(T_u, dtype=float)
    S_u = np.asarray(S_u, dtype=float)

    n = len(T_u)

    if n < 5:
        raise RuntimeError("Too few DRS points are available to detect a window.")

    idx_peak = int(np.argmin(np.abs(T_u - T_peak_original)))
    T_peak = float(T_u[idx_peak])

    dS, d2S = compute_second_derivative(T_u, S_u)
    abs_d2S = np.abs(d2S)

    # ------------------------------------------------------------------
    # 1. Real S''(T) sign crossings around T_peak
    # ------------------------------------------------------------------

    left_candidates = _left_crossing_candidates(d2S, idx_peak)
    right_candidates = _right_crossing_candidates(d2S, idx_peak)

    selected = _select_first_valid_crossing_pair(
        left_candidates=left_candidates,
        right_candidates=right_candidates,
        idx_peak=idx_peak,
        min_points=min_points_in_window,
    )

    window_extended = False
    extension_points_left = 0
    extension_points_right = 0

    if selected is not None:
        i_left = selected["i_left"]
        i_right = selected["i_right"]
        left_zero_index = selected["left_zero_index"]
        right_zero_index = selected["right_zero_index"]
        method = "second_derivative_real_zero_crossing"

    else:
        # --------------------------------------------------------------
        # 2. Low-curvature fallback
        # --------------------------------------------------------------

        fallback = _low_curvature_fallback(
            abs_d2S=abs_d2S,
            idx_peak=idx_peak,
            n=n,
            min_points=min_points_in_window,
        )

        i_left = fallback["i_left"]
        i_right = fallback["i_right"]
        left_zero_index = fallback["left_zero_index"]
        right_zero_index = fallback["right_zero_index"]
        method = "second_derivative_low_curvature_fallback"

    # ------------------------------------------------------------------
    # 3. Minimum-point extension
    # ------------------------------------------------------------------

    n_points_window_raw = int(i_right - i_left + 1)

    if n_points_window_raw < min_points_in_window and extend_if_needed:
        window_extended = True

        (
            i_left,
            i_right,
            extension_points_left,
            extension_points_right,
        ) = _extend_window_to_min_points(
            i_left=i_left,
            i_right=i_right,
            idx_peak=idx_peak,
            n=n,
            min_points=min_points_in_window,
        )

        method = method + "_extended"

    n_points_window = int(i_right - i_left + 1)

    T_left = float(T_u[i_left])
    T_right = float(T_u[i_right])

    valid_window = bool(n_points_window >= min_points_in_window)

    n_points_original = int(
        np.sum((T0 >= T_left) & (T0 <= T_right))
    )

    max_abs_d2S = float(np.nanmax(abs_d2S)) if abs_d2S.size > 0 else np.nan

    # ------------------------------------------------------------------
    # 4. S''(T) values at the window boundaries
    # ------------------------------------------------------------------

    S2_diag = _extract_S2_boundary_diagnostics(
        T_u=T_u,
        d2S=d2S,
        i_left=i_left,
        i_right=i_right,
        left_zero_index=left_zero_index,
        right_zero_index=right_zero_index,
    )

    return dict(
        T_left=T_left,
        T_right=T_right,
        T_peak=T_peak,

        idx_peak=idx_peak,
        idx_left=int(i_left),
        idx_right=int(i_right),
        idx_peak_original=idx_peak_original,

        left_zero_index=left_zero_index,
        right_zero_index=right_zero_index,

        T_u=T_u,

        # Preferred S(T) names.
        S_u=S_u,
        S_s=S_u.copy(),
        dS=dS,
        d2S=d2S,

        # Backward-compatible names kept for older pipeline and plotting code.
        R_u=S_u,
        R_s=S_u.copy(),
        dR=dS,
        d2R=d2S,

        method=method,
        window_method=method,

        regridded=bool(regridded),

        n_points_orig=int(n_points_original),
        n_points_window=int(n_points_window),
        n_points_window_raw=int(n_points_window_raw),
        min_points_in_window=int(min_points_in_window),

        valid_window=bool(valid_window),
        window_extended=bool(window_extended),
        extension_points_left=int(extension_points_left),
        extension_points_right=int(extension_points_right),

        max_abs_d2S=max_abs_d2S,
        max_abs_d2R=max_abs_d2S,

        zero_tol=np.nan,
        zero_tol_rel=np.nan,
        threshold=np.nan,
        mask_good=None,

        S2_left=S2_diag["S2_left"],
        S2_left_outer=S2_diag["S2_left_outer"],
        S2_right=S2_diag["S2_right"],
        S2_right_outer=S2_diag["S2_right_outer"],
        T_left_outer=S2_diag["T_left_outer"],
        T_right_outer=S2_diag["T_right_outer"],
        idx_left_outer=S2_diag["idx_left_outer"],
        idx_right_outer=S2_diag["idx_right_outer"],
    )


def detect_window(
    T_in,
    RESP_in,
    min_points_in_window=DEFAULT_MIN_POINTS_IN_WINDOW,
    zero_tol_rel=None,
    extend_if_needed=True,
    plot=False,
):
    """
    Public window-detection entry point used by the SPECduration pipeline.
    """
    return detect_window_by_second_derivative(
        T_in=T_in,
        RESP_in=RESP_in,
        min_points_in_window=min_points_in_window,
        zero_tol_rel=zero_tol_rel,
        extend_if_needed=extend_if_needed,
        plot=plot,
    )

import numpy as np
from scipy.signal import savgol_filter, hilbert
from scipy.interpolate import interp1d


DEFAULT_sg_frac = 0.12
DEFAULT_sg_poly = 3
DEFAULT_wide_frac = 0.20
DEFAULT_rel_env_fallback = 0.85
DEFAULT_min_points_in_window = 7


def ensure_odd(x: int) -> int:
    x = int(x)
    if x % 2 == 0:
        x += 1
    return max(3, x)


def preprocess_T_RESP(T_in, RESP_in):
    T = np.asarray(T_in, dtype=float)
    R = np.asarray(RESP_in, dtype=float)
    mask = np.isfinite(T) & np.isfinite(R)
    T, R = T[mask], R[mask]
    if T.size == 0:
        raise RuntimeError("T vacío tras filtrar NaN/Inf")
    idx = np.argsort(T)
    return T[idx], R[idx]


def to_uniform_grid_if_needed(T, R, tol_rel=0.2):
    if len(T) < 3:
        return T, R, False

    d = np.diff(T)
    med = np.median(d)
    if med <= 0:
        return T, R, False

    if np.max(np.abs(d - med)) / (med + 1e-12) < tol_rel:
        return T, R, False

    dx = med
    Tnew = np.arange(T[0], T[-1] + dx * 0.5, dx)
    f = interp1d(T, R, kind="cubic", bounds_error=False, fill_value="extrapolate")
    return Tnew, f(Tnew), True


def detect_window_by_derivative(
    T_in,
    RESP_in,
    alpha=0.05,
    sg_frac=DEFAULT_sg_frac,
    sg_poly=DEFAULT_sg_poly,
    wide_frac=DEFAULT_wide_frac,
    rel_env_fallback=DEFAULT_rel_env_fallback,
    min_points_in_window=DEFAULT_min_points_in_window,
    plot=False,
):
    """
    Detecta ventana alrededor del pico REAL del DRS usando criterio:
    |dR/dT| <= alpha * max|dR/dT| (en vecindad del pico).
    """

    # 1) Preproceso
    T0, R0 = preprocess_T_RESP(T_in, RESP_in)

    # 2) Pico REAL del DRS crudo
    idx_peak_real = int(np.nanargmax(R0))
    T_peak_real = float(T0[idx_peak_real])

    # 3) Rejilla uniforme si es necesario
    T_u, R_u, regridded = to_uniform_grid_if_needed(T0, R0)
    T_u = np.asarray(T_u, dtype=float)
    R_u = np.asarray(R_u, dtype=float)
    n = len(T_u)
    if n < 5:
        raise RuntimeError("Pocos puntos en rejilla uniforme (n < 5)")

    # 4) Suavizado SOLO para derivada estable
    wl = ensure_odd(max(5, int(sg_frac * n)))
    if wl >= n:
        wl = ensure_odd(max(3, n - 1))
    poly = min(sg_poly, max(1, wl - 2))
    dt = float(np.mean(np.diff(T_u)))

    try:
        R_s = savgol_filter(R_u, window_length=wl, polyorder=poly, mode="interp")
        dR = savgol_filter(R_u, window_length=wl, polyorder=poly, deriv=1, delta=dt, mode="interp")
    except Exception:
        R_s = R_u.copy()
        dR = np.gradient(R_u, T_u)

    # 5) Índice pico en rejilla uniforme
    idx_peak = int(np.argmin(np.abs(T_u - T_peak_real)))
    T_peak = T_peak_real

    # 6) Vecindad alrededor del pico
    lo = T_peak * (1 - wide_frac)
    hi = T_peak * (1 + wide_frac)
    inds_search = np.where((T_u >= lo) & (T_u <= hi))[0]
    if inds_search.size < 5:
        i0 = max(0, idx_peak - 10)
        i1 = min(n - 1, idx_peak + 10)
        inds_search = np.arange(i0, i1 + 1)

    # 7) Umbral
    max_abs_dR = float(np.nanmax(np.abs(dR[inds_search])))
    threshold = float(alpha) * max_abs_dR
    mask_good = np.abs(dR) <= threshold

    # 8) Expandir desde pico
    iL = idx_peak
    while iL > 0 and mask_good[iL]:
        iL -= 1
    if not mask_good[iL]:
        iL += 1

    iR = idx_peak
    while iR < n - 1 and mask_good[iR]:
        iR += 1
    if not mask_good[iR]:
        iR -= 1

    # 9) Validación / fallback
    if iL <= iR and (iR - iL + 1) >= min_points_in_window:
        T_left = float(T_u[iL])
        T_right = float(T_u[iR])
        method = "derivative_fixed_peak"
    else:
        try:
            env = np.abs(hilbert(R_s))
            peak_env = env[idx_peak]

            li = None
            for i in range(idx_peak, -1, -1):
                if env[i] < rel_env_fallback * peak_env:
                    li = i + 1
                    break
            if li is None:
                li = max(0, idx_peak - int(np.ceil(min_points_in_window / 2)))

            ri = None
            for i in range(idx_peak, n):
                if env[i] < rel_env_fallback * peak_env:
                    ri = i - 1
                    break
            if ri is None:
                ri = min(n - 1, idx_peak + int(np.ceil(min_points_in_window / 2)))

            if ri - li + 1 < min_points_in_window:
                half = max(min_points_in_window // 2, 1)
                li = max(0, idx_peak - half)
                ri = min(n - 1, idx_peak + half)

            T_left = float(T_u[li])
            T_right = float(T_u[ri])
            method = "envelope_fallback_fixed_peak"
        except Exception:
            half = max(min_points_in_window // 2, 3)
            i0 = max(0, idx_peak - half)
            i1 = min(n - 1, idx_peak + half)
            T_left = float(T_u[i0])
            T_right = float(T_u[i1])
            method = "min_symmetric_fixed_peak"

    # 10) info extra
    n_points_orig = int(np.sum((T0 >= T_left) & (T0 <= T_right)))

    return dict(
        T_left=T_left,
        T_right=T_right,
        T_peak=T_peak,
        idx_peak=idx_peak,
        T_u=T_u,
        R_u=R_u,
        R_s=R_s,
        dR=dR,
        mask_good=mask_good,
        method=method,
        regridded=regridded,
        n_points_orig=n_points_orig,
        threshold=threshold,
    )


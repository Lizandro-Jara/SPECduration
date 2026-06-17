import numpy as np

from .newmark import Newmark


# ----------------------------------------------------------------------
# Displacement response spectrum
# ----------------------------------------------------------------------

def maxi(x):
    """
    Return the maximum absolute value of a response time history.
    """
    return float(np.max(np.abs(x)))


def compute_drs(
    ug,
    dt,
    zeta=0.05,
    Tmin=None,
    Tmax=15.0,
    dT=0.01,
):
    """
    Compute the displacement response spectrum, S(T), using Newmark-beta
    time integration.

    Parameters
    ----------
    ug : array-like
        Ground acceleration time history, in m/s².
    dt : float
        Time step of the input record, in seconds.
    zeta : float, optional
        Damping ratio. Default is 0.05.
    Tmin : float, optional
        Minimum oscillator period. If None, it is set automatically.
    Tmax : float, optional
        Maximum oscillator period. Default is 15.0 s.
    dT : float, optional
        Period step used to compute the spectrum. Default is 0.01 s.

    Returns
    -------
    T_values : ndarray
        Period vector, in seconds.
    RESP : ndarray
        Displacement response spectrum values, S(T).
    """
    ug = np.asarray(ug, dtype=float)

    if ug.size < 2:
        raise ValueError("At least two acceleration points are required to compute the DRS.")

    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("The time step dt must be a positive finite value.")

    if not np.isfinite(Tmax) or Tmax <= 0:
        raise ValueError("Tmax must be a positive finite value.")

    if not np.isfinite(dT) or dT <= 0:
        raise ValueError("dT must be a positive finite value.")

    if Tmin is None:
        Tmin = max(10.0 * dt, 0.01)

    if not np.isfinite(Tmin) or Tmin <= 0:
        raise ValueError("Tmin must be a positive finite value.")

    if Tmin >= Tmax:
        raise ValueError("Tmin must be smaller than Tmax.")

    T_values = np.arange(Tmin, Tmax + dT, dT, dtype=float)
    RESP = np.zeros_like(T_values)

    for i, Tn in enumerate(T_values):
        try:
            u, up, upp = Newmark(
                ug,
                0.0,
                0.0,
                0.5,
                0.25,
                dt,
                zeta,
                float(Tn),
            )
            RESP[i] = maxi(u)

        except Exception:
            RESP[i] = np.nan

    return T_values, RESP

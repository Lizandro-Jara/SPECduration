import numpy as np
from .newmark import Newmark

def maxi(x):
    return float(np.max(np.abs(x)))

def compute_drs(
    ug,
    dt,
    zeta=0.05,
    Tmin=None,
    Tmax=15.0,
    dT=0.01
):
    """
    Calcula el Espectro de Respuesta de Desplazamiento (DRS).

    Retorna:
      T_values (np.ndarray)
      RESP (np.ndarray)
    """

    ug = np.asarray(ug, dtype=float)

    if Tmin is None:
        Tmin = max(10*dt, 0.01)

    T_values = np.arange(Tmin, Tmax + dT, dT, dtype=float)
    RESP = np.zeros_like(T_values)

    for i, Tn in enumerate(T_values):
        try:
            u, up, upp = Newmark(ug, 0.0, 0.0, 0.5, 0.25, dt, zeta, float(Tn))
            RESP[i] = maxi(u)
        except Exception:
            RESP[i] = np.nan

    return T_values, RESP


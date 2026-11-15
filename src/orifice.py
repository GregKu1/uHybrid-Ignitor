import numpy as np
from CoolProp.CoolProp import PropsSI

FLUID = "NitrousOxide"

def getMassFlux(P, h0, s0):
    h = PropsSI('H', 'S', s0, 'P', P, FLUID)
    rho = PropsSI('D', 'S', s0, 'P', P, FLUID)
    V = np.sqrt(max(2*(h0-h), 0.0))
    G = rho*V
    return G

def orifice_mdot(Pc_bar, T_tank, D, Cd):
    A = np.pi*D**2 / 4
    P0 = PropsSI('P', 'T', T_tank, 'Q', 1, FLUID)
    h0 = PropsSI('H', 'T', T_tank, 'Q', 1, FLUID)
    s0 = PropsSI('S', 'T', T_tank, 'Q', 1, FLUID)

    P_sweep = np.linspace(0.3*P0, 0.999*P0, 200)
    G_sweep = np.array([getMassFlux(P, h0, s0) for P in P_sweep])

    idx_max = np.argmax(G_sweep)
    P_star = P_sweep[idx_max]
    G_max = G_sweep[idx_max]

    if Pc_bar * 1e5 <= P_star:    # choked
        G = G_max
    else:   # not choked
        G = getMassFlux(Pc, h0, s0)

    mdot = G * Cd * A
    mdot_choked = G_max * Cd * A
    return mdot, mdot_choked


if __name__ == "__main__":
    mdot, mdot_choked = orifice_mdot(Pc_bar=5, T_tank=273.15+15, D=0.5e-3, Cd=0.8)
    # print(f"P0 = {P0/1e5:.3f}bar")
    # print(f"P* = {P_star/1e5:.3f} bar")
    # print(f"P0/P* = {P0/P_star:.3f}")
    # print(f"G_max = {G_max} kg/m^2/s")

    print(f"choked mdot = {mdot_choked*1000:.3f} g/s")
    print(f"mdot = {mdot*1000:.3f} g/s")
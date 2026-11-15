from orifice import orifice_mdot
from combustion import create_cea_object, run_design
import numpy as np

Pc = 5
D_orifice = 0.5e-3
Cd_orifice = 0.8
T_tank = 273.15+5
rho_ABS = 1040
dt = 1e-3
OD = 20e-3


a = 0.038e-3  # wikipedia https://en.wikipedia.org/wiki/Hybrid_rocket_fuel_regression
n = 0.5

cea = create_cea_object()
results = run_design(cea, Pc, D_orifice, 0.8, T_tank)

mdot_ox = results["mdot_ox"]
mdot_fuel = results["mdot_f"]

d0 = 10e-3
G0 = 4 * mdot_ox / (np.pi*d0**2)

r0 =  a * G0**n
L = mdot_fuel / (rho_ABS * np.pi * d0 * r0)
print(f"length of grain = {L*1000:.2f} mm")

d = d0
t_burn = 0.0
fuel_mass = 0.0
while d <= OD:
    G = 4* mdot_ox / (np.pi * d**2)
    r = a * G**n
    mdot_fuel = rho_ABS * np.pi * d * L * r
    fuel_mass += mdot_fuel*dt
    d = d + 2 * r * dt
    t_burn += dt

print(f"burn time = {t_burn:.2f} s")
print(f"fuel mass burnt = {fuel_mass*1000:.2f} g ")


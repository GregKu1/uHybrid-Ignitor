from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt
from orifice import orifice_mdot


# ---------------------------------------------------------------------------
# constant data
# ---------------------------------------------------------------------------

# formation enthalpy [kJ/mol]
hf = {"CO": -110.5, "CO2": -393.5, "H2O": -241.8, "N2": 0.0, "H2": 0.0,
      "O2": 0.0, "H": 218.0, "O": 249.2, "OH": 37.4, "NO": 91.3,
      "HCN": 135.1, "CH4": -74.6, "N": 472.7, "HO2": 12.26, "NH3": -45.94}
# molar masses [g/mol]
M  = {"CO": 28.01, "CO2": 44.01, "H2O": 18.02, "N2": 28.01, "H2": 2.016,
      "O2": 32.00, "H": 1.008, "O": 16.00, "OH": 17.01, "NO": 30.01,
      "HCN": 27.03, "CH4": 16.04, "N": 14.007, "HO2": 33.006, "NH3": 17.031}

abs_card = """
fuel ABS  C 3.85 H 4.85 N 0.15  wt%=100.
h,kj/mol=62.  t(k)=298.15
"""


# ---------------------------------------------------------------------------
# CEA
# ---------------------------------------------------------------------------

def create_cea_object():
    # register the custom fuel with CEA, then build the CEA object
    add_new_fuel("ABS", abs_card)
    cea = CEA_Obj(oxName="N2O", fuelName="Isopropanol", pressure_units="Bar", temperature_units="degK", cstar_units="m/s", enthalpy_units="J/kg")
    return cea


# ---------------------------------------------------------------------------
# CEA lookups
# ---------------------------------------------------------------------------

def sweep_T_comb(cea, OF, Pc):
    T_comb = cea.get_Tcomb(Pc, OF)
    return T_comb

def sweep_C_star(cea, OF, Pc):
    C_star = cea.get_Cstar(Pc, OF)
    return C_star

def sweep_species(cea, OF, Pc):
    _, massfrac = cea.get_SpeciesMassFractions(Pc, OF, 10)
    return massfrac

def h_sensible(cea, Pc, OF, eps=3.34):
    h_chamber = cea.get_Enthalpies(Pc, OF, eps)[0]  # get enthalpies of compounds in the chamber
    _, mfrac = cea.get_SpeciesMassFractions(Pc, OF, eps)    # get mass fractions of those compounds

    mfrac = {k.lstrip("*"): v[0] for k, v in mfrac.items()} # some of the keys have * in the name, so rebuild each "key": value pair without *
    total = sum(mfrac.values()) # get the sum of all mass fractions for normalisation

    h_form = 0.0
    for s in mfrac:
        if s in hf:
            h_form += mfrac[s] / total * hf[s] * 1e6 / M[s] # h_form = sum( w_i · hf_i / M_i )
        else:
            print("no hf for", s, mfrac[s])

    return h_chamber - h_form


# ---------------------------------------------------------------------------
# processing
# ---------------------------------------------------------------------------

def build_species_fractions(species_mass_fractions):
    species = set()
    for d in species_mass_fractions:
        for key in d:
            species.add(key)
    species = sorted(species)

    fractions = {}
    for compound in species:
        values = []
        for d in species_mass_fractions:
            values.append(d.get(compound, [0.0])[0])
        fractions[compound] = np.array(values)

    major_species = []
    for compound in species:
        if fractions[compound].max() > 0.01:
            major_species.append(compound)

    return species, fractions, major_species


# ---------------------------------------------------------------------------
# main calculation
# ---------------------------------------------------------------------------

def run_design(cea, Pc, D_orifice, Cd_orifice, T_ambient, OF_start=1, OF_stop=10, OF_step=0.1):
    OF = np.arange(OF_start, OF_stop, OF_step)

    T_comb = np.array([sweep_T_comb(cea, MR, Pc) for MR in OF])
    C_star = np.array([sweep_C_star(cea, MR, Pc) for MR in OF])
    species_mass_fractions = [sweep_species(cea, MR, Pc) for MR in OF]

    idx_C_star_max = np.argmax(C_star)
    C_star_max = C_star[idx_C_star_max]
    T = T_comb[idx_C_star_max]
    design_OF = OF[idx_C_star_max]
    actual_OF = 2

    mdot_ox = 0.48695
    mdot_f = 0.24347

    P_kW = (mdot_ox + mdot_f) * h_sensible(cea, Pc, actual_OF) / 1e3

    species, fractions, major_species = build_species_fractions(species_mass_fractions)

    results = {
        "Pc": Pc,
        "OF": OF,
        "T_comb": T_comb,
        "C_star": C_star,
        "species_mass_fractions": species_mass_fractions,
        "idx_C_star_max": idx_C_star_max,
        "C_star_max": C_star_max,
        "T": T,
        "design_OF": design_OF,
        "mdot_ox": mdot_ox,
        "mdot_f": mdot_f,
        "P_kW": P_kW,
        "species": species,
        "fractions": fractions,
        "major_species": major_species,
    }
    return results


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

def print_results(results):
    C_star_max = results["C_star_max"]
    T = results["T"]
    design_OF = results["design_OF"]
    P_kW = results["P_kW"]

    print(f"maximum c* = {C_star_max:.2f} m/s , combustion temperature = {T:.2f} K at OF = {design_OF:.2f}")
    print(f"power output = {P_kW:.2f} kW")


def plot_results(results):
    OF = results["OF"]
    C_star = results["C_star"]
    T_comb = results["T_comb"]
    idx_C_star_max = results["idx_C_star_max"]
    major_species = results["major_species"]
    fractions = results["fractions"]

    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 8),
                                      layout="constrained")
    top.plot(OF, C_star, color="tab:blue")
    top.axvline(OF[idx_C_star_max], linestyle="dashed", color="black")
    top.set_ylabel("c* [m/s]", color="tab:blue")
    top2 = top.twinx()
    top2.plot(OF, T_comb, color="tab:orange")
    top2.set_ylabel("T_comb [K]", color="tab:orange")

    for c in major_species:
        bottom.plot(OF, fractions[c], label=c.lstrip("*"))
    bottom.axvline(OF[idx_C_star_max], linestyle="dashed", color="black")
    bottom.set_xlabel("O/F")
    bottom.set_ylabel("chamber mass fraction")
    bottom.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))

    plt.show()


# ---------------------------------------------------------------------------
# script entry point
# ---------------------------------------------------------------------------

def main():
    Pc = 20
    #ignore for janky baldrick version
    D_orifice = 0.5e-3
    Cd_orifice = 0.8
    T_ambient = 273.15+20

    cea = create_cea_object()
    results = run_design(cea, Pc, D_orifice, Cd_orifice, T_ambient)
    print_results(results)
    plot_results(results)


if __name__ == "__main__":
    main()
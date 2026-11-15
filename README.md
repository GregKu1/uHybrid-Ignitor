# Theory crafting of a hybrid ignitor for liquid rocket engines
For N2O/x rockets, a static vent is needed for safety. Because you are always venting nitrous, this nitrous could be used in a small scale hybrid ignitor to light the main engine. The biggest challenge here is oxidiser massflow, as the static vent is small on purpose, to minimise nitrous losses.  
- orifice.py calculates massflow through the static vent orifice at varying tank conditions.  
- combustion.py runs CEA to N2O/ABS to find optimal OF and therefore fuel mass flow. It also calculates the heat power output, to compare with the main engine and see if it can even be lit.  
- fuel.py sizes the grain according to grain regression formula.
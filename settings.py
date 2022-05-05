"""
  File to set general settings of the simulation
"""

# Global simulation settings
R = 8.314  # J/molK
NA = 6.022 * 1e23  # mole/mol
kB = 1.38064852 * 1e-23  # m2kgs-2K-1 Boltzmann constant PV/TN
h = 6.62607004 * 1e-34  # m2kgs-1 Planck's constant
T = 433  # K

# Transient simulation settings
N_TIME_STEPS = 100000

# --- Arbitrary calculations ---

C = 0.3  # mol/L
V = 1 * 1e-20  # L
FRUCTOSE_INIT_NUM_MOLES = V * C * NA

# ------------------------------

# Initial number of molecules for each species
initialNumberOfMolecules = [
    FRUCTOSE_INIT_NUM_MOLES,  # fructose
    0.0,  # int1
    0.0,  # int2
    0.0,  # int3
    0.0,  # int4
    0.0,  # int5
    0.0,  # int6
    0.0,  # int7
    0.0,  # HMF
]
N_SPECIES = len(initialNumberOfMolecules)

# Reactions (are there coef-s in the chemical equations?)
reactions = [
    {  # reaction index = 0
        "reactants": [0],
        "products": [1],
        "Ea_f": 12 * 4.18,
        "Ea_b": 9 * 4.18,
    },
    {  # reaction index = 1
        "reactants": [1],
        "products": [2],
        "Ea_f": 6 * 4.18,
        "Ea_b": 8.2 * 4.18,
    },
    {  # reaction index = 2
        "reactants": [2],
        "products": [3],
        "Ea_f": 26.3 * 4.18,
        "Ea_b": 19.6 * 4.18,
    },
    {
        "reactants": [3],
        "products": [4],
        "Ea_f": 17.2 * 4.18,
        "Ea_b": 8.5 * 4.18,
    },
    {
        "reactants": [4],
        "products": [5],
        "Ea_f": 10.5 * 4.18,
        "Ea_b": 36.5 * 4.18,
    },
    {
        "reactants": [5],
        "products": [6],
        "Ea_f": 32 * 4.18,
        "Ea_b": 19.6 * 4.18,
    },
    {
        "reactants": [6],
        "products": [7],
        "Ea_f": 32.2 * 4.18,
        "Ea_b": 66.4 * 4.18,
    },
    {
        "reactants": [7],
        "products": [8],
        "Ea_f": 10.8 * 4.18,
        "Ea_b": 14.4 * 4.18,
    },
]

# Define probabilities settings
# Values in the arrays are the reaction indexes
# the first index defines which reaction rate goes
# in the numerator
probabilitiesMap = [
    [0],  # P1
    [1],  # P2
    [2],
    [3],
    [4],
    [5],
    [6],
    [7],
]
N_PROBABILITIES = len(probabilitiesMap)

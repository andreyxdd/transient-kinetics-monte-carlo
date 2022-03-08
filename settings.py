# Global simulation settings
R  = 8.314 # J/molK
Na = 6.022*1e23 # mole/mol
kBoltzman = 1.38064852*1e-23 # m2kgs-2K-1 Boltzmann constant PV/TN
h  = 6.62607004*1e-34 # m2kgs-1 Planck's constant
T  = 433 # K
h  = 1

# Transient simulation settings
tEnd = 10 # seconds

# Species conentrations object
species = [
  'DHH',    # 0
  'prot1',  # 1
  'int1',   # 2
  'eno112', # 3
  'eno123', # 4
  'prot2',  # 5
  'int2',   # 6
  'eno145', # 7
  'eno156', # 8
]

# Reactions (are there coef-s in the chemical equations?)
reactions = [
  { # reaction index = 0
    'reactants': [0],
    'products': [1],
    'Ea_f': 1.0,
    'Ea_b': 1.0,
  },
  { # reaction index = 1
    'reactants': [0],
    'products': [5],
    'Ea_f': 1.0,
    'Ea_b': 1.0, 
  },
  { # reaction index = 2
    'reactants': [2],
    'Ea_f': 1.0,
  },
  {
    'reactants': [2],
    'Ea_f': 1.0,
  },
  {
    'reactants': [6],
    'Ea_f': 1.0,
  },
  {
    'reactants': [6],
    'Ea_f': 1.0,
  }
]

# Define probabilities settings
# Values in the arrays are the reaction indexes
# the first index defines which reaction rate goes
# in the numerator
probabilities = {
  'P1': [0, 2],
  'P2': [0, 1],
  'P3': [4, 5],
  'P4': [2, 3],
  'P5': [6, 7]
}

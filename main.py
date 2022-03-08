import numpy as np
from scipy.optimize import fsolve
from settings import kBoltzman, h, T, R, reactions

def computeRate(Ea, C, backward = False):
  sign = 1
  if backward:
    sign = - sign
  
  return sign * ( (kBoltzman * T / h) * (np.exp(-Ea / (R * T))) ) * C;

# x - vector species concentrations
def assembleSystem(x):
  reactionRates = np.zeros(len(x))
  for idx, reaction in enumerate(reactions):
    reactionRate = 0.0;
    
    # if a reaction proceeds forward:
    if 'Ea_f' in reaction:
      # iterating over reactants:
      for r in reaction['reactants']:
        reactionRate += computeRate(
          reaction['Ea_f'],
          x[r],
          True
        )
    
    # if a reaction proceeds backwards:
    elif 'Ea_b' in reaction:
      # iterating over products:
      for p in reaction['products']:
        reactionRate += computeRate(
          reaction['Ea_b'],
          x[p]
        )
    
    reactionRates[idx] = reactionRate

  return reactionRates

if __name__ == "__main__":
  initialNumberOfMolecules = [
    200.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
  ]
  
  result = fsolve(assembleSystem, initialNumberOfMolecules)
  reactiionRates = assembleSystem(result)

  # printing new vector of number of molecules  
  print(result)
  
  # printing vector of reaction rates
  print(reactiionRates)
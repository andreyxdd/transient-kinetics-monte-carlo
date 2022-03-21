import numpy as np
from scipy.optimize import fsolve
from settings import (
  kBoltzman, h, T, R,
  initialNumberOfMolecules, reactions, probabilitiesMap,
  tEnd, nTimeSteps )
from runMonteCarlo import runMonteCarlo

def computeRate(Ea, C, backward = False):
  sign = 1
  if backward:
    sign = - sign
  
  return sign * ( (kBoltzman * T / h) * (np.exp(-Ea / (R * T))) ) * C;

# x - vector with number of molecules (or species concentrations)
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
          x[r]
        )
    
    # if a reaction proceeds backwards:
    elif 'Ea_b' in reaction:
      # iterating over products:
      for p in reaction['products']:
        reactionRate += computeRate(
          reaction['Ea_b'],
          x[p],
          True
        )
    
    reactionRates[idx] = reactionRate

  return reactionRates

if __name__ == "__main__":
  # resulting matrix:
  #  - each row represents a time step
  #  - each column  represents a species (number of molecules)
  dimensions = (nTimeSteps, len(initialNumberOfMolecules))
  numberOfMoleculesWithTime = np.zeros(dimensions)
  
  # fill the first row with initial data
  numberOfMoleculesWithTime[0] = initialNumberOfMolecules
  
  # array of time steps:
  # tArray = np.linspace(0, tEnd, nTimeSteps)
  # TODO: HOW EXACTLY DOES THE TIME STEP AFFECTS THE SIMULATION?
  # TODO: use another criteri to go out of the loop
  
  for i in range(1, nTimeSteps):
    
    # ---
    # solving system of ODEs by using the species
    # vector for the previous time step 'i - 1' 
    currentNumberOfMolecules = fsolve(
      assembleSystem,
      numberOfMoleculesWithTime[i-1]
    )
    # ---
    
    # ---
    # computing the reaction rates using the obtained
    # solution vector
    reacttionRates = assembleSystem(currentNumberOfMolecules)
    # ---
    
    # ---
    # computing reactions' probabilites by using
    # scheme set in the settings file
    probabilities = np.array([])
    for pMap in probabilitiesMap:

      reactionRatesSum = 1.0e-20
      for idx in pMap:
        reactionRatesSum += reacttionRates[idx]
      
      np.append(
        probabilities,
        reacttionRates[pMap[0]]/reactionRatesSum
      )
    # ---
    
    # ---
    # Running Monte-Carlo simulation by passing the
    # current vectors with species and probabilities
    newNumberOfMolecules = runMonteCarlo(
      currentNumberOfMolecules,
      probabilities
    )
    # ---
    
    # ---
    # adjusting the values for the current step 'i'
    numberOfMoleculesWithTime[i] = newNumberOfMolecules
    # ---
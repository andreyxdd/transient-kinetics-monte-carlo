import numpy as np
import matplotlib.pyplot as plt
from settings import (
  kB, h, T, R,
  initialNumberOfMolecules, reactions, probabilitiesMap,
  tEnd, nTimeSteps, nSpecies )
from runMonteCarlo import runMonteCarlo

def computeRate(Ea, C):
  return ( (kB * T / h) * (np.exp(-Ea / (R * T))) ) * C;

if __name__ == "__main__":
  # resulting matrix:
  #  - each row represents a time step
  #  - each column  represents a species (number of molecules)
  dimensions = (nTimeSteps, len(initialNumberOfMolecules))
  numberOfMoleculesWithTime = np.zeros(dimensions)
  
  # fill the first row with initial data
  numberOfMoleculesWithTime[0] = initialNumberOfMolecules
  
  # array of time steps:
  cutIndex = 0 # last index designating the total number of iterations
  for i in range(1, nTimeSteps):

    currentNumberOfMolecules = numberOfMoleculesWithTime[i-1]

    reactionRates = np.zeros((len(reactions), 2))
    for idx, r in enumerate(reactions):
      for rt in r['reactants']:
        reactionRates[idx, 0] += computeRate(
          r['Ea_f'],
          currentNumberOfMolecules[rt]
        )
      for pt in r['products']:
        reactionRates[idx, 1] += computeRate(
          r['Ea_b'],
          currentNumberOfMolecules[pt]
        )

    # ---
    # computing reactions' probabilites by using
    # scheme set in the settings file
    probabilities = np.zeros(len(probabilitiesMap))
    for idxP, pMap in enumerate(probabilitiesMap):
      reactionRatesSum = 1.0e-25
      for idxR in pMap:
        # only pMap[0] is taken into account
        reactionRatesSum += reactionRates[idxR, 0] + reactionRates[idxR, 1] 

      probabilities[idxP] = reactionRates[pMap[0], 0]/reactionRatesSum
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

    cutIndex = i
    
    # TODO: Think about condition to finish the loop

  #--- Plotting and outputting the image
  speciesLabel = ''
  for j in range(nSpecies):
    if j == 0:
      speciesLabel = 'Fructose'
    elif j == nSpecies-1:
      speciesLabel = 'HMF'
    else:
      speciesLabel = 'int-'+str(j)

    plt.plot(numberOfMoleculesWithTime[:cutIndex+1,j], label=speciesLabel)

  plt.title(r'Species number of molecules', fontsize='16')
  plt.ylabel(r'Number of molecules',fontsize='13')
  plt.xlabel("# iterations",fontsize='13')
  plt.legend()
  plt.grid()
  plt.savefig('species.png')
  # plt.show()
  # ---
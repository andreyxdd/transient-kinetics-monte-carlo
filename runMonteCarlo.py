import numpy as np
from settings import reactions, probabilitiesMap

nProbabilities = len(probabilitiesMap)

def runMonteCarlo(speciesVector, probabilities, currProbabNum = 0):
  if(currProbabNum == nProbabilities - 1):
    return speciesVector

  if(probabilities[currProbabNum] >= np.random.random()):
    reactionIds = probabilitiesMap[currProbabNum]

    for reactionId in reactionIds:
      reactantsIds = reactions[reactionId]['reactants']
      for reactantId in reactantsIds:
        speciesVector[reactantId] -= 1

      productsIds = reactions[reactionId]['products']
      for productId in productsIds:
        speciesVector[productId] += 1

  currProbabNum += 1
  return runMonteCarlo(speciesVector, probabilities, currProbabNum)
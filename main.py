"""
  Main script to run simulation
"""
import numpy as np
import matplotlib.pyplot as plt
from settings import (
    kB,
    h,
    T,
    R,
    initialNumberOfMolecules,
    reactions,
    probabilitiesMap,
    N_TIME_STEPS,
    N_SPECIES,
)
from run_monte_carlo import run_monte_carlo


def compute_rate(activation_energy, concnetration):
    """
    Function to compute reaction rate (in one direction) accroding
    to the Arrhenius law

    Args:
        Ea (number): activation energy for the given reaction
        C (number): the species concentration

    Returns:
        number: reaction rate
    """
    return ((kB * T / h) * (np.exp(-activation_energy / (R * T)))) * concnetration


if __name__ == "__main__":
    # resulting matrix:
    #  - each row represents a time step
    #  - each column  represents a species (number of molecules)
    dimensions = (N_TIME_STEPS, len(initialNumberOfMolecules))
    numberOfMoleculesWithTime = np.zeros(dimensions)

    # fill the first row with initial data
    numberOfMoleculesWithTime[0] = initialNumberOfMolecules

    # array of time steps:
    CUT_INDEX = 0  # last index designating the total number of iterations
    for i in range(1, N_TIME_STEPS):

        currentNumberOfMolecules = numberOfMoleculesWithTime[i - 1]

        reactionRates = np.zeros((len(reactions), 2))
        for idx, r in enumerate(reactions):
            for rt in r["reactants"]:
                reactionRates[idx, 0] += compute_rate(
                    r["Ea_f"], currentNumberOfMolecules[rt]
                )
            for pt in r["products"]:
                reactionRates[idx, 1] += compute_rate(
                    r["Ea_b"], currentNumberOfMolecules[pt]
                )

        # ---
        # computing reactions' probabilites by using
        # scheme set in the settings file
        probabilities = np.zeros(len(probabilitiesMap))
        for idxP, pMap in enumerate(probabilitiesMap):
            RATE_SUM = 1.0e-25
            for idxR in pMap:
                # only pMap[0] is taken into account
                RATE_SUM += reactionRates[idxR, 0] + reactionRates[idxR, 1]

            probabilities[idxP] = reactionRates[pMap[0], 0] / RATE_SUM
        # ---

        # ---
        # Running Monte-Carlo simulation by passing the
        # current vectors with species and probabilities
        newNumberOfMolecules = run_monte_carlo(currentNumberOfMolecules, probabilities)
        # ---

        for idx, n in enumerate(newNumberOfMolecules):
            if n <= 0:
                newNumberOfMolecules[idx] = 0.0

        # ---
        # adjusting the values for the current step 'i'
        numberOfMoleculesWithTime[i] = newNumberOfMolecules
        # ---

        CUT_INDEX = i

    # --- Plotting and outputting the image
    SPECIES_LABEL = ""
    for j in range(N_SPECIES):
        if j == 0:
            SPECIES_LABEL = "Fructose"
        elif j == N_SPECIES - 1:
            SPECIES_LABEL = "HMF"
        else:
            SPECIES_LABEL = "int-" + str(j)

        plt.plot(numberOfMoleculesWithTime[: CUT_INDEX + 1, j], label=SPECIES_LABEL)

    plt.title(r"Species number of molecules vs reaction direction", fontsize="16")
    plt.ylabel(r"Number of molecules", fontsize="13")
    plt.xlabel("# iterations", fontsize="13")
    plt.legend(loc="upper right")
    plt.grid()
    plt.savefig("species.png")
    # plt.show()
    # ---

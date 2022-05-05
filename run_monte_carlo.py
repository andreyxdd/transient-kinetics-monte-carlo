"""
    File for function implementing the Monte-Carlo method
"""
import numpy as np
from settings import reactions, probabilitiesMap, N_PROBABILITIES


def run_monte_carlo(species_vector, probabilities, curr_probab_num=0):
    """
    Recursive function based on the Monte-Carlo method for chemical reaction mechanisms

    Args:
        species_vector (vector<number>): vector with all species number of molecules
        probabilities (vector<number>): array of the reaction rate probabilites
        curr_probab_num (int, optional): recursive incrementor. Defaults to 0.

    Returns:
        vector<number>: updated vector with all species number of molecules
    """
    if curr_probab_num == N_PROBABILITIES:
        return species_vector

    if probabilities[curr_probab_num] >= np.random.random():
        reaction_ids = probabilitiesMap[curr_probab_num]

        for reaction_id in reaction_ids:
            reactants_ids = reactions[reaction_id]["reactants"]
            for reactant_id in reactants_ids:
                species_vector[reactant_id] -= 1

            products_ids = reactions[reaction_id]["products"]
            for product_id in products_ids:
                species_vector[product_id] += 1

    return run_monte_carlo(species_vector, probabilities, curr_probab_num + 1)

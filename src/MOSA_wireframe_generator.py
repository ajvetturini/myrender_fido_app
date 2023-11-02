"""
AJ Vetturini
IDIG and MMBL
Carnegie Mellon University
Advised By: Jon Cagan and Rebecca Taylor

This script is in charge of running either individual optimization problem setups or running a batch job of files using
the Multi-Objective Simulated Annealing (MOSA) algorithm.

Note to self: Copy code over from "wireframe_generator_3D.py" to incorporate the multi-processing!
"""
# Import Modules:
def find_min_and_maxes(l1: list, l2: list) -> list:
    """
    This function finds the min and max values for the repulsion, buckling, and faces between the two lists so that the
    range is FIXED in both plots to understand how MOSA searches the design space
    """
    buckling, repulsion, faces = [list(col) for col in zip(*l1)]
    buckling2, repulsion2, faces2 = [list(col) for col in zip(*l2)]

    b_max, b_min = max([max(buckling), max(buckling2)]), min([min(buckling), min(buckling2)])
    r_max, r_min = max([max(repulsion), max(repulsion2)]), min([min(repulsion), min(repulsion2)])
    f_max, f_min = max([max(faces), max(faces2)]), min([min(faces), min(faces2)])

    return [b_max, b_min, r_max, r_min, f_max, f_min]



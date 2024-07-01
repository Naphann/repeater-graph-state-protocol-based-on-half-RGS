#! /usr/bin/python3

from enum import Enum
from typing import Self

import numpy as np

Pauli = Enum("Pauli", ["I", "X", "Y", "Z"])

"""The updates of measurement trees occur from the following:

1. side effect from generations (toggling of pauli Z tracking)
    1.1 push-out operations -- affecting all but leaf nodes of inner physical qubits
    1.2 XX measurements (fusing outer and inner qubits) -- affecting outer qubit and 1st level vertices
    1.3 XX measurements (fusing two half-RGSs) -- affecting 1st level vertices 
2. The update of 
"""


class MeasurementTree:

    class Node:

        def __init__(self, parent_node: Self):
            self.measurement_result: bool | None = None
            self.measurement_basis: Pauli | None = None
            self.eigenvalue: bool | None = None
            self.parent: Self | None = parent_node
            self.children: list[Self] = []
            self.is_lost = False
            self.has_z = False

    def __init__(self, m: int, branching_vector: list[int]):
        # generate the tree structure with node
        pass

    def reset(self):
        # is called to reset measurement results and reinitialized before each RGS trial run
        pass

    def parse_postorder_side_effect(self, zs: list[bool]):
        # read and assign the has_z
        pass

    def parse_postorder_lost(self, photon_losts: list[bool]):
        pass


    
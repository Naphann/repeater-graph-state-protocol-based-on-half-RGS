#! usr/bin/python3

from enum import Enum

import numpy as np

Pauli = Enum("Pauli", ["I", "X", "Y", "Z"])


def num_qubits_per_rgs_arm(bvec: list[int]) -> int:
    num_in_layers = bvec[:]
    for i in range(1, len(num_in_layers)):
        num_in_layers[i] *= num_in_layers[i - 1]
    num_qubits_per_arm = 1 + np.sum(num_in_layers)
    return num_qubits_per_arm


class Node:
    """Represent a node in the tree where the actual root of the tree has a special meaning.
    The root corresponds to the physical outer qubit while all other nodes are the physical qubits
    of the tree code encoding the inner qubits."""

    def __init__(self, qubit_index=-1, parent_index=-1):
        self.qubit_index = qubit_index
        self.parent_index = parent_index  # for debugging
        self.measurement_result: bool | None = None  # this should be True and False if the qubit has been measured indicating the raw measurement result
        self.eigenvalue: bool | None = None  # use this to track the decoded measurement (after taking side effect and raw results into account)
        self.measurement_basis: Pauli | None = None
        self.children: list[Node] = []
        self.is_lost = False  # this is used to denote whether the qubit is lost in the fiber or not
        self.has_z = False  # this is used to denote whether the qubit has Z side effect from the emission process or not (only from the emission!!)

    def get_level_traversal(self):
        return_list = [self]
        queue = [self]
        while len(queue) > 0:
            new_queue = []
            for u in queue:
                return_list.append(u)
                new_queue.extend(u.children)
            queue = new_queue
        return return_list

    def get_indices_from_level(self, k: int) -> list[int]:
        """get all the indices from the nodes"""
        cur_level = 0
        queue = [self]
        while len(queue) > 0:
            # arrive at the correct level, return the qubit indices
            if cur_level == k:
                return [v.qubit_index for v in queue]
            # need to go to the next level
            num_nodes_in_level = len(queue)
            for _ in range(num_nodes_in_level):
                u = queue.pop(0)
                queue.extend(u.children)
            cur_level += 1
        # should change this to an error
        return []  # empty list indicating the level specified is out of range

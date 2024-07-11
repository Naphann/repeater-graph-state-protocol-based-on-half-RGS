#! /usr/bin/python

from enum import Enum
from typing import Self

import numpy as np
import stim

Pauli = Enum("Pauli", ["I", "X", "Y", "Z"])


class Node:
    """Represent a node in the tree where the actual root of the tree has a special meaning.
    The root corresponds to the physical outer qubit while all other nodes are the physical qubits
    of the tree code encoding the inner qubits."""

    def __init__(self):
        self.measurement_result: bool | None = None  # raw measurement results, None if not measured
        self.eigenvalue: bool | None = None  # use to compute parity, corrected results after taking side effects into account
        self.measurement_basis: Pauli | None = None
        self.children: list[Node] = []
        self.is_lost = False  # denote whether the qubit is lost in the fiber
        self.has_z = False  # denote whether the qubit has Z side effect from the emission process

    def get_postorder_traversal(self) -> list[Self]:
        ret: list[Self] = []

        def __inner_postorder_recurse(u: Node):
            for v in u.children:
                __inner_postorder_recurse(v)
            ret.append(u)

        __inner_postorder_recurse(self)
        return ret

    def reset(self):
        self.measurement_result = None
        self.eigenvalue = None
        self.measurement_basis = None
        self.is_lost = False
        self.has_z = False

        for u in self.children:
            u.reset()


class RgsConfig:

    def __init__(self, number_of_hops: int, m: int, bv: list[int], loss_probability: float, tab_sim: stim.TableauSimulator):
        self.rng = np.random.default_rng()
        self.t = tab_sim

        self.m = m
        self.bv = bv
        self.loss_probability = loss_probability
        self.number_of_hops = number_of_hops

        # # minimal circuit
        # self.alice = 0
        # self.bob = 3
        # self.outer_emitter_left = 1
        # self.outer_emitter_right = 4
        # self.photon = 6
        # self.photon_left = 5
        # self.photon_right = 6
        # self.emitters = [2 + i for i in range(len(bv))]

        # default
        self.alice = 0
        self.bob = 1
        self.anchor_left = 2
        self.anchor_right = 3
        self.outer_emitter_left = 4
        self.outer_emitter_right = 5
        self.photon = 6
        self.photon_left = 6
        self.photon_right = 7
        self.emitters = [8 + i for i in range(len(bv))]

        # data structures for data
        self.measurement_trees: list[list[Node]] = [[Node() for _ in range(m)] for _ in range(2 * number_of_hops)]
        self.logical_results: list[list[None | bool]] = [[None for _ in range(m)] for _ in range(2 * number_of_hops)]
        self.succeeded_bsm_arm_indices = [-1 for _ in range(number_of_hops)]

        # debugging variables
        self.lost_photons = 0
        self.total_photons = 0
        self.correct_bell_pair_count = 0
        self.incorrect_bell_pair_count = 0
        self.unentangled_pair_count = 0

        # debugging circuit
        self.circuit = stim.Circuit()

        # initialization of measurement trees
        for arms in self.measurement_trees:
            for root in arms:
                queue = [root]
                for bi in bv:
                    temp_queue = []
                    for u in queue:
                        for _ in range(bi):
                            v = Node()
                            u.children.append(v)
                            temp_queue.append(v)
                    queue = temp_queue

    def reset(self):
        self.t.reset(*range(self.emitters[-1] + 1))
        self.t.h(0, 1, 2, 3, 4, 5, *self.emitters)
        # self.t.reset(*range(8))
        # self.t.h(0, 1, 3, 4, *self.emitters)
        self.logical_results: list[list[None | bool]] = [[None for _ in range(self.m)] for _ in range(2 * self.number_of_hops)]
        self.succeeded_bsm_arm_indices = [-1 for _ in range(self.number_of_hops)]
        self.circuit = stim.Circuit("""
        H 0 1 2 3 4 5
        """)
        self.circuit.append("H", self.emitters)
        # self.circuit = stim.Circuit("""
        # H 0 1 3 4
        # """)
        # self.circuit.append("H", self.emitters)

        for arms in self.measurement_trees:
            for root in arms:
                root.reset()

    def reset_debug_statistics(self):
        self.lost_photons = 0
        self.total_photons = 0
        self.correct_bell_pair_count = 0
        self.incorrect_bell_pair_count = 0
        self.unentangled_pair_count = 0

#! usr/bin/python3

from typing import Self
from node_qubit import Pauli

import numpy as np
import stim


class SideEffectNode:
    """Represent a node in the side effect tree where the actual root of the tree has a special meaning.
    The root corresponds to the physical outer qubit while all other nodes are the physical qubits
    of the tree code encoding the inner qubits."""

    def __init__(self):
        self.children: list[SideEffectNode] = []
        self.has_z = False  # this is used to denote whether the qubit has Z side effect or not

    def reset(self):
        self.has_z = False
        for u in self.children:
            u.reset()


class RGSS:

    def __init__(
        self,
        m: int,
        bv: list[int],
        loss_probability: float,
        emitters: list[int],
        left_anchor: int,
        right_anchor: int,
        left_outer_emitter: int,
        right_outer_emitter: int,
        left_photon: int,
        right_photon: int,
        tableau_simulator: stim.TableauSimulator,
        rng: np.Generator,
    ):
        self.m = m
        self.bv = bv
        self.loss_probability = loss_probability
        self.t = tableau_simulator
        self.rng = rng
        self.emitters = emitters
        self.left_anchor = left_anchor
        self.right_anchor = right_anchor
        self.left_outer_emitter = left_outer_emitter
        self.right_outer_emitter = right_outer_emitter
        self.left_photon = left_photon
        self.right_photon = right_photon
        self.left_side_effects: list[SideEffectNode] = [SideEffectNode() for _ in range(m)]
        self.right_side_effects: list[SideEffectNode] = [SideEffectNode() for _ in range(m)]

        # create the side effect tree
        for root in [*self.left_side_effects, *self.right_side_effects]:
            queue = [root]
            for bi in bv:
                temp_queue = []
                for u in queue:
                    for _ in range(bi):
                        v = SideEffectNode()
                        u.children.append(v)
                        temp_queue.append(v)
                queue = temp_queue

    def reset(self):
        for root in [*self.left_side_effects, *self.right_side_effects]:
            root.reset()

    def emit_left_outer_photon(self) -> int:
        """return the index of generated photon"""
        self.t.reset(self.left_photon)
        self.t.cx(self.left_outer_emitter, self.left_photon)
        self.t.h(self.left_photon)
        return self.left_photon
    
    def emit_left_inner_photon(self, logical_basis: Pauli) -> SideEffectNode:
        



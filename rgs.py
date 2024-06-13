#! usr/bin/python3

import numpy as np
import random
import stim
from node_qubit import Pauli, Node
from test_helper import verify_vertex_stabilizer
from tree_code_helper import decode_tree_logical_x, decode_tree_logical_z


def helper_assign_qubit_indices(root: Node, bv: list[int], starting_index: int) -> int:
    cur_index = starting_index
    root.qubit_index = cur_index
    cur_index += 1
    queue = [root]
    for bi in bv:
        temp_queue = []
        for u in queue:
            for _ in range(bi):
                v = Node(cur_index, u.qubit_index)
                cur_index += 1
                u.children.append(v)
                temp_queue.append(v)
        queue = temp_queue
    # return the next unused index
    return cur_index


def helper_initialize_rgs_arm(t: stim.TableauSimulator, root: Node, anchor: int, outer_emitter: int, root_ancilla: int) -> bool:
    # return whether the anchor should be flipped or not (side effects to the anchor)
    # generate outer qubit
    t.reset(outer_emitter, root_ancilla)
    t.h(root.qubit_index, outer_emitter)
    t.cz(root.qubit_index, outer_emitter)

    # generate inner qubit tree
    t.h(root_ancilla)
    queue = root.children  # nodes in the first level
    for u in queue:
        t.h(u.qubit_index)
        t.cz(root_ancilla, u.qubit_index)
    # assuming that the anchor is already has Hadamard applied
    while len(queue) > 0:
        temp_queue = []
        for u in queue:
            for v in u.children:
                t.h(v.qubit_index)
                t.cz(u.qubit_index, v.qubit_index)
                temp_queue.append(v)
        queue = temp_queue

    # add random side effects to nodes in the tree except the leaves
    queue = root.children
    while len(queue) > 0:
        temp_queue = []
        for u in queue:
            if len(u.children) == 0:
                break
            if random.random() < 0.5:
                t.z(u.qubit_index)
                u.has_z = not u.has_z
            temp_queue.extend(u.children)
        queue = temp_queue

    # verify anchor stabilizer
    verify_vertex_stabilizer(t, root_ancilla, [u.qubit_index for u in root.children], 1)

    # join inner and outer qubits
    t.cz(anchor, outer_emitter)
    t.cz(root_ancilla, outer_emitter)
    t.h(outer_emitter, root_ancilla)
    meas_outer = t.measure(outer_emitter)
    meas_root = t.measure(root_ancilla)
    if meas_outer:
        # flip first level qubits
        for u in root.children:
            u.has_z = not u.has_z

    if meas_root:
        # flip outer qubits (and return the flip to the anchor)
        root.has_z = not root.has_z

    assert meas_root == root.has_z
    verify_vertex_stabilizer(t, root.qubit_index, [u.qubit_index for u in root.children], -1 if root.has_z else 1)

    return meas_root


def helper_process_photon_loss(t: stim.TableauSimulator, root: Node, loss_probability: float, rng: np.random.Generator):
    """traverse the tree and apply loss probability to all qubits
    If a qubit is lost, randomly select Pauli X, Y, or Z to apply followed by a measurement in the Z basis.
    """
    queue = [root]
    while len(queue) > 0:
        temp_queue = []
        for u in queue:
            temp_queue.extend(u.children)
            if rng.random() < loss_probability:
                q = u.qubit_index
                u.is_lost = True
                pauli_op = random.choice(list(Pauli))
                if pauli_op == Pauli.X:
                    t.x(q)
                elif pauli_op == Pauli.Y:
                    t.y(q)
                elif pauli_op == Pauli.Z:
                    t.z(q)
                t.measure(q)
        queue = temp_queue


def helper_update_eigenvalue_with_side_effect(root: Node):
    if (not root.is_lost) and (root.measurement_basis == Pauli.X) and (root.has_z):
        root.eigenvalue = not root.eigenvalue
    for v in root.children:
        helper_update_eigenvalue_with_side_effect(v)


class HalfRGS:
    """Currently this is only used at end nodes"""

    def __init__(self, m: int, branching_params: list[int], anchor_index: int):
        if len(branching_params) == 0:
            raise ValueError("branching parameters cannot be an empty list")
        self.m = m
        self.bv = branching_params
        self.arms = [Node() for _ in range(m)]
        self.anchor = anchor_index
        self.measurement_bases: list[Pauli | None] = [None for _ in range(m)]
        self.successful_arm_index = -1
        self.logical_results: list[bool | None] = [None for _ in range(m)]

    def assign_qubit_indices(self, starting_index: int) -> int:
        """Assign qubits to half RGS
        Return: next unused qubit index"""
        cur_index = starting_index
        for root in self.arms:
            cur_index = helper_assign_qubit_indices(root, self.bv, cur_index)
        return cur_index

    def initialize_quantum_state(self, t: stim.TableauSimulator, outer_emitter: int, root_ancilla: int):
        anchor_has_z = False
        t.h(self.anchor)
        for root in self.arms:
            anchor_has_z = anchor_has_z ^ helper_initialize_rgs_arm(t, root, self.anchor, outer_emitter, root_ancilla)
        if anchor_has_z:
            t.z(self.anchor)
        first_level_qubits = []
        for root in self.arms:
            first_level_qubits.extend([u.qubit_index for u in root.children])
        verify_vertex_stabilizer(t, self.anchor, first_level_qubits, 1)

    def get_bsm_arm(self) -> Node:
        return self.arms[self.successful_arm_index]

    def process_photon_loss(self, t: stim.TableauSimulator, loss_probability: float, rng: np.random.default_rng):
        for root in self.arms:
            helper_process_photon_loss(t, root, loss_probability, rng)

    def update_measurement_with_side_effects(self):
        """using the side effect stored in .has_z to update .eigenvalues"""
        for u in self.arms:
            helper_update_eigenvalue_with_side_effect(u)

    def decode_logical_results(self) -> bool:
        """decoding the logical measurements of the inner qubits
        Returns False when the decoding fail and the trial needs retried"""
        flag = True
        for i, u in enumerate(self.arms):
            if i == self.successful_arm_index:
                self.logical_results[i] = decode_tree_logical_x(u)
                if decode_tree_logical_x(u) is None:
                    flag = False
            else:
                self.logical_results[i] = decode_tree_logical_z(u)
                if decode_tree_logical_z(u) is None:
                    flag = False
        if all(map(lambda res: res is not None, self.logical_results)) != flag:
            raise RuntimeError("Fuck!")
        return all(map(lambda res: res is not None, self.logical_results))

    def count_lost_photons(self) -> tuple[int, int]:
        """return (lost_photons, total_photons)"""
        total_photons = 0
        lost_photons = 0
        for root in self.arms:
            queue = [root]
            while len(queue) > 0:
                temp_queue = []
                for u in queue:
                    total_photons += 1
                    lost_photons += 1 if u.is_lost else 0
                    temp_queue.extend(u.children)
                queue = temp_queue
        return lost_photons, total_photons


class RGS:
    def __init__(self, m: int, branching_params: list[int]):
        if len(branching_params) == 0:
            raise ValueError("branching parameters cannot be an empty list")
        self.m = m
        self.bv = branching_params
        self.left_arms = [Node() for _ in range(m)]
        self.right_arms = [Node() for _ in range(m)]
        self.measurement_bases_left: list[Pauli | None] = [None for _ in range(m)]
        self.measurement_bases_right: list[Pauli | None] = [None for _ in range(m)]
        self.successful_left_arm_index = -1
        self.successful_right_arm_index = -1
        self.left_logical_results: list[bool | None] = [None for _ in range(m)]
        self.right_logical_results: list[bool | None] = [None for _ in range(m)]

    def assign_qubit_indices(self, starting_index: int) -> int:
        """Assign qubits to be used for the photonic qubits of the RGS
        Return: next unused qubit index"""
        cur_index = starting_index
        for root in self.left_arms:
            cur_index = helper_assign_qubit_indices(root, self.bv, cur_index)
        for root in self.right_arms:
            cur_index = helper_assign_qubit_indices(root, self.bv, cur_index)
        return cur_index

    def initialize_quantum_state(self, t: stim.TableauSimulator, anchor_left: int, anchor_right: int, outer_emitter: int, root_ancilla: int):
        # make sure the qubits are properly initialized
        t.reset(anchor_left, anchor_right)

        # generate the left arms
        anchor_has_z = False
        t.h(anchor_left)
        for root in self.left_arms:
            anchor_has_z = anchor_has_z ^ helper_initialize_rgs_arm(t, root, anchor_left, outer_emitter, root_ancilla)
        if anchor_has_z:
            # we fix the anchor as should be done by the RGSS during the generation process
            t.z(anchor_left)

        # generate the right arms
        anchor_has_z = False
        t.h(anchor_right)
        for root in self.right_arms:
            anchor_has_z = anchor_has_z ^ helper_initialize_rgs_arm(t, root, anchor_right, outer_emitter, root_ancilla)
        if anchor_has_z:
            t.z(anchor_right)

        # join the two halves
        t.cz(anchor_left, anchor_right)
        t.h(anchor_left, anchor_right)
        meas_left = t.measure(anchor_left)
        meas_right = t.measure(anchor_right)

        # tracking the side effects (toggling first level nodes of all arms)
        if meas_left:
            for root in self.right_arms:
                for u in root.children:
                    u.has_z = not u.has_z
        if meas_right:
            for root in self.left_arms:
                for u in root.children:
                    u.has_z = not u.has_z

    def process_photon_loss(self, t: stim.TableauSimulator, loss_probability: float, rng: np.random.default_rng):
        for root in self.left_arms:
            helper_process_photon_loss(t, root, loss_probability, rng)
        for root in self.right_arms:
            helper_process_photon_loss(t, root, loss_probability, rng)

    def update_measurements_with_side_effect(self):
        for u in [*self.left_arms, *self.right_arms]:
            helper_update_eigenvalue_with_side_effect(u)

    def decode_logical_results(self) -> bool:
        """decoding the logical measurements of the inner qubits
        Returns False when the decoding fail and the trial needs retried"""
        for i, u in enumerate(self.left_arms):
            if i == self.successful_left_arm_index:
                self.left_logical_results[i] = decode_tree_logical_x(u)
            else:
                self.left_logical_results[i] = decode_tree_logical_z(u)

        for i, u in enumerate(self.right_arms):
            if i == self.successful_right_arm_index:
                self.right_logical_results[i] = decode_tree_logical_x(u)
            else:
                self.right_logical_results[i] = decode_tree_logical_z(u)

        return all(map(lambda res: res is not None, [*self.left_logical_results, *self.right_logical_results]))

    def count_lost_photons(self) -> tuple[int, int]:
        """return (lost_photons, total_photons)"""
        total_photons = 0
        lost_photons = 0
        for root in [*self.left_arms, *self.right_arms]:
            queue = [root]
            while len(queue) > 0:
                temp_queue = []
                for u in queue:
                    total_photons += 1
                    lost_photons += 1 if u.is_lost else 0
                    temp_queue.extend(u.children)
                queue = temp_queue
        return lost_photons, total_photons

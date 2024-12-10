#! usr/bin/python3

import random

import numpy as np
import stim

from node_qubit import Node, Pauli
from test_helper import verify_vertex_stabilizer
from tree_code_helper import decode_tree_logical_x, decode_tree_logical_z


def measurements_at_absa(t: stim.TableauSimulator, m: int, left_halfs: list[Node], right_halfs: list[Node]) -> int:
    """measurement of all qubits in the RGS (step 1) and return the index of the arm that has a successful BSM"""
    if m != len(left_halfs) or m != len(right_halfs):
        ValueError(f"number of arms {m} does not equal the input length of list of two halves {len(left_halfs)}, {len(right_halfs)}")

    # BSM part (outer qubit measurements)
    success_arm_index = -1
    for i in range(m):
        unode = left_halfs[i]
        vnode = right_halfs[i]

        if unode.is_lost or vnode.is_lost:
            continue

        u = unode.qubit_index
        v = vnode.qubit_index
        t.cz(u, v)
        t.h(u, v)
        unode.measurement_result = unode.eigenvalue = t.measure(u)
        vnode.measurement_result = vnode.eigenvalue = t.measure(v)
        unode.measurement_basis = vnode.measurement_basis = Pauli.X

        # simulating linear optics; consider +1/-1 and -1/+1 to be the two case ABSAs can distinguish
        if success_arm_index == -1 and (unode.measurement_result != vnode.measurement_result):
            success_arm_index = i

    # inner qubits measurements
    for i in range(m):
        unode = left_halfs[i]
        vnode = right_halfs[i]
        if i == success_arm_index:
            tree_code_physical_measure(t, unode, Pauli.X)
            tree_code_physical_measure(t, vnode, Pauli.X)
        else:
            tree_code_physical_measure(t, unode, Pauli.Z)
            tree_code_physical_measure(t, vnode, Pauli.Z)

    return success_arm_index


def update_tree_with_outer_qubits(left_tree_root: Node, right_tree_root: Node):
    """update in place with the BSM results; toggling 1st level results with root of another tree (step 2)"""
    if left_tree_root.is_lost or right_tree_root.is_lost:
        RuntimeError("trying to update tree with outer qubits that were lost!")

    if left_tree_root.eigenvalue:
        for u in right_tree_root.children:
            if u.is_lost:
                continue
            u.eigenvalue = not u.eigenvalue
    if right_tree_root.eigenvalue:
        for u in left_tree_root.children:
            if u.is_lost:
                continue
            u.eigenvalue = not u.eigenvalue


def compute_parity_for_end_nodes(m: int, left_logical_results: list[bool], right_logical_results: list[bool], successful_bsm_index: int) -> tuple[bool, bool]:
    """apply the parity at end nodes (step 3)
    return tuple of parity to be sent to the left and right respectively"""
    # parity multiplied together of Z of left (right) is sent to right (left),
    # while X of left (right) is sent to left (right).
    left_par, right_par = False, False
    for i in range(m):
        if i == successful_bsm_index:
            continue
        left_par ^= left_logical_results[i]
        right_par ^= right_logical_results[i]
    left_par ^= right_logical_results[successful_bsm_index]
    right_par ^= left_logical_results[successful_bsm_index]
    return left_par, right_par

total_photons = 0
lost_photons = 0

rng = np.random.default_rng()


def experiment_setup(
    number_of_hops: int,
    m: int,
    branching_parameters: list[int],
    loss_probability: float = 0,
    # photon_error_probability: float = 0,
    # emitter_error_probability: float = 0,
) -> tuple[bool, int, int]:
    """One run of the biclique RGS protocol
    Return: success-or-failure of the trial (bool), expectation value of ZX, expectation value of XZ at the end between Alice and Bob"""
    global total_photons, lost_photons
    # Ancilla qubits we require (total 4)
    #   temporary anchor for tree encoding: 1 (ancilla[0])
    #   emitter for outer qubit: 1 (ancilla[1])
    #   anchor for the half-RGS: 2 (ancilla[2-3])
    # qubit range
    alice = 0
    bob = 1
    # ancilla qubits
    anchor_left = 2
    anchor_right = 3
    outer_emitter = 4
    root_id = 5
    # starting index of unused qubit
    next_id = 6

    # we need (hop - 1) RGS
    rgss = [RGS(m, branching_parameters) for _ in range(number_of_hops - 1)]
    half_alice = HalfRGS(m, branching_parameters, alice)
    half_bob = HalfRGS(m, branching_parameters, bob)

    # assign qubit indices
    for rgs in rgss:
        next_id = rgs.assign_qubit_indices(next_id)
    next_id = half_alice.assign_qubit_indices(next_id)
    next_id = half_bob.assign_qubit_indices(next_id)

    # RGS creation
    t = stim.TableauSimulator()
    for rgs in rgss:
        rgs.initialize_quantum_state(t, anchor_left, anchor_right, outer_emitter, root_id)
    half_alice.initialize_quantum_state(t, outer_emitter, root_id)
    half_bob.initialize_quantum_state(t, outer_emitter, root_id)

    # process photon loss
    for rgs in rgss:
        rgs.process_photon_loss(t, loss_probability, rng)
    half_alice.process_photon_loss(t, loss_probability, rng)
    half_bob.process_photon_loss(t, loss_probability, rng)

    # Debugging, check how many photon got lost
    for rgs in rgss:
        lost_ph, total_ph = rgs.count_lost_photons()
        lost_photons += lost_ph
        total_photons += total_ph
    lost_ph, total_ph = half_alice.count_lost_photons()
    lost_photons += lost_ph
    total_photons += total_ph
    lost_ph, total_ph = half_bob.count_lost_photons()
    lost_photons += lost_ph
    total_photons += total_ph

    # (Protocol step 1) ABSA measurements
    success_bsm_indices = [-1] * number_of_hops  # number of ABSAs in the repeater chain
    for i in range(len(rgss) - 1):
        success_bsm_indices[i + 1] = measurements_at_absa(t, m, rgss[i].right_arms, rgss[i + 1].left_arms)
        rgss[i].successful_right_arm_index = rgss[i + 1].successful_left_arm_index = success_bsm_indices[i + 1]
    if len(rgss) > 0:
        success_bsm_indices[0] = measurements_at_absa(t, m, half_alice.arms, rgss[0].left_arms)
        success_bsm_indices[-1] = measurements_at_absa(t, m, rgss[-1].right_arms, half_bob.arms)
        rgss[0].successful_left_arm_index = half_alice.successful_arm_index = success_bsm_indices[0]
        rgss[-1].successful_right_arm_index = half_bob.successful_arm_index = success_bsm_indices[-1]
    else:
        # special case for 1 hop (no RGSS source nodes)
        success_bsm_indices[0] = measurements_at_absa(t, m, half_alice.arms, half_bob.arms)
        half_alice.successful_arm_index = half_bob.successful_arm_index = success_bsm_indices[0]

    # print(success_bsm_indices)
    if any(map(lambda id: id == -1, success_bsm_indices)):
        return False, None, None, None, None

    # (Protocol Step 1) Update measurements tree by assigning eigenvalues to the nodes taking side effects into account
    # imitating the classical messages received from RGSSs to ABSAs
    # we need to take note of the successful BSM arm index to denote the arm that undergone logical X measurements of inner qubits
    half_alice.update_measurement_with_side_effects()
    half_bob.update_measurement_with_side_effects()
    for i, rgs in enumerate(rgss):
        rgs.update_measurements_with_side_effect()

    # (Protocol Step 2) Propagating side effects of BSMs of outer qubits into their connected inner qubits
    bsm_arm_pairs = [half_alice.arms[success_bsm_indices[0]]]
    for i in range(len(rgss)):
        bsm_arm_pairs.append(rgss[i].left_arms[success_bsm_indices[i]])
        bsm_arm_pairs.append(rgss[i].right_arms[success_bsm_indices[i + 1]])
    bsm_arm_pairs.append(half_bob.arms[success_bsm_indices[-1]])
    for i in range(0, len(bsm_arm_pairs), 2):
        update_tree_with_outer_qubits(bsm_arm_pairs[i], bsm_arm_pairs[i + 1])

    # (Protocol Step 2/3?) Decoding logical measurements
    is_trial_successful = half_alice.decode_logical_results()
    is_trial_successful &= half_bob.decode_logical_results()
    for rgs in rgss:
        is_trial_successful &= rgs.decode_logical_results()
    if not is_trial_successful:
        return False, None, None, None, None

    # (Protocol Step 3) Compute parity at each ABSA for Pauli frame corrections
    parities: list[tuple[bool, bool]] = []
    if number_of_hops == 1:
        parities.append(compute_parity_for_end_nodes(m, half_alice.logical_results, half_bob.logical_results, success_bsm_indices[0]))
    else:
        parities.append(compute_parity_for_end_nodes(m, half_alice.logical_results, rgss[0].left_logical_results, success_bsm_indices[0]))
        for i in range(len(rgss) - 1):
            parities.append(compute_parity_for_end_nodes(m, rgss[i].right_logical_results, rgss[i + 1].left_logical_results, success_bsm_indices[i + 1]))
        parities.append(compute_parity_for_end_nodes(m, rgss[-1].right_logical_results, half_bob.logical_results, success_bsm_indices[-1]))

    # (Protocol Step 4) Combining all the parities from all ABSAs and correct at end nodes
    total_parity = (False, False)
    for l, r in parities:
        total_parity = total_parity[0] ^ l, total_parity[1] ^ r
    if total_parity[0]:
        t.z(alice)
    if total_parity[1]:
        t.z(bob)

    return True, t.peek_observable_expectation(stim.PauliString("XZ")), t.peek_observable_expectation(stim.PauliString("ZX")), t.canonical_stabilizers(), total_parity
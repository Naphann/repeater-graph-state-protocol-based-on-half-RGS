#! usr/bin/python3

from helper import Pauli, Node
import stim


def tree_code_physical_measure(t: stim.TableauSimulator, root: Node, is_starting_with_x: bool):
    queue = root.children
    is_x_basis = is_starting_with_x
    while len(queue) > 0:
        temp_queue = []
        for u in queue:
            if is_x_basis:
                t.h(u.qubit_index)
            if not u.is_lost:
                # if the qubit is lost, don't record or assign measurement result
                u.eigenvalue = u.measurement_result = t.measure(u.qubit_index)
                u.measurement_basis = Pauli.X if is_x_basis else Pauli.Z
            temp_queue.extend(u.children)
        is_x_basis = not is_x_basis
        queue = temp_queue


def __get_z_result(root: Node) -> bool | None:
    # TODO: implement majority vote mechanism
    # direct result available
    if not root.is_lost:
        return root.eigenvalue
    # indirect measurement
    for u in root.children:
        if u.is_lost:
            continue
        next_level_zs = [__get_z_result(v) for v in u.children]
        if not all([z is not None for z in next_level_zs]):
            continue
        parity = u.eigenvalue
        [parity := parity ^ z for z in next_level_zs]
        return parity
    # measurement result cannot be determined neither direct nor indirect
    return None


def decode_tree_logical_z(root: Node) -> bool | None:
    """find the XOR (parity) of the first level nodes"""
    zs = [__get_z_result(u) for u in root.children]
    if any(map(lambda z: z is None, zs)):
        return None
    logical_z = False
    [logical_z := logical_z ^ z for z in zs]
    return logical_z


def decode_tree_logical_x(root: Node) -> bool | None:
    """find successful X in the first level and return the parity of it with all its children's Z results
    i.e., the parity of the X_i Z_children_of_i of any first level node i"""
    # TODO: implement majority vote mechanism
    for u in root.children:
        if u.is_lost:
            continue
        zs = [__get_z_result(v) for v in u.children]
        if len(zs) > 0 and any(map(lambda z: z is None, zs)):
            continue
        logical_x = u.eigenvalue
        [logical_x := logical_x ^ z for z in zs]
        return logical_x
    return None

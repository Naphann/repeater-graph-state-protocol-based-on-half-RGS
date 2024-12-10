#! usr/bin/python3

from rgs_config import Node
import numpy as np


def __get_z_result(root: Node) -> bool | None:
    """Get measurement result of this node qubit.
    Preferentially select the majority vote of indirect measurements is available."""
    # indirect measurement
    indirect_zs: list[bool] = []
    for u in root.children:
        if u.is_lost:
            continue
        next_level_zs = [__get_z_result(v) for v in u.children]
        if not all([z is not None for z in next_level_zs]):
            continue
        parity = u.eigenvalue
        [parity := parity ^ z for z in next_level_zs]  # type: ignore
        indirect_zs.append(parity)
    # majority vote for the result guessing
    indirect_meas_count = len(indirect_zs)
    minus_count = np.sum(indirect_zs)
    plus_count = indirect_meas_count - minus_count

    if indirect_meas_count == 0:
        # fall back to direct measurement result if available
        return root.eigenvalue if not root.is_lost else None

    if plus_count == minus_count:
        # TODO: can we do better than randomly select? Choose one with least confidence to remove?
        #       i.e., the one where we get results from direct measurements down below
        if not root.is_lost:
            return root.eigenvalue
        else:
            print("Z: randomly chosen")
            return np.random.choice([True, False])
    else:
        return minus_count > plus_count


def decode_tree_logical_z(root: Node) -> bool | None:
    """find the XOR (parity) of the first level nodes"""
    zs = [__get_z_result(u) for u in root.children]
    if len(zs) == 0:
        raise RuntimeError("We should not encounter this at all!")
    if any(map(lambda z: z is None, zs)):
        return None
    logical_z = False
    [logical_z := logical_z ^ z for z in zs]  # type: ignore
    return logical_z


def decode_tree_logical_x(root: Node) -> bool | None:
    """find successful X in the first level and return the parity of it with all its children's Z results
    i.e., the parity of the X_i Z_children_of_i of any first level node i"""
    results: list[bool] = []
    for u in root.children:
        if u.is_lost:
            continue
        zs = [__get_z_result(v) for v in u.children]
        if len(zs) > 0 and any(map(lambda z: z is None, zs)):
            continue
        logical_x = u.eigenvalue
        [logical_x := logical_x ^ z for z in zs]  # type: ignore
        results.append(logical_x)

    success_count = len(results)
    minus_count = np.sum(results)
    plus_count = success_count - minus_count

    if success_count == 0:
        return None

    if plus_count == minus_count:
        # TODO: better way of randomly choosing the result?
        print("X: randomly chosen")
        return np.random.choice([True, False])
    else:
        return minus_count > plus_count

#! usr/bin/python3

import numpy as np
import stim


def verify_vertex_stabilizer(t: stim.TableauSimulator, vertex: int, neighbours: list[int], expected_value) -> bool:
    r = max(np.max(neighbours), vertex)
    stabilizer = ["_"] * (r + 1)
    stabilizer[vertex] = "X"
    for v in neighbours:
        stabilizer[v] = "Z"
    if t.peek_observable_expectation(stim.PauliString("".join(stabilizer))) != expected_value:
        error_message = f"stabilizer not correct: {vertex}: {neighbours}\n    Given expected value = {expected_value} but got {t.peek_observable_expectation(stim.PauliString(''.join(stabilizer)))}"
        raise RuntimeError(error_message)

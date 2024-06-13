import numpy as np


def prob_bell(m, p):
    # linear BSM 50% when photons arrive
    return 1 - (1 - p * p / 2) ** m


def prob_logical_measure_z(bv, p):
    n = len(bv)
    if n % 2 == 1:
        mz = False
    else:
        mz = True
    prob = p
    for b in reversed(bv[1:]):
        if mz:
            prob = p + (1 - p) * (1 - (1 - prob) ** b)
        else:
            prob = p * prob**b
        mz = not mz
    prob = prob ** bv[0]
    return prob


def prob_logical_measure_x(bv, p):
    n = len(bv)
    if n % 2 == 1:
        mz = True
    else:
        mz = False
    prob = p
    for b in reversed(bv[1:]):
        if mz:
            prob = p + (1 - p) * (1 - (1 - prob) ** b)
        else:
            prob = p * prob**b
        mz = not mz
    prob = 1 - (1 - prob) ** bv[0]
    return prob


def photon_arrival_probability_from_km_distance(distance, loss_db_per_km=0.2):
    attenuation_distance = 10 / (np.log(10) * loss_db_per_km)
    return np.exp(-distance / attenuation_distance)


def prob_rgs_trial(m, bv, p, number_of_hops):
    """p is photon arrival probability"""
    p_one_hop = prob_bell(m, p) * (prob_logical_measure_x(bv, p) ** (2)) * (prob_logical_measure_z(bv, p) ** (2 * m - 2))
    return p_one_hop**number_of_hops


# def trial_success_probability(number_of_hops: int, m: int, bvec: list[int], loss_probability: float) -> float:
#     p_bell = 0.5
#     p_ph = 1 - loss_probability

#     def s(k):
#         return p_ph * p

#     def r(k):
#         return 1 - (1 - s(k)) ** bvec[k]


#     def prob_logical_x():
#         nonlocal number_of_hops, m, bvec, loss_probability
#         pass

#     def prob_logical_z():
#         nonlocal number_of_hops, m, bvec, loss_probability

#     def p_single_hop():
#         nonlocal number_of_hops, m, bvec, loss_probability, p_bell
#         (1 - (1 - p_bell) ** m) * prob_logical_x() ** 2 * prob_logical_z() ** (2*m - 2)

#     return p_single_hop() ** number_of_hops

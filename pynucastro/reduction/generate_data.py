#!/usr/bin/env python3

import numpy as np

from pynucastro import Composition, Nucleus
from pynucastro.reduction.load_network import load_network


def dataset(network, n=10, permute=True, b_rho=None, b_T=None, b_Z=None):
    """
    Generate a dataset with *n* datapoints. Will either be returned as a sequence of tuples,
    each with order (composition, density, temperature) if *permute* is *True* (default),
    or the transpose of that if *permute* is *False*. The parameters *b_rho*, *b_T*, and
    *b_Z* are tuples giving the bounds on density, temperature, and metallicity respectively.
    """

    if isinstance(n, int):
        n = np.ones(3, dtype=np.int32) * n
    else:
        n = np.array(list(n), dtype=np.int32)

    # Bounds on each variable
    if b_rho is None:
        b_rho = (1e2, 1e6)  # density (g/cm^3)
    if b_T is None:
        b_T = (8e6, 1.5e9)  # temperature (K)
    if b_Z is None:
        b_Z = (0.02, 0.2)  # metallicity

    rho = np.geomspace(*b_rho, num=n[1])
    T = np.geomspace(*b_T, num=n[2])
    comp_list = []

    for Z in np.linspace(*b_Z, num=n[0]):

        comp = Composition(network.get_nuclei())
        # 75/25 H1/He4 ratio
        comp.X[Nucleus("p")] = (1. - Z) * 0.75
        comp.X[Nucleus("he4")] = (1. - Z) * 0.25
        # 50% of remaining stuff in carbon-12
        comp.X[Nucleus("c12")] = Z * 0.5
        # Split the rest across all other nuclei
        omit = {Nucleus("p"), Nucleus("he4"), Nucleus("c12")}
        rem = 0.5 * Z / (len(comp.X) - len(omit))

        for nuc in comp.X:
            if nuc not in omit:
                comp.X[nuc] = rem

        comp.normalize()
        comp_list.append(comp)

    if not permute:
        yield comp_list
        yield rho
        yield T
        return

    for rho_i in rho:
        for T_i in T:
            for comp in comp_list:
                yield (comp, rho_i, T_i)


def main():
    network = load_network(Nucleus("ni56"))
    conds_list = list(dataset(network))
    rho = sorted({conds[1] for conds in conds_list})
    T = sorted({conds[2] for conds in conds_list})
    comp = {conds[0] for conds in conds_list}

    print("ρ")
    print(rho)
    print()

    print("T")
    print(T)
    print()

    def comp_converter(comp):

        X = comp.X[Nucleus("p")]
        Y = comp.X[Nucleus("he4")]
        ZC12 = comp.X[Nucleus("c12")]
        Zmin = min(comp.X.values())
        return (X, Y, ZC12, Zmin)

    print("X")
    print(np.array(list(map(comp_converter, comp))))


if __name__ == "__main__":
    main()

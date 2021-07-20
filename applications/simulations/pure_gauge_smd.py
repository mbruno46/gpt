#!/usr/bin/env python3
#
# Authors: Christoph Lehner, Mattia Bruno 2021
#
# HMC for phi^4 scalar theory
#
import gpt as g
import sys, os
import numpy

beta = g.default.get_float("--beta", 5.96)

g.default.set_verbose("omf4")

grid = g.grid([8, 8, 8, 16], g.double)
rng = g.random("smd-pure-gauge")

U = g.qcd.gauge.unit(grid)
rng.normal_element(U)

# conjugate momenta
mom = g.group.cartesian(U)
rng.normal_element(mom)

# Log
g.message(f"Lattice = {grid.fdimensions}")
g.message("Actions:")
# action for conj. momenta
a0 = g.qcd.scalar.action.mass_term()
g.message(f" - {a0.__name__}")

# wilson action
a1 = g.qcd.gauge.action.wilson(beta)
g.message(f" - {a1.__name__}")


def hamiltonian():
    return a0(mom) + a1(U)


# molecular dynamics
sympl = g.algorithms.integrator.symplectic

ip = sympl.update_p(mom, lambda: a1.gradient(U, U))
iq = sympl.update_q(U, lambda: a0.gradient(mom, mom))

# integrator
mdint = sympl.OMF4(1, ip, iq)
g.message(f"Integration scheme:\n{mdint}")

# metropolis
def on_reject(fields, previous_fields):
    for i in range(len(fields)):
        fields[i] @= (1 if i < 4 else -1) * previous_fields[i]


metro = g.algorithms.markov.metropolis(rng, on_reject)

g.message("SMD smd")
eps = 0.2
gamma = 0.3
c1 = numpy.exp(-gamma * eps)
c2 = (1 - c1 * c1) ** 0.5
g.message(f"eps = {eps}; gamma = {gamma}")


def momentum_rotation(mom):
    nu = g.copy(mom)
    rng.normal_element(nu)
    for m, n in zip(mom, nu):
        m @= c1 * m + c2 * n


def smd(eps, mom):
    accrej = metro(U + mom)
    momentum_rotation(mom)
    h0 = hamiltonian()
    mdint(eps)
    h1 = hamiltonian()
    return [accrej(h1, h0), h1 - h0]


# thermalization
ntherm = 100
for i in range(1, 11):
    h = []
    timer = g.timer("hmc")
    for _ in range(ntherm // 10):
        timer("trajectory")
        h += [smd(eps, mom)]
    h = numpy.array(h)
    timer()
    g.message(f"{i*10} % of thermalization completed")
    g.message(
        f'Average time per trajectory = {timer.dt["trajectory"]/ntherm*10:g} secs'
    )
    g.message(
        f"Plaquette = {g.qcd.gauge.plaquette(U)}, Acceptance = {numpy.mean(h[:,0]):.2f}, |dH| = {numpy.mean(numpy.abs(h[:,1])):.4e}"
    )

# production
history = []
data = []
n = 100
dtrj = 10
for i in range(n):
    for k in range(dtrj):
        history += [smd(eps, mom)]
    data += [g.qcd.gauge.plaquette(U)]
    g.message(f"Trajectory {i}, {history[-1]}")

history = numpy.array(history)
g.message(f"Acceptance rate = {numpy.mean(history[:,0]):.2f}")
g.message(f"<|dH|> = {numpy.mean(numpy.abs(history[:,1])):.4e}")

data = numpy.array(data)
g.message(f"<plaq>   = {numpy.mean(data[:,0])}")

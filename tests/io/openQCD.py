#!/usr/bin/env python3
#
# Authors: Christoph Lehner 2020
#
# Desc.: Illustrate core concepts and features
#
import gpt as g
import numpy as np
import sys
import time

# load configuration
U = g.load("/hpcgpfs01/work/clehner/configs/openQCD/A250t000n54")

# do everything in single-precision
U = g.convert(U, g.single)

# use the gauge configuration grid
grid=U[0].grid

# quark
w=g.qcd.fermion.wilson_clover(U,{
    "kappa" : 0.13565,
    "csw_r" : 2.0171,
    "csw_t" : 2.0171,
    "xi_0" : 1,
    "nu" : 1,
    "isAnisotropic" : False,
    "boundary_phases" : [ 1.0, 1.0, 1.0, 1.0 ]
})

# create point source
src=g.mspincolor(grid)
g.create.point(src, [0,0,0,0])

# build solver
s=g.qcd.fermion.solver
cg=g.algorithms.iterative.cg({
    "eps" : 1e-6,
    "maxiter" : 1000
})
slv=s.propagator(s.eo_ne(g.qcd.fermion.preconditioner.eo2(w), cg))

# propagator
dst=g.mspincolor(grid)
slv(src,dst)

# two-point
correlator=g.slice(g.trace(dst*g.adj(dst)),3)

# output
for t,c in enumerate(correlator):
    g.message(t,c.real)

#
#    GPT - Grid Python Toolkit
#    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Note: we use the proper order of the chebyshev_t
#       in contrast to current Grid
#
import gpt as g
import numpy as np
from gpt.params import params_convention


class coarse_deflate:
    @params_convention(block=32)
    def __init__(self, inverter, cevec, basis, fev, params):
        self.inverter = inverter
        self.basis = basis
        self.params = params

        def noop(matrix):
            def noop_mat(dst, src):
                pass

            return noop_mat

        self.cdefl = g.algorithms.eigen.deflate(noop, cevec, fev, params)
        self.cgrid = cevec[0].grid

    def __call__(self, matrix):

        otype = self.basis[0].otype
        grid = self.basis[0].grid
        cb = self.basis[0].checkerboard()

        cdefl = self.cdefl(matrix)
        inverter = self.inverter(matrix)
        b = g.block.map(self.cgrid, self.basis)

        def inv(dst, src):
            verbose = g.default.is_verbose("deflate")
            # |dst> = sum_n 1/ev[n] |n><n|src>

            # temporaries
            template = self.cdefl.evec[0]
            csrc = [g.lattice(template) for x in src]
            cdst = [g.lattice(template) for x in src]

            t0 = g.time()
            b.project(csrc, src)
            t1 = g.time()
            # g.default.push_verbose("deflate", False)
            cdefl(cdst, csrc)
            # g.default.pop_verbose()
            t2 = g.time()
            b.promote(dst, cdst)
            t3 = g.time()
            if verbose:
                g.message(
                    "Coarse-grid deflated %d vector(s) in %g s (project %g s, coarse deflate %g s, promote %g s)"
                    % (len(src), t3 - t0, t1 - t0, t2 - t1, t3 - t2)
                )
            return inverter(dst, src)

        return g.matrix_operator(
            mat=inv, inv_mat=matrix, otype=otype, grid=grid, cb=cb, accept_list=True
        )

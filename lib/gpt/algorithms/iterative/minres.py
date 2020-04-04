#
#    GPT - Grid Python Toolkit
#    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)
#                        Mattia Bruno
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
# MinRes algorithm: solution to M x = b
#
# x=0, p=b, r=b
# for k in 0,1,...
#   q = Ar
#   a = (q, r)/|q|^2
#   x <- x + a*p
#   r <- r - a*q
#
import gpt as g
import numpy as np
from time import time

class minres:
    
    def __init__(self, params):
        self.params = params
        self.nmr = params["nmr"]
        self.res = params["res"]
        
    def __call__(self, matrix, eta, psi):
        verbose=g.default.is_verbose("minres")
        t0=time()
        
        q=g.lattice(eta)
        x=g.lattice(eta)
        x[:]=0
        r=g.copy(eta)
        psi[:]=0
        
        rn=np.sqrt(g.norm2(eta))
        tol=rn*self.res
        
        for i in range(self.nmr):
            matrix(r,q)
            rn=g.norm2(q)
            if (rn==0.0):
                continue
            a = g.innerProduct(q,r)*(1.0/rn)
            x @= x + a*r
            r @= r - a*q
            
            rn=np.sqrt(g.norm2(r))
            if verbose:
                g.message("minres: %d, |rho| = %g, |psi| = %g" % (i,rn,np.sqrt(g.norm2(x))))
            if (rn<=tol):
                break
                
        if verbose:
            t1=time()
            if (rn<=tol):
                g.message("minres: converged in %g sec " % (t1-t0))
            else:
                g.message("minres: did not converge in %g sec " % (t1-t0))
                
        psi @= x
        return r
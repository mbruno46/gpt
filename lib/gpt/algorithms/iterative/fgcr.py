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
# Pseudo code and algorithm from 
#    M. Luescher, Solution of the Dirac equation in lattice QCD 
#                 using a domain decomposition method,
#                 https://arxiv.org/pdf/hep-lat/0310048.pdf
# 
# input arguments:
#   preconditioner M = approx. solution single prec
#   operator D = Dirac operator. in double or single prec
#   eta, psi = source, solution
#
# while not converged
#   GCR CYCLE
#   rho = eta
#   for k=0...
#     phi[k] = M rho[k]
#     chi[k] = D phi[k]
#     for l=0...k-1
#       a[l,k] = (chi[l], chi[k])
#       chi[k] = chi[k] - a[l,k]*chi[l]
#     b[k] = |chi[k]|
#     chi[k]=chi[k]/b[k]
#     c[k] = (chi[k],rho[k])
#     rho[k+1] = rho[k] - c[k]*chi[k]
#
# UPDATE SOLUTION
# psi = sum_l alpha[l] * phi[l]
# rho = eta - D psi
#
# b_l alpha_l + sum_{i=l+1}^k a[l,i]*alpha[i] = c[l] for l=k,k-1,...,0
#
import gpt as g
import numpy as np
from time import time

class fgcr:
    
    def __init__(self, params):
        self.params = params
        #self.eps = params["eps"]
        self.maxiter = params["maxiter"]
        self.nkv = params["nkv"]
        self.res = params["res"]
        self.gcr_prec = 1e-6
   
    def get_alpha(self,k, a, b, c):
        alpha=[0.]*(k+1)
        alpha[k]=c[k]/b[k]
        for l in range(k-1,-1,-1):
            alpha[l]= (c[l] - np.dot(a[l,l+1:k+1],alpha[l+1:k+1]))/b[l]
        return alpha
   
    def __call__(self, M, D, eta, psi):
        verbose=g.default.is_verbose("fgcr")
        t0=time()
        
        a=np.zeros((self.nkv,self.nkv),dtype=np.complex)
        b=[0.]*self.nkv #np.zeros((self.nkv,))
        c=[0.]*self.nkv #np.zeros((self.nkv,))
        
        rho=g.vspincolor(M.F_grid)
        rho=g.convert(rho, eta)
        psi[:]=0
        wsd=g.lattice(psi)
        chi=[g.lattice(rho)]*self.nkv
        phi=[g.lattice(rho)]*self.nkv

        rn=np.sqrt(g.norm2(eta))
        tol=rn*self.res
        it=0
        while True:
            rn0=rn
            
            for k in range(self.nkv):
                # gcr step
                chi[k]@=M(rho,phi[k])
                
                for l in range(k):
                    a[l,k] = g.innerProduct(chi[l],chi[k])
                    chi[k] @= chi[k] - a[l,k]*chi[l];
                
                b[k]=np.sqrt(g.norm2(chi[k]))
                chi[k] *= (1.0/b[k])
                c[k] = g.innerProduct(chi[k],rho)
                rho -= c[k]*chi[k]
                # end gcr step
                
                it+=1
                rn=np.sqrt(g.norm2(rho))
                if (it>self.maxiter) or (rn<=tol):
                    break
                if (rn<rn0*self.gcr_prec):
                    break
            
            # update psi
            alpha=self.get_alpha(k,a,b,c)
            rho[:]=0
            for l in range(k,-1,-1):
                rho += phi[l]*alpha[l]
            g.convert(wsd, rho)
            psi @= psi + wsd
            D(psi,wsd)
            wsd @= wsd-eta
            rho = g.convert(rho, wsd)
            
            rn=np.sqrt(g.norm2(rho))
            if verbose:
                g.message("fgcr: %d, |rho| = %g, |psi| = %g" % (it,rn,np.sqrt(g.norm2(psi))))
            if (rn<=tol):
                if verbose:
                    t1=time()
                    g.message("fgcr converged in %g sec " % (t1-t0))
                break
            else:
                if (it>self.maxiter):
                    if verbose:
                        t1=time()
                        g.message("fgcr not converged in %g sec " % (t1-t0))
                    break

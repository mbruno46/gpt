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
import gpt as g
import cgpt
import numpy as np
import time

# ix = x[3] + lat[3]*(x[2] + lat[2]*x[1]+..))
def index2coor(ix,lat):
    x=[0]*len(lat)
    for mu in range(len(lat)-1,-1,-1):
        x[mu] = int(ix % lat[mu])
        ix = np.floor(ix/lat[mu])
    return x

class sap_blocks:
    
    def __init__(self,U,bs,ns):
        Nd = len(bs)
        self.blk_size = np.array(bs)
        self.grid = U[0].grid
        if np.any(self.blk_size>self.grid.ldimensions):
            raise
        self.blk_dims = np.floor_divide(self.grid.ldimensions,self.blk_size)
        self.nblocks = int(np.prod(self.blk_dims))
        if (self.nblocks<2) or ((self.nblocks % 2)!=0):
            raise
        self.blk_vol = int(np.prod(self.blk_size))
        if (self.blk_vol<4):
            raise
        
        # extended block coor system
        exbs = [bs[0]*int(self.nblocks/2)] + bs[1:]
        self.blk_grid = g.grid(exbs,g.double) 
        self.Ublk = [[g.mcolor(self.blk_grid)]*4]*2
        for mu in range(4):
            self.Ublk[0][mu][:]=0
            self.Ublk[1][mu][:]=0
        
        self.nmodes = ns
        self.sfld = [[g.vspincolor(self.blk_grid)]*self.nmodes]*2
        for i in range(self.nmodes):
            self.sfld[0][i][:]=0
            self.sfld[1][i][:]=0
            
        self.lcoor = []
        self.bcoor = []
        self.beo = []
        idx=[-1,-1]
        for ib in range(self.nblocks):
            blk_coor = index2coor(ib,self.blk_dims) # coor of block in block-layout
            blk_ofs = self.blk_size*blk_coor   # origin of block in full-grid
            
            eo = int(np.sum(blk_coor) % 2) # even-odd block
            self.beo.append( eo )
            idx[eo]+=1
            exofs = np.array([bs[0]*idx[eo]] + bs[1:]) # ofs in extended block-coord system
            
            self.lcoor.append(self.get_coor(blk_ofs,blk_ofs+self.blk_size))
            self.bcoor.append(self.get_coor(exbs, exbs+self.blk_size))
            
            lcoor_obc=self.get_coor(blk_ofs,[blk_ofs[mu]+self.blk_size[mu]-1 for mu in range(Nd)])
            bcoor_obc=self.get_coor(exofs,[exofs[mu]+self.blk_size[mu]-1 for mu in range(Nd)])
            
            for mu in range(4):
                self.Ublk[eo][mu][bcoor_obc] = U[mu][lcoor_obc]
    
    def get_coor(self,top,bottom):
        return cgpt.coordinates_form_cartesian_view([int(t) for t in top],[int(b) for b in bottom])
    
    def lat2blk(self, src, eo, imode):
        for ib in range(self.nblocks):
            if self.beo[ib]==eo:
                self.sfld[eo][imode][self.bcoor[ib]] = src[self.lcoor[ib]]



class sap:
    
    def __init__(self, params, U):
        self.blk_size = params["blk_size"]
        self.nmr = params["nmr"]
        self.ncy = params["ncy"]
        self.mres = g.algorithms.iterative.minres({
            "res" : 1e-6,
            "nmr" : self.nmr
        })
        
        self.blk = sap_blocks(U,self.blk_size,3)
        self.D_blk = []
        
        #tmp_params = {}
        #for k in D.params:
        #    assert(not k in [ "U_grid", "U_grid_rb", "F_grid", "F_grid_rb", "U" ])
        #    tmp_params[k] = D.params[k]
            
        #if "wilson" in D.name:
        #    tmp_params["csw_r"],tmp_params["csw_t"] = 0.,0.
        #    for eo in [0,1]:
        #        self.D.append( g.qcd.fermion.wilson_exp_clover(self.blk.Ublk[eo], tmp_params) )
        #        #cgpt.update_block_clover
        #else:
        #    for eo in [0,1]:
        #        if D.name=="mobius":
        #            self.D.append( g.qcd.fermion.mobius(self.blk.Ublk[eo], tmp_params) )
        #        elif D.name=="zmobius":
        #            self.D.append( g.qcd.fermion.zmobius(self.blk.Ublk[eo], tmp_params) )

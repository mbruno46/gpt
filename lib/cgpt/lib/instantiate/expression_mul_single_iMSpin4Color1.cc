// This file is automatically generated, do not modify!
/*
    GPT - Grid Python Toolkit
    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#include "../lib.h"

#include "../expression/matmul.h"
#include "../expression/mul.h"

template<>
cgpt_Lattice_base* cgpt_lattice_gammamul(cgpt_Lattice_base* dst, bool ac, int unary_a, Lattice< iMSpin4Color1<vComplexF> >& la, Gamma::Algebra gamma, int unary_expr, bool rev) {
  if (rev) {
    return lattice_unary_rmul(dst, ac, unary_a, la, Gamma(gamma), unary_expr);
  } else {
    return lattice_unary_mul(dst, ac, unary_a, la, Gamma(gamma), unary_expr);
  }
}

template<>
cgpt_Lattice_base* cgpt_lattice_matmul(cgpt_Lattice_base* dst, bool ac, int unary_a, Lattice< iMSpin4Color1<vComplexF> >& la, PyArrayObject* b, std::string& bot, int unary_b, int unary_expr, bool rev) {
  typedef vComplexF vtype;
  if (unary_b == 0) {
    _MM_COMPATIBLE_RL_(iMSpin4);
    _MM_COMPATIBLE_RL_(iMColor1);
    _MM_COMPATIBLE_RL_(iMSpin4Color1);
    _MM_COMPATIBLE_L_(iVSpin4Color1);
  }
  ERR("Not implemented");
}

template<>
cgpt_Lattice_base* cgpt_lattice_mul(cgpt_Lattice_base* dst, bool ac, int unary_a, Lattice< iMSpin4Color1<vComplexF> >& la,int unary_b, cgpt_Lattice_base* b, int unary_expr) {
  typedef vComplexF vtype;
  _COMPATIBLE_(iSinglet);
  _COMPATIBLE_(iMSpin4);
  _COMPATIBLE_(iMColor1);
  _COMPATIBLE_(iVSpin4Color1);
  _COMPATIBLE_(iMSpin4Color1);
  ERR("Not implemented");
}

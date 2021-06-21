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
import gpt as g


def line_search_quadratic(s, x, dv0, df, step):
    dv1 = df(g(g.group.compose(step * s, x)))
    # ansatz: f(x) = a + b*(x-c)^2, then solve for c from dv1 and dv0
    sv0 = dv0.otype.inner_product(s, dv0)
    sv1 = dv0.otype.inner_product(s, dv1)
    r = sv0 / sv1
    if abs(r - 1.0) < 1e-15:
        return 1.0
    c = r / (r - 1.0)
    return c


def line_search_none(s, x, dv0, df, step):
    return 1.0
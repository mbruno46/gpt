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


def projected_convert(x, otype):
    return g.project(g.convert(x, otype), "defect")


def compose(left, right):
    as_list = isinstance(left, list)

    left = g.util.to_list(g(left))
    right = g.util.to_list(g(right))

    left_type = left[0].otype
    right_type = right[0].otype

    if left_type.__name__ == right_type.__name__:
        # if both are of same type, use common group compose and return
        dst = [left_type.compose(l, r) for l, r in zip(left, right)]

    else:
        # if they are not, see if either is a cartesian element of the other type
        left_type_cartesian = left_type.cartesian()
        right_type_cartesian = right_type.cartesian()

        if left_type.__name__ == right_type_cartesian.__name__:
            dst = [
                right_type.compose(projected_convert(l, right_type), r)
                for l, r in zip(left, right)
            ]

        elif left_type_cartesian.__name__ == right_type.__name__:
            dst = [
                left_type.compose(l, projected_convert(r, left_type))
                for l, r in zip(left, right)
            ]

        else:
            raise TypeError(
                f"{left_type.__name__} and {right_type.__name__} are not composable"
            )

    if as_list:
        return dst

    return dst[0]


def approximate_gradient(x, functional, *site, epsilon=1e-7):
    # This helper function allows for quick checks of gradient implementations on single sites
    c = x.otype.cartesian()
    grid = x.grid
    epsilon = complex(epsilon)

    # move to neutral element of group (\vec{0} in cartesian space)
    t = g.lattice(grid, c)
    t[:] = 0
    r = t[site]

    # generators of cartesian space
    gen = c.generators(grid.precision.complex_dtype)

    # functional at neutral element
    f0 = functional(x)

    for i, gg in enumerate(gen):
        t[site] += epsilon * gg
        r += ((functional(g(g.group.compose(t, x))) - f0) / epsilon) * gg
        t[site] -= epsilon * gg

    return r
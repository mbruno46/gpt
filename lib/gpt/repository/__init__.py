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
from .repo import repo


_repository = {'Qlattice': repo('https://github.com/waterret/Qlattice','master')}
_repository['Qlattice'].addfile('examples/propagators/sample-results/test-4nt8/results=1000/psrc-prop-0.field')
_repository['Qlattice'].addfile('examples/propagators/sample-results/test-4nt8/results=1000/pion-corr.txt')

def repository(key):
    if key in _repository:
        return _repository[key]

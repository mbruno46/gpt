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
#  Reader for the file format of the Qlattice library written by Luchang Jin
#    https://github.com/waterret/Qlattice
#
#  For file format documentation see
#    https://github.com/waterret/Qlattice/blob/master/docs/file-format.md
#
import cgpt, gpt, numpy, os, sys


class qlat_io:
    def __init__(self, path):
        self.path = path
        self.ldimensions = []
        self.fdimensions = []
        self.bytes_header = -1
        self.verbose = gpt.default.is_verbose("io")
        self.size = 1
        gpt.barrier()

    def read_header(self):

        # make sure this is a file
        if not os.path.isfile(self.path):
            return False

        with open(self.path, "rb") as f:
            line = self.getline(f)
            if line != "BEGIN_FIELD_HEADER":
                return False

            # need to be mute before this line since this is used to autodetect the file format
            if self.verbose:
                gpt.message(f"Qlattice file format; reading {self.path}")
                gpt.message(f"   {line}")

            for i in range(8):
                line = self.getline(f)
                if self.verbose:
                    gpt.message(f"\t{line}")

                [field, val] = line.split(" = ")
                if field == "field_version":
                    self.version = val
                    assert val == "1.0"
                elif field[:10] == "total_site":
                    self.fdimensions.append(int(val))
                elif field == "multiplicity":
                    if int(val) == 288:
                        self.otype = gpt.ot_matrix_spin_color(4, 3)
                    self.size *= int(val)
                elif field == "sizeof(M)":
                    if int(val) == 4:
                        self.precision = gpt.single
                    elif int(val) == 8:
                        self.precision = gpt.double
                    self.size *= int(val)
                elif field == "field_crc32":
                    self.crc_exp = val

            line = self.getline(f)
            if self.verbose:
                gpt.message(f"   {line}")
            assert line == "END_HEADER"
            self.bytes_header = f.tell()

        self.cv_desc = "[" + ",".join(["1"] * len(self.fdimensions)) + "]"
        self.ldimensions = [fd for fd in self.fdimensions]
        return True

    def getline(self, f):
        line = []
        while True:
            try:
                c = str(f.read(1), "utf-8")
            except UnicodeDecodeError:
                # this should not happen and likely indicates that we read an incompatible file
                break
            if c == "\n":
                break
            else:
                line += [c]
            if len(line) > 1024:
                break
        return "".join(line)

    def swap(self, data):
        if sys.byteorder == "big":
            return

        if self.precision == gpt.single:
            size = 4
        elif self.precision == gpt.double:
            size = 8
        cgpt.munge_byte_order(data, data, size)

    def read_lattice_single(self):
        if self.bytes_header < 0:
            raise

        # define grid from header
        g = gpt.grid(self.fdimensions, self.precision)
        # create lattice
        l = gpt.lattice(g, self.otype)

        # performance
        dt_distr, dt_crc, dt_read, dt_misc = 0.0, 0.0, 0.0, 0.0
        szGB = 0.0
        crc_comp = 0
        g.barrier()
        t0 = gpt.time()

        # single file: each rank opens it and reads it all
        g.barrier()
        dt_read -= gpt.time()

        cv = gpt.cartesian_view(
            gpt.rank(), self.cv_desc, g.fdimensions, g.cb, l.checkerboard()
        )
        pos = gpt.coordinates(cv)

        if gpt.rank() == 0:
            f = gpt.FILE(self.path, "rb")
            f.seek(self.bytes_header, 0)
            sz = self.size * int(numpy.prod(g.fdimensions))
            data = memoryview(bytearray(f.read(sz)))
            f.close()

            dt_crc -= gpt.time()
            crc_comp = gpt.crc32(data)
            crc_comp = f"{crc_comp:8X}"
            assert crc_comp == self.crc_exp
            dt_crc += gpt.time()

            dt_misc -= gpt.time()
            self.swap(data)
            dt_misc += gpt.time()

            sys.stdout.flush()
            szGB += len(data) / 1024.0 ** 3.0
        else:
            assert len(pos) == 0
            data = None

        g.barrier()
        dt_read += gpt.time()

        # distributes data accordingly
        dt_distr -= gpt.time()
        l[pos] = data
        g.barrier()
        dt_distr += gpt.time()

        g.barrier()
        t1 = gpt.time()

        szGB = g.globalsum(szGB)
        if self.verbose and dt_crc != 0.0:
            gpt.message(
                "Read %g GB at %g GB/s (%g GB/s for distribution, %g GB/s for reading + checksum, %g GB/s for checksum, %d views per node)"
                % (
                    szGB,
                    szGB / (t1 - t0),
                    szGB / dt_distr,
                    szGB / dt_read,
                    szGB / dt_crc,
                    1,
                )
            )
        return l


def load(filename, p={}):

    qlat = qlat_io(filename)

    # check if this is right file format from header
    if not qlat.read_header():
        raise NotImplementedError()

    return qlat.read_lattice_single()

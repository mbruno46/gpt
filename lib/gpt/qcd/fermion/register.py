# warning, this file is automatically generated, do not modify!
def register(reg,op):
  reg.M = lambda dst, src: op.unary(2001,dst,src)
  reg.Mdag = lambda dst, src: op.unary(2002,dst,src)
  reg.Meooe = lambda dst, src: op.unary(2003,dst,src)
  reg.MeooeDag = lambda dst, src: op.unary(2004,dst,src)
  reg.Mooee = lambda dst, src: op.unary(2005,dst,src)
  reg.MooeeDag = lambda dst, src: op.unary(2006,dst,src)
  reg.MooeeInv = lambda dst, src: op.unary(2007,dst,src)
  reg.MooeeInvDag = lambda dst, src: op.unary(2008,dst,src)
  reg.Mdiag = lambda dst, src: op.unary(2009,dst,src)
  reg.Dminus = lambda dst, src: op.unary(2010,dst,src)
  reg.DminusDag = lambda dst, src: op.unary(2011,dst,src)
  reg.ImportPhysicalFermionSource = lambda dst, src: op.unary(2012,dst,src)
  reg.ImportUnphysicalFermion = lambda dst, src: op.unary(2013,dst,src)
  reg.ExportPhysicalFermionSolution = lambda dst, src: op.unary(2014,dst,src)
  reg.ExportPhysicalFermionSource = lambda dst, src: op.unary(2015,dst,src)

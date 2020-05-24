/*************************************************************************************

    Copyright (C) 2019

    Author: Mattia Bruno

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
/*  END LEGAL */

#pragma once

#include <Grid/Grid.h>

NAMESPACE_BEGIN(Grid);

/************************************************************************
 * Wilson Exp Clover
 *
 * exp-term = exp[sum_{mu,nu} csw/(4+m0) * i/4 * sigma_{mu,nu} F_{mu,nu} ]
 *
 * Full Dirac operator
 *
 * (Dee + Doo) psi =  (4+m0) * exp-term psi
 *
 ***********************************************************************/

template <class Impl>
class WilsonExpCloverFermion : public WilsonFermion<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
  template <typename vtype> using iImplClover = iScalar<iMatrix<iMatrix<vtype, Impl::Dimension>, Ns>>;
  typedef iImplClover<Simd> SiteCloverType;
  typedef Lattice<SiteCloverType> CloverFieldType;

public: 
  typedef WilsonFermion<Impl> WilsonBase;

  virtual int    ConstEE(void)     { return 0; };
  virtual void Instantiatable(void) {};
  // Constructors
  WilsonExpCloverFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                         GridRedBlackCartesian &Hgrid,
                         const RealD _mass,
                         const RealD _csw_r,
                         const RealD _csw_t = -1.0,
                         int isw = 0,
                         const WilsonAnisotropyCoefficients &clover_anisotropy = WilsonAnisotropyCoefficients(),
                         const ImplParams &impl_p = ImplParams()) : WilsonFermion<Impl>(_Umu, Fgrid, Hgrid,
                                                                                        _mass, impl_p, clover_anisotropy),
                                                                    CloverTerm(&Fgrid),
                                                                    CloverTermInv(&Fgrid),
                                                                    CloverTermEven(&Hgrid),
                                                                    CloverTermOdd(&Hgrid),
                                                                    CloverTermInvEven(&Hgrid),
                                                                    CloverTermInvOdd(&Hgrid),
                                                                    CloverTermDagEven(&Hgrid),
                                                                    CloverTermDagOdd(&Hgrid),
                                                                    CloverTermInvDagEven(&Hgrid),
                                                                    CloverTermInvDagOdd(&Hgrid)
   {
      assert(Nd == 4); // require 4 dimensions
      assert(!clover_anisotropy.isAnisotropic);
      
      diag_mass = 4.0 + _mass;
      csw_r = _csw_r * 0.5;
      if (_csw_t==-1.0)
         csw_t = csw_r;
      else
         csw_t = _csw_t * 0.5;
      
      if (csw_r!=0.0)
      {
         std::cout << GridLogMessage << "WilsonExpCloverFermion (csw_r="<<csw_r<<", csw_t="<<csw_t<<", mass="<<_mass<<")\n";
         ImportGauge(_Umu,isw);
         std::cout << GridLogMessage << "WilsonExpCloverFermion initialized, CloverTerm computed\n";
      }
      else
      {
         std::cout << GridLogMessage << "WilsonExpCloverFermion (csw_r="<<csw_r<<", csw_t="<<csw_t<<", mass="<<_mass<<")\n";
         ImportGauge_noClover(_Umu);
         std::cout << GridLogMessage << "WilsonExpCloverFermion initialized, CloverTerm set to zero\n";
      }
   }

  void ImportClover(WilsonExpCloverFermion &D)
  {
     CloverTerm=D.CloverTerm;
     CloverTermInv=D.CloverTermInv;

     // Separate the even and odd parts
     pickCheckerboard(Even, CloverTermEven, CloverTerm);
     pickCheckerboard(Odd, CloverTermOdd, CloverTerm);

     pickCheckerboard(Even, CloverTermDagEven, adj(CloverTerm));
     pickCheckerboard(Odd, CloverTermDagOdd, adj(CloverTerm));

     pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
     pickCheckerboard(Odd, CloverTermInvOdd, CloverTermInv);
   
     pickCheckerboard(Even, CloverTermInvDagEven, adj(CloverTermInv));
     pickCheckerboard(Odd, CloverTermInvDagOdd, adj(CloverTermInv));
  }

  // *NOT* EO
  //virtual RealD M(const FermionField &in, FermionField &out);
  RealD M(const FermionField &in, FermionField &out)
  {
     FermionField temp(out.Grid());

     // Wilson term
     out.Checkerboard() = in.Checkerboard();
     this->Dhop(in, out, DaggerNo);

     // Clover+mass term
     Mooee(in, temp);

     out += temp;
     return norm2(out);
  }

  RealD Mdag(const FermionField &in, FermionField &out)
  {
     FermionField temp(out.Grid());

     // Wilson term
     out.Checkerboard() = in.Checkerboard();
     this->Dhop(in, out, DaggerYes);

     // Clover+mass term
     MooeeDag(in, temp);

     out += temp;
     return norm2(out);
  }
  
  void Mooee(const FermionField &in, FermionField &out) {this->MooeeInternal(in, out, DaggerNo, InverseNo);}
  void MooeeDag(const FermionField &in, FermionField &out) {this->MooeeInternal(in, out, DaggerYes, InverseNo);}
  void MooeeInv(const FermionField &in, FermionField &out) {this->MooeeInternal(in, out, DaggerNo, InverseYes);}
  void MooeeInvDag(const FermionField &in, FermionField &out) {this->MooeeInternal(in, out, DaggerYes, InverseYes);}
  
  virtual void MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv)
  {
     out.Checkerboard() = in.Checkerboard();
     CloverFieldType *Clover;
     assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);

     if (dag)
     {
        if (in.Grid()->_isCheckerBoarded)
        {
           if (in.Checkerboard() == Odd)
           {
              Clover = (inv) ? &CloverTermInvDagOdd : &CloverTermDagOdd;
           }
           else
           {
              Clover = (inv) ? &CloverTermInvDagEven : &CloverTermDagEven;
           }
           out = *Clover * in;
        }
        else
        {
           Clover = (inv) ? &CloverTermInv : &CloverTerm;
           out = adj(*Clover) * in;
        }
     }
     else
     {
        if (in.Grid()->_isCheckerBoarded)
        {

           if (in.Checkerboard() == Odd)
           {
              //  std::cout << "Calling clover term Odd" << std::endl;
              Clover = (inv) ? &CloverTermInvOdd : &CloverTermOdd;
           }
           else
           {
              //  std::cout << "Calling clover term Even" << std::endl;
              Clover = (inv) ? &CloverTermInvEven : &CloverTermEven;
           }
           out = *Clover * in;
           //  std::cout << GridLogMessage << "*Clover.Checkerboard() "  << (*Clover).Checkerboard() << std::endl;
        }
        else
        {
           Clover = (inv) ? &CloverTermInv : &CloverTerm;
           out = *Clover * in;
        }
     }

  } // MooeeInternal

  
  void ImportGauge(const GaugeField &_Umu,int isw)
  {
     WilsonFermion<Impl>::ImportGauge(_Umu);
     GridBase *grid = _Umu.Grid();
     typename Impl::GaugeLinkField Bx(grid), By(grid), Bz(grid), Ex(grid), Ey(grid), Ez(grid);

     // Compute the field strength terms mu>nu
     WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Zdir, Ydir);
     WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
     WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Ydir, Xdir);
     WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
     WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
     WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

     // Compute the Clover Operator acting on Colour and Spin
     // multiply here by the clover coefficients for the anisotropy
     RealD f_r, f_t;
     if (isw==0)
     {
        f_r = csw_r;
        f_t = csw_t;
     }
     else
     {
        f_r = csw_r / diag_mass * 0.5;
        f_t = csw_t / diag_mass * 0.5;
     }
     /* computes csw/(4+m0) * sum_{mu,nu} (i/4) * sigma_{mu,nu} * F_{mu,nu} */
     CloverTerm  = fillCloverYZ(Bx) * f_r;
     CloverTerm += fillCloverXZ(By) * f_r;
     CloverTerm += fillCloverXY(Bz) * f_r;
     CloverTerm += fillCloverXT(Ex) * f_t;
     CloverTerm += fillCloverYT(Ey) * f_t;
     CloverTerm += fillCloverZT(Ez) * f_t;

     if (isw==0)
     {
        int lvol = _Umu.Grid()->lSites();
        int DimRep = Impl::Dimension;

        Eigen::MatrixXcd EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
        Eigen::MatrixXcd EigenInvCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);

        Coordinate lcoor;
        typename SiteCloverType::scalar_object Qx = Zero(), Qxinv = Zero();

        thread_for(site, lvol, {
           grid->LocalIndexToLocalCoor(site, lcoor);
           EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
           peekLocalSite(Qx, CloverTerm, lcoor);
           Qxinv = Zero();
           for (int j = 0; j < Ns; j++)
              for (int k = 0; k < Ns; k++)
                 for (int a = 0; a < DimRep; a++)
                    for (int b = 0; b < DimRep; b++){
                       auto zz =  Qx()(j, k)(a, b);
                       EigenCloverOp(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
                    }

           EigenInvCloverOp = EigenCloverOp.inverse();
           for (int j = 0; j < Ns; j++)
              for (int k = 0; k < Ns; k++)
                 for (int a = 0; a < DimRep; a++)
                    for (int b = 0; b < DimRep; b++)
                       Qxinv()(j, k)(a, b) = EigenInvCloverOp(a + j * DimRep, b + k * DimRep);
           pokeLocalSite(Qxinv, CloverTermInv, lcoor);
        });
     }
     else
     {     
        RealD R=3.*csw_r/diag_mass;
        int NMAX=this->get_NMAX(CloverTerm,R);
        std::cout << GridLogMessage << "ExponentialClover NMAX = " << NMAX << std::endl;
   
        CloverFieldType tmp(CloverTerm), Clover(CloverTerm);

        setIdentity(CloverTerm);
        setIdentity(CloverTermInv);

        double f=1.0;
        for (int l=1; l<=NMAX; l++)
        {
           CloverTerm = CloverTerm + tmp*f;
           if ((l%2)==0)
              CloverTermInv = CloverTermInv + tmp*f;
           else
              CloverTermInv = CloverTermInv - tmp*f;

           f *= 1./double(l);
           tmp = tmp*Clover;
        }
        CloverTerm *= diag_mass;
        CloverTermInv *= (1.0/diag_mass);
     }

     // Separate the even and odd parts
     pickCheckerboard(Even, CloverTermEven, CloverTerm);
     pickCheckerboard(Odd, CloverTermOdd, CloverTerm);
   
     pickCheckerboard(Even, CloverTermDagEven, adj(CloverTerm));
     pickCheckerboard(Odd, CloverTermDagOdd, adj(CloverTerm));
   
     pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
     pickCheckerboard(Odd, CloverTermInvOdd, CloverTermInv);
   
     pickCheckerboard(Even, CloverTermInvDagEven, adj(CloverTermInv));
     pickCheckerboard(Odd, CloverTermInvDagOdd, adj(CloverTermInv));
  }

  void ImportGauge_noClover(const GaugeField &_Umu) 
  {
     WilsonFermion<Impl>::ImportGauge(_Umu);
     setIdentity(CloverTerm);
     setIdentity(CloverTermInv);
 
     // Separate the even and odd parts
     pickCheckerboard(Even, CloverTermEven, CloverTerm);
     pickCheckerboard(Odd, CloverTermOdd, CloverTerm);

     pickCheckerboard(Even, CloverTermDagEven, adj(CloverTerm));
     pickCheckerboard(Odd, CloverTermDagOdd, adj(CloverTerm));

     pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
     pickCheckerboard(Odd, CloverTermInvOdd, CloverTermInv);
   
     pickCheckerboard(Even, CloverTermInvDagEven, adj(CloverTermInv));
     pickCheckerboard(Odd, CloverTermInvDagOdd, adj(CloverTermInv));
  }

public:
  CloverFieldType CloverTerm, CloverTermInv;                 // Clover term
private:
  RealD csw_r;                                               // Clover coefficient - spatial
  RealD csw_t;                                               // Clover coefficient - temporal
  RealD diag_mass;                                           // Mass term
  CloverFieldType CloverTermEven, CloverTermOdd;             // Clover term EO
  CloverFieldType CloverTermInvEven, CloverTermInvOdd;       // Clover term Inv EO
  CloverFieldType CloverTermDagEven, CloverTermDagOdd;       // Clover term Dag EO
  CloverFieldType CloverTermInvDagEven, CloverTermInvDagOdd; // Clover term Inv Dag EO

  int get_NMAX(RealD prec, RealD R)
  {
     /* compute stop condition for exponential */
     int NMAX=1;
     RealD cond=R*R/2.;

     while (cond*std::exp(R)>prec)
     {
        NMAX++;
        cond*=R/(double)(NMAX+1);
     }
     return NMAX;
  }
  int get_NMAX(Lattice<iImplClover<vComplexD>> &t, RealD R) {return get_NMAX(1e-12,R);}
  int get_NMAX(Lattice<iImplClover<ComplexD>> &t, RealD R) {return get_NMAX(1e-12,R);}
  int get_NMAX(Lattice<iImplClover<vComplexF>> &t, RealD R) {return get_NMAX(1e-6,R);}
  int get_NMAX(Lattice<iImplClover<ComplexF>> &t, RealD R) {return get_NMAX(1e-6,R);}

  template <typename vtype> void setIdentity(Lattice<iImplClover<vtype>> &Clover)
  {
     int os,is,mu,nu,ca,cb;
     Coordinate lcoor;
     typedef typename vtype::scalar_type stype;
     iImplClover<stype> mat;

     for (mu=0;mu<Nd;mu++)
        for (nu=0;nu<Nd;nu++)
           for (ca=0;ca<Nc;ca++)
              for (cb=0;cb<Nc;cb++)
                 mat()(mu,nu)(ca,cb) = (mu==nu && ca==cb) ? stype(1.0,0.0) : stype(0.,0.);

     for (os=0;os<Clover.Grid()->lSites();os++)
     {
        Clover.Grid()->LocalIndexToLocalCoor(os, lcoor);
        pokeLocalSite(mat,Clover,lcoor);
     }
  }

  // eventually these can be compressed into 6x6 blocks instead of the 12x12
  // using the DeGrand-Rossi basis for the gamma matrices
  CloverFieldType fillCloverYZ(const GaugeLinkField &F)
  {
    CloverFieldType T(F.Grid());
    T = Zero();
    auto T_v = T.View();
    auto F_v = F.View();
    thread_for(i, CloverTerm.Grid()->oSites(),
    {
      T_v[i]()(0, 1) = timesMinusI(F_v[i]()());
      T_v[i]()(1, 0) = timesMinusI(F_v[i]()());
      T_v[i]()(2, 3) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 2) = timesMinusI(F_v[i]()());
    });

    return T;
  }

  CloverFieldType fillCloverXZ(const GaugeLinkField &F)
  {
    CloverFieldType T(F.Grid());
    T = Zero();
    
    auto T_v = T.View();
    auto F_v = F.View();
    thread_for(i, CloverTerm.Grid()->oSites(),
    {
      T_v[i]()(0, 1) = -F_v[i]()();
      T_v[i]()(1, 0) = F_v[i]()();
      T_v[i]()(2, 3) = -F_v[i]()();
      T_v[i]()(3, 2) = F_v[i]()();
    });

    return T;
  }

  CloverFieldType fillCloverXY(const GaugeLinkField &F)
  {
    CloverFieldType T(F.Grid());
    T = Zero();

    auto T_v = T.View();
    auto F_v = F.View();
    thread_for(i, CloverTerm.Grid()->oSites(),
    {
      T_v[i]()(0, 0) = timesMinusI(F_v[i]()());
      T_v[i]()(1, 1) = timesI(F_v[i]()());
      T_v[i]()(2, 2) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 3) = timesI(F_v[i]()());
    });

    return T;
  }

  CloverFieldType fillCloverXT(const GaugeLinkField &F)
  {
    CloverFieldType T(F.Grid());
    T = Zero();

    auto T_v = T.View();
    auto F_v = F.View();
    thread_for(i, CloverTerm.Grid()->oSites(),
    {
      T_v[i]()(0, 1) = timesI(F_v[i]()());
      T_v[i]()(1, 0) = timesI(F_v[i]()());
      T_v[i]()(2, 3) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 2) = timesMinusI(F_v[i]()());
    });

    return T;
  }

  CloverFieldType fillCloverYT(const GaugeLinkField &F)
  {
    CloverFieldType T(F.Grid());
    T = Zero();
    
    auto T_v = T.View();
    auto F_v = F.View();
    thread_for(i, CloverTerm.Grid()->oSites(),
    {
      T_v[i]()(0, 1) = -(F_v[i]()());
      T_v[i]()(1, 0) = (F_v[i]()());
      T_v[i]()(2, 3) = (F_v[i]()());
      T_v[i]()(3, 2) = -(F_v[i]()());
    });

    return T;
  }

  CloverFieldType fillCloverZT(const GaugeLinkField &F)
  {
    CloverFieldType T(F.Grid());

    T = Zero();

    auto T_v = T.View();
    auto F_v = F.View();
    thread_for(i, CloverTerm.Grid()->oSites(),
    {
      T_v[i]()(0, 0) = timesI(F_v[i]()());
      T_v[i]()(1, 1) = timesMinusI(F_v[i]()());
      T_v[i]()(2, 2) = timesMinusI(F_v[i]()());
      T_v[i]()(3, 3) = timesI(F_v[i]()());
    });

    return T;
  }

};

typedef WilsonExpCloverFermion<WilsonImplF> WilsonExpCloverFermionF;
typedef WilsonExpCloverFermion<WilsonImplD> WilsonExpCloverFermionD;

NAMESPACE_END(Grid);

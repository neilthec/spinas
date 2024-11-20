
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2023 Neil Christensen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//File:  SPINAS/SM/uuuu.h

namespace spinas {

  class uuuu : public process {
  private:
    ldouble e, gs;//Electric Charge
    ldouble mu, mh, wh, MZ, WZ, MW, SW, CW;//Mass of u, c, Higgs, width of Higgs, Z mass, Z width, W-boson mass and sin and cos of the Weinberg angle.
    particle p1,p2,p3,p4;
    propagator propAG, proph, propZ;
    cdouble pDenSAG, pDenSh, pDenSZ;
    cdouble pDenTAG, pDenTh, pDenTZ;
    sproduct a13a, s13s, a14a, s14s, a23a, s23s, a24a, s24s;
    sproduct s12s, a12a, s34s, a34a;
    ldouble preh, gL, gR, preZ, preZ0;

    
  public:
    uuuu(const ldouble& echarge, const ldouble& gscharge,
	 const ldouble& massu,
	 const ldouble& massh, const ldouble& widthh,
	 const ldouble& massW, const ldouble& sinW, const ldouble& widthZ);

    //Set Masses
    void set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp_gluon_S(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    cdouble amp_rest_S(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    cdouble amp_gluon_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    cdouble amp_rest_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    

    
    
    
  };
  //Tests
  int test_uuuu();
}

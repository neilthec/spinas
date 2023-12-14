
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

//File:  SPINAS/SM/gdZd.h

namespace spinas {

  class gdZd : public process {
  private:
    ldouble e, Qd, gs;//Electric strong Charge
    ldouble md, MW, SW, CW, MZ;//Mass of d and h and W, width of h and Z and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator propd;
    cdouble pDenS, pDenU;
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<14>,[14]
    sproduct s12s, a12a, s23s, a23a, s13s, a13a, s34s, a34a, s14s, a14a;
    //[1421], <1421>
    sproduct s1421s, a1421a;
    ldouble sqrt2;
    ldouble preTU, preh, gL, gR;

    
  public:
    //gdZd();
    gdZd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massd, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    ldouble amp2_gplus();//Positive Helicity Photon only.
    ldouble amp2_gminus();//Minus Helicity Photon only.
    

    
    
  };
  //Tests
  int test_gdZd();
}

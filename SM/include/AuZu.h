
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

//File:  SPINAS/SM/AuZu.h

namespace spinas {

  class AuZu : public process {
  private:
    ldouble e, Qu;//Electric Charge
    ldouble mu, MW, SW, CW, MZ;//Mass of e and h and W, width of h and Z and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator prope;
    cdouble pDenS, pDenU;
    //Spinor Products
    sproduct s12s, s13s, s14s, a23a, a34a, s123a, s132a, s134a, s143a;//Spinor Products
    sproduct a12a, a13a, a14a, s23s, s34s, a123s, a132s, a134s, a143s;
    ldouble sqrt2;
    ldouble preTU, preh, gLu, gRu;

    
  public:
    //AuZu();
    AuZu(const ldouble& echarge, const ldouble& massu, const ldouble& massW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massu, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    ldouble amp2_Aplus();//Positive Helicity Photon only.
    ldouble amp2_Aminus();//Minus Helicity Photon only.
    ldouble amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
    

    
    
  };
  //Tests
  int test_AuZu();
}

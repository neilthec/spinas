
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

//File:  SPINAS/SM/AeZe.h

namespace spinas {

  class AeZe : public process {
  private:
    ldouble e;//Electric Charge
    ldouble me, MW, SW, CW, MZ;//Mass of e and h and W, width of h and Z and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator prope;
    cdouble pDenU, pDenS;
    //Spinor Products
    sproduct s12s, a12a, s13s, a13a, s14s, a14a, a23a, s23s, a34a, s34s, s123a, a123s, s132a, a132s, s134a, a134s, s143a, a143s;
    ldouble sqrt2;
    ldouble preSU, preh, gL, gR;

    
  public:
    //AeZe();
    AeZe(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& masse, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    ldouble amp2_Aplus();//Positive Helicity Photon only.
    ldouble amp2_Aminus();//Negative Helicity Photon only.

    

    
    
  };
  //Tests
  int test_AeZe();
}

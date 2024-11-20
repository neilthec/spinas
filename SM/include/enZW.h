
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

//File:  SPINAS/SM/enZW.h

namespace spinas {

  class enZW : public process {
  private:
    ldouble e;//Electric Charge
    ldouble me, MW, MZ, SW, CW, WW;//Mass of W, width of W and sin(theta_W)
    ldouble preW,pree,pren,gL,gR;
    particle p1,p2,p3,p4;
    propagator propW, prope, propne;
    cdouble pDenS, pDenT, pDenU;
    //<34>,[34],<23>,[14],<24>,[13],<12>,[132>
    sproduct a34a, s34s, a23a, s14s, a24a, s13s, a12a, s132a;
    //<13>,[413>
    sproduct a13a, s413a;
    //[314>
    sproduct s314a;
    

    
  public:
    //enZW();
    enZW(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthW);

    //Set Masses
    void set_masses(const ldouble& masse, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    

    
    
  };
  //Tests
  int test_enZW();

}

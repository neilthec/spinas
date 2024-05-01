
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

//File:  SPINAS/SM/WWWW2.h

namespace spinas {

  class WWWW2 : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mh, wh, MW, WZ, SW, CW, MZ;//Mass of h, W and Z, width of h and W and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator proph, propZ, propA;
    cdouble pDenTh, pDenSh, pDenTZ, pDenSZ, pDenTA, pDenSA;
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    sproduct s14s, a14a, a32a, s32s, s42s, a42a, s43s, a43a, a12a, s12s, a13a, s13s;
    sproduct s231a, s132a, s123a, s321a;
    ldouble preh, preA, preZ;

    
  public:
    WWWW2(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthZ, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    

    
    
  };
  //Tests
  int test_WWWW2();

}

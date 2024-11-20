
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

//File:  SPINAS/SM/ZhWW.h

namespace spinas {

  class ZhWW : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mh, MW, MZ, SW, CW, WW, WZ;//Mass of h, Z and W, width of Z and W and sin(theta_W) 
    particle p1,p2,p3,p4;
    propagator propW, propZ;
    cdouble pDenS, pDenT, pDenU;
    //<34>,[34],<14>,[14],<13>,[13]
    sproduct a34a, s34s, a14a, s14s, a13a, s13s;
    //[131>,[141>,[431>,[424>,[341>,[323>,[123>,[124>,[321>,[421>
    sproduct s131a, s141a, s431a, s424a, s341a, s323a, s123a, s124a, s321a, s421a;
    ldouble sqrt2;
    ldouble pre, preS, preTU;

    
  public:
    //ZhWW();
    ZhWW(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW, const ldouble& widthZ);

    //Set Masses
    void set_masses(const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    

    
    
  };
  //Tests
  int test_ZhWW();

}


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

//File:  SPINAS/SM/udAW.h

namespace spinas {

  class udAW : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mu, md, MW, SW, WW;//Mass of W, width of W and sin(theta_W)
    ldouble pre;
    particle p1,p2,p3,p4;
    propagator propW, propu, propd;
    cdouble pDenS, pDenT, pDenU, prop, propAp, propAm;
    //<42>,[41],[34],[31],[3123]
    sproduct a42a,s41s,s34s,s31s,s3123s;
    sproduct a3123a,a34a,a32a;
    

    
  public:
    //udAW();
    udAW(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthW);

    //Set Masses
    void set_masses(const ldouble& masse, const ldouble& massd, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    

    
    
  };
  //Tests
  int test_udAW();

}

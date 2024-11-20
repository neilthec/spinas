
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

//File:  SPINAS/SM/enAW.h

namespace spinas {

  class enAW : public process {
  private:
    ldouble e;//Electric Charge
    ldouble me, MW, SW, WW;//Mass of W, width of W and sin(theta_W)
    ldouble pre;
    particle p1,p2,p3,p4;
    propagator propW, prope;
    cdouble pDenS, pDenT;
    //[34],<12>,<24>,[13],[314>
    sproduct s34s, a12a, a24a, s13s, s314a;
    //<23>,[14],<34>,[413>
    sproduct a23a,s14s,a34a,s413a;
    

    
  public:
    //enAW();
    enAW(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthW);

    //Set Masses
    void set_masses(const ldouble& masse, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    

    
    
  };
  //Tests
  int test_enAW();

}

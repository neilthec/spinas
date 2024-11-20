
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

//File:  SPINAS/SM/AWZW.h

namespace spinas {

  class AWZW : public process {
  private:
    ldouble e;//Electric Charge
    ldouble MW, MZ, SW, CW, WW;//Mass of W and Z, width of W and sin(theta_W) 
    particle p1,p2,p3,p4;
    propagator propW;
    cdouble pDenS, pDenU;
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    sproduct s12s, a12a, s24s, a24a, s23s, a23a, a34a, s34s, a14a, s14s, a13a, s13s;
    ldouble sqrt2;
    ldouble pre;

    
  public:
    //AWZW();
    AWZW(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthW);

    //Set Masses
    void set_masses(const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    ldouble amp2_Aplus();
    ldouble amp2_Aplus_params(const ldouble params[25]);
    

    
    
  };
  //Tests
  int test_AWZW();

}

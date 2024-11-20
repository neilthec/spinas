
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

//File:  SPINAS/SM/include/gggg.h

namespace spinas {
  
  class gggg : public process {
  private:
    ldouble gs;//Electric Charge
    particle p1,p2,p3,p4;
    sproduct s12s, a12a, s13s, a13a, s41s, a41a, s24s, a24a, s23s, a23a, s34s, a34a;
    sproduct s1432s, a1432a, s3214s, a3214a;
    sproduct a3124a;
    propagator prop;
    ldouble pDenS, pDenT, pDenU;

    
    
  public:
    //gggg();
    gggg(const ldouble& gcharge);
    
    //Set Masses
    void set_masses(const ldouble& masse);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
    
    
    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2_Feynman();
    ldouble amp2();
    //Helicity Cases for Tests
    ldouble amp2_gplus_gplus();
    ldouble amp2_gplus_gminus();
    ldouble amp2_gminus_gplus();
    ldouble amp2_gminus_gminus();
    
  };
  //Tests
  int test_gggg();
}

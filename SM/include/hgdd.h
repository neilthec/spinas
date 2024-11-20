
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

//File:  SPINAS/SM/hgdd.h

namespace spinas {

  class hgdd : public process {
  private:
    ldouble e, gs;//Electric Charge, strong coupling constant
    ldouble md,mh,MW,SW;//Mass of d and h and W and sin thetaW
    ldouble pre;//prefactor
    particle p1,p2,p3,p4;
    propagator prop;
    cdouble pDenT, pDenU;
    sproduct s24s, a24a, s23s, a23a, s34s, a34a, s213a, a213s, s214a, a214s, s2312s, a2312a;
    ldouble sqrt2;

    
  public:
    //hgdd();
    hgdd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    ldouble amp2_gplus();//Positive Helicity Gluon only.

    

    
    
  };
  //Tests
  int test_hgdd();
}

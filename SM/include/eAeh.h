
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

//File:  SPINAS/SM/eAeh.h

namespace spinas {

  class eAeh : public process {
  private:
    ldouble e;//Electric Charge
    ldouble me,mh;//Mass of e and h
    ldouble v;//vev
    particle p1,p2,p3,p4;
    propagator prop;
    cdouble pDenS, pDenU;
    sproduct s12s, a12a, s23s, a23a, s13s, a13a, s243a, a243s, s241a, a241s, s2342s, a2342a;
    ldouble sqrt2;

    
  public:
    //eAeh();
    eAeh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& vev);

    //Set Masses
    void set_masses(const ldouble& masse, const ldouble& massh);
    //Set v
    void set_v(const ldouble& vev);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3);//Double the spins
    ldouble amp2();
    ldouble amp2_Aplus();//Positive Helicity Photon only.
    ldouble amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
    

    
    
  };
  //Tests
  int test_eAeh();
}

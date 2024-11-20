
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

//File:  SPINAS/SM/hhhh.h

namespace spinas {

  class hhhh : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mh, wh, MW, SW;//Mass of e, Higgs, width of Higgs, mass of W, sin of Weinberg angle.
    particle p1,p2,p3,p4;
    propagator proph;
    cdouble pDenS, pDenT, pDenU;
    ldouble preh, prehSTU;

    
  public:
    hhhh(const ldouble& echarge,
	 const ldouble& massh, const ldouble& widthh,
	 const ldouble& massW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp();
    ldouble amp2();
    

    
    
    
  };
  //Tests
  int test_hhhh();
}

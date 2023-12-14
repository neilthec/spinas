
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

//File:  SPINAS/SM/dnnd.h

namespace spinas {

  class dnnd : public process {
  private:
    ldouble e;//Electric Charge
    ldouble md, MZ, WZ, MW, SW, CW;//Mass of u, s, Higgs, width of Higgs, Z mass, Z width, W-boson mass and sin and cos of the Weinberg angle.
    particle p1,p2,p3,p4;
    propagator propZ;
    cdouble pDenUZ;
    sproduct s34s, a12a, s13s, a24a;
    ldouble preh, gLd, gRd, gLn, gRn, preZ, preZ0;

    
  public:
    dnnd(const ldouble& echarge,
	 const ldouble& massd,
	 const ldouble& massW, const ldouble& sinW, const ldouble& widthZ);

    //Set Masses
    void set_masses(const ldouble& massu, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2);//Double the spins
    ldouble amp2();
    

    
    
    
  };
  //Tests
  int test_dnnd();
}

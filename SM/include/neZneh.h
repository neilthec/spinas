
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

//File:  SPINAS/SM/neZneh.h

namespace spinas {

  class neZneh : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mh, MW, SW, CW, MZ, WZ;//Mass of h and W, width of h and Z and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator propZ;
    cdouble pDenT;
    //[23], <13>
    sproduct s23s, a12a;
    ldouble sqrt2;
    ldouble preZ;

    
  public:
    //neZneh();
    neZneh(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ);

    //Set Masses
    void set_masses(const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds3);//Double the spins
    ldouble amp2();
    ldouble amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
    

    
    
  };
  //Tests
  int test_neZneh();
}

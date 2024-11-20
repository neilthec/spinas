
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

//File:  SPINAS/SM/enWh.h

namespace spinas {

  class enWh : public process {
  private:
    ldouble e;//Electric Charge
    ldouble me, mh, MW, WW, SW;//Mass of e and h and W, width ofW and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator prope, propW;
    cdouble pDenWS, pDeneU;
    sproduct a23a, s13s, a12a, s343a, s341a;
    ldouble pree, preW;

    
  public:
    //enWh();
    enWh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds3);//Double the spins
    ldouble amp2();
    
    
  };
  //Tests
  int test_enWh();
}

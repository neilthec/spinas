
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2024 Neil Christensen

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

//File:  SPINAS/SM/udWh.h

namespace spinas {

  class udWh : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mu, md, mh, MW, WW, SW;//Mass of e and h and W, width ofW and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator propu, propd, propW;
    cdouble pDenS, pDenU, pDenT;
    //<13>,[23],<12>,[12],[343>
    sproduct a13a, s23s, a12a, s12s, s343a;
    //<13>,[23],[342>
    sproduct s342a;
    //[23],<13>,[143>
    sproduct s143a;
    ldouble preu, pred, preW;

    
  public:
    //udWh();
    udWh(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3);//Double the spins
    ldouble amp2();
    
    
  };
  //Tests
  int test_udWh();
}

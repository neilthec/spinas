
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

//File:  SPINAS/SM/ZWWh.h

namespace spinas {

  class ZWWh : public process {
  private:
    ldouble e;//Electric Charge
    ldouble mh, MW, MZ, SW, CW, WW, WZ;//Mass of h, Z and W, width of Z and W and sin(theta_W) 
    particle p1,p2,p3,p4;
    propagator propW, propZ;
    cdouble pDenS, pDenT, pDenU;
    //<23>,[23],<12>,[12],<13>,[13]
    sproduct a23a, s23s, a12a, s12s, a13a, s13s;
    //[131>,[121>,[231>,[242>,[321>,[343>,[143>,[142>,[341>,[241>
    sproduct s131a, s121a, s231a, s242a, s321a, s343a, s143a, s142a, s341a, s241a;
    ldouble sqrt2;
    ldouble pre, preS, preTU;

    
  public:
    //ZWWh();
    ZWWh(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW, const ldouble& widthZ);

    //Set Masses
    void set_masses(const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3);//Double the spins
    ldouble amp2();
    

    
    
  };
  //Tests
  int test_ZWWh();

}

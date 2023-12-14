
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

//File:  SPINAS/SM/uuAh.h

namespace spinas {

  class uuAh : public process {
  private:
    ldouble e, Qu;//Electric Charge
    ldouble mu,mh, MW, SW;//Mass of u and h, W and sin theta_W
    ldouble pre;//prefactor
    particle p1,p2,p3,p4;
    propagator prop;
    cdouble pDenT, pDenU;
    sproduct s13s, a13a, s23s, a23a, s12s, a12a, s342a, a342s, s341a, a341s, s3243s, a3243a;
    ldouble sqrt2;

    
  public:
    //uuAh();
    uuAh(const ldouble& echarge, const ldouble& massu, const ldouble& massh, const ldouble& massW, const ldouble& sinW);

    //Set Masses
    void set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3);//Double the spins
    ldouble amp2();
    ldouble amp2_Aplus();//Positive Helicity Photon only.
    ldouble amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
    

    
    
  };
  //Tests
  int test_uuAh();
}

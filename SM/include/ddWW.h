
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

//File:  SPINAS/SM/ddWW.h

namespace spinas {

  class ddWW : public process {
  private:
    ldouble e, Qf;//Electric Charge
    ldouble mu, md, mh, wh, MW, SW, CW, MZ, WZ;//Mass of e and h and W, width of h and Z and sin(theta_W)
    particle p1,p2,p3,p4;
    propagator proph, propu, propZ, propA;
    cdouble pDenhS, pDendT, pDendU, pDenZS, pDenAS;
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    sproduct s12s, a12a, s23s, a23a, s13s, a13a, s34s, a34a, s24s, a24a, s14s, a14a;
    //[413>,[314>,[231>,[132>
    sproduct s413a, s314a, s231a, s132a;
    ldouble sqrt2;
    ldouble pred, preh, preZ, gL, gR;

    
  public:
    //ddWW();
    ddWW(const ldouble& echarge, const ldouble& massd, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ);

    //Set Masses
    void set_masses(const ldouble& massd, const ldouble& massu, const ldouble& massh, const ldouble& massW);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
    ldouble amp2();
    ldouble amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
    
    
  };
  //Tests
  int test_ddWW();
}

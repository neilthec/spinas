
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

//File:  SPINAS/user-dir/eemm.h

class eemm : public spinas::process {
 private:
  ldouble e;//Electric Charge
  ldouble me,mm;//Mass of e and mu
  spinas::particle p1,p2,p3,p4;
  spinas::propagator prop;
  cdouble pDenS;
  spinas::sproduct a13a, s13s, a14a, s14s, a23a, s23s, a24a, s24s;
  
  
 public:
  eemm(const ldouble& echarge, const ldouble& masse, const ldouble& massmu);
  
  //Set Masses
  void set_masses(const ldouble& masse, const ldouble& massmu);
  //Set Momenta
  void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);
  
  //Amplitude
  cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4);//Double the spins
  ldouble amp2();
  
  
};

//Tests
int test_eemm();



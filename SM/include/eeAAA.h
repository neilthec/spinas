
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2023-2026 Neil Christensen, Nick Majestic

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


#include <random>


namespace spinas {

  class eeAAA : public process {
  private:
    ldouble e;//Electric Charge
    ldouble me;//Mass of e
    particle p1,p2,p3,p4,p5;
    propagator prop;
    cdouble pDenS13, pDenS14, pDenS15, pDenS23, pDenS24, pDenS25;
    sproduct s34s, a34a, s12s, a12a, s13s, a13a, s15s, a15a, s24s, a24a, s23s, a23a, s25s, a25a, s14s, a14a, s35s, a35a, s45s, a45a;
    sproduct s134a, s135a, s143a, s145a, s153a, s154a, s234a, s243a, s235a, s245a, s253a, s254a, s314a, s315a, s324a, s325a, s413a, s414a, s415a, s423a, s425a, s453a, s434a, s513a, s514a, s523a, s524a, s525a;
    sproduct s3145s, s3154s, s3123s, a3123a, s4124s, a4124a, s5125s, a5125a, s4135s, s4153s, s5134s, s5143s;
    ldouble sqrt2;

    
  public:
    //eeAAA();
    eeAAA(const ldouble& echarge, const ldouble& masse);

    //Set Masses
    void set_masses(const ldouble& masse);
    //Set Momenta
    void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4], const ldouble mom5[4]);

    bool check_phase_space(
      const ldouble p1[4],
      const ldouble p2[4],
      const ldouble p3[4],
      const ldouble p4[4],
      const ldouble p5[4],
      ldouble m1,
      ldouble m2,
      ldouble m3,
      ldouble m4,
      ldouble m5,
      ldouble tol
    );

    void make_phase_space(
      ldouble energy,
      ldouble me,
      ldouble E3,
      ldouble theta3,
      ldouble phi3,
      ldouble theta4,
      ldouble phi4,
      ldouble p1[4],
      ldouble p2[4],
      ldouble p3[4],
      ldouble p4[4],
      ldouble p5[4]
    );

    bool make_random_phase_space(
      ldouble energy,
      ldouble me,
      ldouble p1[4],
      ldouble p2[4],
      ldouble p3[4],
      ldouble p4[4],
      ldouble p5[4],
      std::mt19937& rng
    );

    //Amplitude
    cdouble amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5);//Double the spins
    cdouble amp_permutation(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5);
    cdouble amp_feynman(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5);
    cdouble amp_feynman_r(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5);
    ldouble amp2();
    ldouble amp2_feynman();
    ldouble amp2_permutation();
    ldouble amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4], const ldouble mom5[5]);
    

    
    
  };
  //Tests
  int test_eeAAA();
}

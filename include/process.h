
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

//File:  SPINAS/include/process.h

#pragma once

#include <functional>

#include "types.h"

namespace spinas {

  class process{
  private:
    

  public:
    process();

    //Prototype set_momenta
    //Used in test_2to2_amp2 tests.
    virtual void set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]);

    //Functions to determine the number of loops, the normalization factor and the spin indices for massive spin 1.
    //Number of loops
    int get_num_spin_loops(const int& dsA);
    int get_num_spin_loops(const int& dsA, const int& dsB);
    int get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC);
    int get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC, const int& dsD);
    int get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE);
    int get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE, const int& dsF);
    int get_num_spin_loops(const int* dsList, const int length);
    
  //Normalization Factor
    ldouble get_spin_normalization(const int& dsA);
    ldouble get_spin_normalization(const int& dsA, const int& dsB);
    ldouble get_spin_normalization(const int& dsA, const int& dsB, const int& dsC);
    ldouble get_spin_normalization(const int& dsA, const int& dsB, const int& dsC, const int& dsD);
    ldouble get_spin_normalization(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE);
    ldouble get_spin_normalization(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE, const int& dsF);
    ldouble get_spin_normalization(const int* dsList, const int length);

    //Spinor Spins
    void get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& iterator);
    void get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& iterator);
    void get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& iterator);
    void get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& dsD, int& dsDa, int& dsDb, const int& iterator);
    void get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& dsD, int& dsDa, int& dsDb, const int& dsE, int& dsEa, int& dsEb, const int& iterator);
    void get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& dsD, int& dsDa, int& dsDb, const int& dsE, int& dsEa, int& dsEb, const int& dsF, int& dsFa, int& dsFb, const int& iterator);
    void get_spinor_spins(const int* dsList[], int* dsLista[], int* dsListb[], int length, const int& iterator);

    
    //Parameter Scan
    typedef std::function<ldouble(const ldouble*)> amp2ScanFunc;
    ldouble score_2to2_amp2(amp2ScanFunc amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20], const ldouble params[10]);

    //Test against expected amp2
    typedef std::function<ldouble()> amp2Func;
    int test_2to2_amp2(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]);
    int test_2to2_amp2_rotations(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]);
    int test_2to2_amp2_boosts(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]);
    int test_2to2_amp2_boosts_and_rotations(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]);
    void set_test_initializations(const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble &Pin,
					   ldouble &En1, ldouble &En2, ldouble &Pout, ldouble &En3, ldouble &En4) const;
    void set_test_momenta(ldouble p1[4], ldouble p2[4], ldouble p3[4], ldouble p4[4], ldouble En1, ldouble En2,
			  const ldouble &Pin, const ldouble &En3, const ldouble &En4, const ldouble &Pout, const ldouble &cost) const;
    void rotate_random_test_momenta(ldouble p1[4], ldouble p2[4], ldouble p3[4], ldouble p4[4]) const;
    void boost_random_test_momenta(ldouble p1[4], ldouble p2[4], ldouble p3[4], ldouble p4[4]) const;
    void print_test_message(const char *frame_string, const ldouble &m1, const ldouble &m2, const ldouble &m3, const ldouble &m4,
			    const ldouble &ampSquared, const ldouble &amp2_data, const ldouble &cost) const;

  };

  
}

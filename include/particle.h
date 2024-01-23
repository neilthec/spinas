
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

//File:  SPINAS/include/particle.h

#pragma once

#include "types.h"

namespace spinas {

  class particle{
  private:
    ldouble m;//Mass
    ldouble p[4];//Momentum
    ldouble pmag, pxymag;//Magnitude of momentum
    ldouble theta, phi;//Spherical angles
    cdouble c, s, sc;//cos(theta/2), sin(theta/2)*exp(i*phi), s*
    cdouble sqrtEpP, sqrtEmP; //sqrt(E+p), sqrt(E-p)

    //The matrices and spinors
    bool upMatCalculated = false, loMatCalculated = false;
    cmatrix upMat, loMat;
    //Helicity Spinors
    bool m0rangleCalculated = false, m0langleCalculated = false;
    bool m0rsquareCalculated = false, m0lsquareCalculated = false;
    cvector m0rangle, m0langle, m0rsquare, m0lsquare;
    //Spin Spinors
    //rangle
    bool rangleUpperP1Calculated = false, rangleUpperM1Calculated = false;
    cvector rangleUpperP1, rangleUpperM1;
    bool rangleLowerP1Calculated = false, rangleLowerM1Calculated = false;
    cvector rangleLowerP1, rangleLowerM1;
    //langle
    bool langleUpperP1Calculated = false, langleUpperM1Calculated = false;
    cvector langleUpperP1, langleUpperM1;
    bool langleLowerP1Calculated = false, langleLowerM1Calculated = false;
    cvector langleLowerP1, langleLowerM1;
    //lsquare
    bool lsquareUpperP1Calculated = false, lsquareUpperM1Calculated = false;
    cvector lsquareUpperP1, lsquareUpperM1;
    bool lsquareLowerP1Calculated = false, lsquareLowerM1Calculated = false;
    cvector lsquareLowerP1, lsquareLowerM1;
    //rsquare
    bool rsquareUpperP1Calculated = false, rsquareUpperM1Calculated = false;
    cvector rsquareUpperP1, rsquareUpperM1;
    bool rsquareLowerP1Calculated = false, rsquareLowerM1Calculated = false;
    cvector rsquareLowerP1, rsquareLowerM1;


  public:
    particle();
    particle(const ldouble& mass);
    particle(const ldouble momentum[4], const ldouble& mass);
    void set_mass(const ldouble& mass);
    void set_momentum(const ldouble momentum[4]);
    const ldouble mass() const;
    const ldouble get_mass() const;
    const ldouble get_momentum(const int& mu) const;
    void update();
    bool test_angles() const;


    //dot: p1.p2
    ldouble dot(const particle& p2) const;

    //Matrices
    cmatrix umat();//Upper
    cmatrix lmat();//Lower

    //Spinors
    //Massless
    cvector rangle();
    cvector lsquare();
    cvector langle();
    cvector rsquare();

    //Massive
    cvector rangle(const int& spin2);
    cvector rangle(const int& spin2, const bool& upper);
    cvector langle(const int& spin2);
    cvector langle(const int& spin2, const bool& upper);
    cvector lsquare(const int& spin2);
    cvector lsquare(const int& spin2, const bool& upper);
    cvector rsquare(const int& spin2);
    cvector rsquare(const int& spin2, const bool& upper);

    //Massive cmatrix form with both spins
    cmatrix rangle_matrix();
    cmatrix rangle_matrix(const bool& upper);
    cmatrix langle_matrix();
    cmatrix langle_matrix(const bool& upper);
    
    cmatrix lsquare_matrix();
    cmatrix lsquare_matrix(const bool& upper);
    cmatrix rsquare_matrix();
    cmatrix rsquare_matrix(const bool& upper);


    //Error message:
    void usage(const char* msg) const;

    //Generators of Lorentz Rotations with no dot (acting on left chiral angle brackets)
    cmatrix lorentz_j3_lu() const;
    cmatrix lorentz_jp_lu() const;
    cmatrix lorentz_jm_lu() const;

    cmatrix lorentz_j3_ul() const;
    cmatrix lorentz_jp_ul() const;
    cmatrix lorentz_jm_ul() const;
    
    //Generators of Lorentz Rotations with a dot (acting on right chiral square brackets)
    cmatrix lorentz_j3_lu_dot() const;
    cmatrix lorentz_jp_lu_dot() const;
    cmatrix lorentz_jm_lu_dot() const;

    cmatrix lorentz_j3_ul_dot() const;
    cmatrix lorentz_jp_ul_dot() const;
    cmatrix lorentz_jm_ul_dot() const;
    
    //Generators of spin
    cmatrix spin_j3_lu() const;
    cmatrix spin_jp_lu() const;
    cmatrix spin_jm_lu() const;
    
    cmatrix spin_j3_ul() const;
    cmatrix spin_jp_ul() const;
    cmatrix spin_jm_ul() const;
    

  };

  int test_particle();
}



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
    ldouble sqrt2;
    cdouble eppz, empz;//E+p, E-p
    cdouble pxppy, pxmpy;//px+ipy, px-ipy

    // 2-dimensional 
    //The matrices and spinors
    bool upMat2dimCalculated = false, loMat2dimCalculated = false;
    cmatrix upMat2dim, loMat2dim;
    //Helicity Spinors
    bool m0rangle2dimCalculated = false, m0langle2dimCalculated = false;
    bool m0rsquare2dimCalculated = false, m0lsquare2dimCalculated = false;
    cvector m0rangle2dim, m0langle2dim, m0rsquare2dim, m0lsquare2dim;
    //Spin Spinors
    //rangle
    bool rangleUpperP12dimCalculated = false, rangleUpperM12dimCalculated = false;
    cvector rangleUpperP12dim, rangleUpperM12dim;
    bool rangleLowerP12dimCalculated = false, rangleLowerM12dimCalculated = false;
    cvector rangleLowerP12dim, rangleLowerM12dim;
    //langle
    bool langleUpperP12dimCalculated = false, langleUpperM12dimCalculated = false;
    cvector langleUpperP12dim, langleUpperM12dim;
    bool langleLowerP12dimCalculated = false, langleLowerM12dimCalculated = false;
    cvector langleLowerP12dim, langleLowerM12dim;
    //lsquare
    bool lsquareUpperP12dimCalculated = false, lsquareUpperM12dimCalculated = false;
    cvector lsquareUpperP12dim, lsquareUpperM12dim;
    bool lsquareLowerP12dimCalculated = false, lsquareLowerM12dimCalculated = false;
    cvector lsquareLowerP12dim, lsquareLowerM12dim;
    //rsquare
    bool rsquareUpperP12dimCalculated = false, rsquareUpperM12dimCalculated = false;
    cvector rsquareUpperP12dim, rsquareUpperM12dim;
    bool rsquareLowerP12dimCalculated = false, rsquareLowerM12dimCalculated = false;
    cvector rsquareLowerP12dim, rsquareLowerM12dim;

    // 3-dimensional 
    //The matrices and spinors
    bool upMat3dimCalculated = false, loMat3dimCalculated = false;
    cmatrix upMat3dim, loMat3dim;
    //Helicity Spinors
    bool m0rangle3dimCalculated = false, m0langle3dimCalculated = false;
    bool m0rsquare3dimCalculated = false, m0lsquare3dimCalculated = false;
    cvector m0rangle3dim, m0langle3dim, m0rsquare3dim, m0lsquare3dim;
    //Spin Spinors
    //rangle
    bool rangleUpperP13dimCalculated = false, rangleUpperM13dimCalculated = false;
    cvector rangleUpperP13dim, rangleUpperM13dim;
    bool rangleLowerP13dimCalculated = false, rangleLowerM13dimCalculated = false;
    cvector rangleLowerP13dim, rangleLowerM13dim;
    //langle
    bool langleUpperP13dimCalculated = false, langleUpperM13dimCalculated = false;
    cvector langleUpperP13dim, langleUpperM13dim;
    bool langleLowerP13dimCalculated = false, langleLowerM13dimCalculated = false;
    cvector langleLowerP13dim, langleLowerM13dim;
    //lsquare
    bool lsquareUpperP13dimCalculated = false, lsquareUpperM13dimCalculated = false;
    cvector lsquareUpperP13dim, lsquareUpperM13dim;
    bool lsquareLowerP13dimCalculated = false, lsquareLowerM13dimCalculated = false;
    cvector lsquareLowerP13dim, lsquareLowerM13dim;
    //rsquare
    bool rsquareUpperP13dimCalculated = false, rsquareUpperM13dimCalculated = false;
    cvector rsquareUpperP13dim, rsquareUpperM13dim;
    bool rsquareLowerP13dimCalculated = false, rsquareLowerM13dimCalculated = false;
    cvector rsquareLowerP13dim, rsquareLowerM13dim;


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
    cmatrix umat(const int& dim);//Upper
    cmatrix lmat(const int& dim);//Lower

    //Spinors
    //Massless
    cvector rangle(const int& dim);
    cvector lsquare(const int& dim);
    cvector langle(const int& dim);
    cvector rsquare(const int& dim);

    //Massive
    cvector rangle(const int& spin2, const int& dim);
    cvector rangle(const int& spin2, const bool& upper, const int& dim);
    cvector langle(const int& spin2, const int& dim);
    cvector langle(const int& spin2, const bool& upper, const int& dim);
    cvector lsquare(const int& spin2, const int& dim);
    cvector lsquare(const int& spin2, const bool& upper, const int& dim);
    cvector rsquare(const int& spin2, const int& dim);
    cvector rsquare(const int& spin2, const bool& upper, const int& dim);

    //Massive cmatrix form with both spins
    cmatrix rangle_matrix(const int& dim);
    cmatrix rangle_matrix(const bool& upper, const int& dim);
    cmatrix langle_matrix(const int& dim);
    cmatrix langle_matrix(const bool& upper, const int& dim);
    
    cmatrix lsquare_matrix(const int& dim);
    cmatrix lsquare_matrix(const bool& upper, const int& dim);
    cmatrix rsquare_matrix(const int& dim);
    cmatrix rsquare_matrix(const bool& upper, const int& dim);


    //Error message:
    void usage(const char* msg) const;

    //Generators of Lorentz Rotations with no dot (acting on left chiral angle brackets)
    cmatrix lorentz_j3_lu(const int& dim) const;
    cmatrix lorentz_jp_lu(const int& dim) const;
    cmatrix lorentz_jm_lu(const int& dim) const;

    cmatrix lorentz_j3_ul(const int& dim) const;
    cmatrix lorentz_jp_ul(const int& dim) const;
    cmatrix lorentz_jm_ul(const int& dim) const;
    
    //Generators of Lorentz Rotations with a dot (acting on right chiral square brackets)
    cmatrix lorentz_j3_lu_dot(const int& dim) const;
    cmatrix lorentz_jp_lu_dot(const int& dim) const;
    cmatrix lorentz_jm_lu_dot(const int& dim) const;

    cmatrix lorentz_j3_ul_dot(const int& dim) const;
    cmatrix lorentz_jp_ul_dot(const int& dim) const;
    cmatrix lorentz_jm_ul_dot(const int& dim) const;
    
    //Generators of spin
    cmatrix spin_j3_lu(const int& dim) const;
    cmatrix spin_jp_lu(const int& dim) const;
    cmatrix spin_jm_lu(const int& dim) const;
    
    cmatrix spin_j3_ul(const int& dim) const;
    cmatrix spin_jp_ul(const int& dim) const;
    cmatrix spin_jm_ul(const int& dim) const;
    

  };

}


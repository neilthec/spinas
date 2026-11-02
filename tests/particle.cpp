


//File:  SPINAS/tests/particle.cpp

#include "test_utils.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "types.h"
//#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"


constexpr int NTESTS = 100;

using namespace spinas;


bool test_cdouble_equality(const cdouble& z1, const cdouble& z2) {
  const ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  cdouble z3 = z1 - z2;
  if (std::abs(z3) > epsilon) return false;
  return true;
}


SPINAS_TEST_SUITE(particle_tests)

// Read every cached representation, so later updates must invalidate all of it.
static std::vector<cdouble> cached_components(particle& p) {
  std::vector<cdouble> values;
  const auto vector = [&](const cvector& v) {
    for(int i=0;i<2;++i) values.push_back(v.get(i));
  };
  const auto matrix = [&](const cmatrix& m) {
    for(int i=0;i<2;++i)
      for(int j=0;j<2;++j) values.push_back(m.get(i,j));
  };
  matrix(p.umat()); matrix(p.lmat());
  if(p.mass()==0) {
    vector(p.rangle()); vector(p.langle());
    vector(p.rsquare()); vector(p.lsquare());
  } else {
    for(const bool upper : {false,true}) {
      for(const int spin : {-1,1}) {
        vector(p.rangle(spin,upper)); vector(p.langle(spin,upper));
        vector(p.rsquare(spin,upper)); vector(p.lsquare(spin,upper));
      }
      matrix(p.rangle_matrix(upper)); matrix(p.langle_matrix(upper));
      matrix(p.rsquare_matrix(upper)); matrix(p.lsquare_matrix(upper));
    }
  }
  return values;
}

static void same_cached_components(particle& actual, particle& expected) {
  const std::vector<cdouble> a = cached_components(actual);
  const std::vector<cdouble> b = cached_components(expected);
  SPINAS_CHECK_EQUAL(a.size(),b.size());
  for(std::size_t i=0;i<std::min(a.size(),b.size());++i)
    SPINAS_CHECK_EQUAL(a[i],b[i]);
}

SPINAS_TEST_CASE(cached_state_updates_and_independent_copies) {
  const ldouble originalMomentum[4] = {std::sqrt(ldouble(18)),1,2,3};
  particle reused(originalMomentum,2), original(originalMomentum,2);
  same_cached_components(reused,original);
  particle copied(reused), assigned;
  assigned = reused;

  // Exercise ordinary updates and both directions of the massless transition.
  const ldouble masses[] = {2,5,0,0,3};
  for(int step=0;step<5;++step) {
    const ldouble mass = masses[step];
    ldouble momentum[4] = {0,ldouble(0.75L)*(step+1),
      ldouble(-0.5L)*(step+2),ldouble(0.25L)*(2*step+1)};
    momentum[0] = std::sqrt(mass*mass + momentum[1]*momentum[1] +
      momentum[2]*momentum[2] + momentum[3]*momentum[3]);
    reused.set_mass(mass);
    reused.set_momentum(momentum);
    particle fresh(momentum,mass);
    same_cached_components(reused,fresh);
    same_cached_components(copied,original);
    same_cached_components(assigned,original);
  }

  // Updating a copied particle must not alter either the source or another copy.
  const ldouble copyMomentum[4] = {5,-2,4,-1};
  copied.set_momentum(copyMomentum);
  particle freshCopy(copyMomentum,2);
  same_cached_components(copied,freshCopy);
  same_cached_components(assigned,original);
  const std::vector<cdouble> before = cached_components(reused);
  assigned = copied;
  same_cached_components(assigned,freshCopy);
  const std::vector<cdouble> after = cached_components(reused);
  SPINAS_CHECK_EQUAL(before.size(),after.size());
  for(std::size_t i=0;i<before.size();++i)
    SPINAS_CHECK_EQUAL(before[i],after[i]);
}

SPINAS_TEST_CASE(invalid_momenta_and_recovery_after_rejection) {
  const ldouble nan = std::numeric_limits<ldouble>::quiet_NaN();
  const ldouble infinity = std::numeric_limits<ldouble>::infinity();
  const ldouble valid[4] = {std::sqrt(ldouble(18)),1,2,3};
  for(const ldouble invalid : {nan,infinity})
    for(int component=0;component<4;++component) {
      ldouble bad[4] = {valid[0],valid[1],valid[2],valid[3]};
      bad[component] = invalid;
      SPINAS_CHECK_THROW(particle(bad,ldouble(2)),std::runtime_error);
    }
  particle reused(valid,2);
  cached_components(reused);
  const ldouble offShell[4] = {1,2,3,4};
  SPINAS_CHECK_THROW(reused.set_momentum(offShell),std::runtime_error);
  const ldouble next[4] = {5,-2,4,-1};
  reused.set_momentum(next);
  particle fresh(next,2);
  same_cached_components(reused,fresh);
}

SPINAS_TEST_CASE(lorentz_ladder_generators_have_finite_signed_phases) {
  const ldouble tolerance = 128*std::numeric_limits<ldouble>::epsilon();
  const cmatrix shapes[] = {cmatrix(1,1,-1,-1),cmatrix(1,-1,1,-1),
    cmatrix(-1,1,-1,1),cmatrix(-1,-1,1,1)};
  for(const ldouble mass : {ldouble(0),ldouble(1)}) {
    const ldouble energy = std::sqrt(mass*mass+1);
    const ldouble alongX[4] = {energy,1,0,0};
    particle p(alongX,mass);
    const ldouble small = mass==0 ? 0 : (energy-1)/(2*mass);
    const ldouble large = mass==0 ? 1 : (energy+1)/(2*mass);
    // At theta=pi/2 and phi=0, the generator entries are real and explicit.
    const cmatrix expected[] = {small*shapes[0],large*shapes[1],
      small*shapes[2],large*shapes[3],large*shapes[2],small*shapes[3],
      large*shapes[0],small*shapes[1]};
    const cmatrix actual[] = {p.lorentz_jp_lu(),p.lorentz_jm_lu(),
      p.lorentz_jp_ul(),p.lorentz_jm_ul(),p.lorentz_jp_lu_dot(),
      p.lorentz_jm_lu_dot(),p.lorentz_jp_ul_dot(),p.lorentz_jm_ul_dot()};
    for(int generator=0;generator<8;++generator)
      for(int row=0;row<2;++row)
        for(int column=0;column<2;++column)
          SPINAS_CHECK_SMALL(std::abs(actual[generator].get(row,column)-
            expected[generator].get(row,column)),tolerance);

    // Both z axes and all azimuthal quadrants must also avoid invalid radii.
    const ldouble directions[][3] = {{0,0,1},{0,0,-1},{1,0,0},
      {-1,0,0},{0,1,0},{0,-1,0}};
    for(const auto& direction : directions) {
      const ldouble momentum[4] = {energy,direction[0],direction[1],direction[2]};
      p.set_momentum(momentum);
      const cmatrix matrices[] = {p.lorentz_j3_lu(),p.lorentz_j3_ul(),
        p.lorentz_j3_lu_dot(),p.lorentz_j3_ul_dot(),p.lorentz_jp_lu(),
        p.lorentz_jm_lu(),p.lorentz_jp_ul(),p.lorentz_jm_ul(),
        p.lorentz_jp_lu_dot(),p.lorentz_jm_lu_dot(),
        p.lorentz_jp_ul_dot(),p.lorentz_jm_ul_dot()};
      for(const cmatrix& matrix : matrices)
        for(int row=0;row<2;++row)
          for(int column=0;column<2;++column) {
            const cdouble value = matrix.get(row,column);
            SPINAS_CHECK(std::isfinite(value.real()) && std::isfinite(value.imag()));
          }
    }
  }
}

SPINAS_TEST_CASE(constructor_test) {
  SPINAS_TEST_MESSAGE("Testing particle: ");
  SPINAS_TEST_MESSAGE("\t* Constructor");
  ldouble mom1[4], mom2[4], mass2;
  particle p1;

  for(int i=0;i<NTESTS;i++){
    //Choose Mom1
    choose_random_momentum(mom1, -50, 50);
    mass2 = mom1[0]*mom1[0];
    for(int j=1;j<4;j++) mass2 -= mom1[j]*mom1[j];
    p1 = particle(mom1, std::sqrt(mass2));
    SPINAS_CHECK(p1.test_angles());
    
    //Choose Mom2 with same mass
    choose_random_momentum(mom2, -50, 50);
    mom2[0] = mass2;
    for(int j=1;j<4;j++) mom2[0] += mom2[j]*mom2[j];
    mom2[0] = std::sqrt(mom2[0]);
    p1.set_momentum(mom2);
    SPINAS_CHECK(p1.test_angles());
  }
}

// A boosted light particle must retain its declared mass in spinor products.
// Rounding external momenta to double is normal when exchanging phase-space
// points with another numerical code; subtracting E-|p| loses this information.
SPINAS_TEST_CASE(boosted_spinor_mass_normalization) {
  const ldouble mass = ldouble(0.0005L);
  const ldouble relativeTolerance = 64*std::numeric_limits<ldouble>::epsilon();
  const ldouble energies[] = {mass, ldouble(1), ldouble(50), ldouble(200)};
  for(const ldouble energy : energies){
    const ldouble magnitude = static_cast<ldouble>(
      std::sqrt(static_cast<double>(energy*energy-mass*mass)));
    for(int direction=0; direction<2; ++direction){
      const ldouble momenta[4] = {energy,
        direction==0 ? ldouble(0) : static_cast<ldouble>(static_cast<double>(ldouble(0.6L)*magnitude)),
        ldouble(0),
        direction==0 ? magnitude : static_cast<ldouble>(static_cast<double>(ldouble(0.8L)*magnitude))};
      particle boosted(momenta,mass);
      // <pp>^{-+}=-m and [pp]^{-+}=m at every boost and direction.
      SPINAS_CHECK_SMALL(std::abs((boosted.langle(-1)*boosted.rangle(1))/mass+ldouble(1)), relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs((boosted.lsquare(-1)*boosted.rsquare(1))/mass-ldouble(1)), relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs((boosted.langle(1)*boosted.rangle(-1))/mass-ldouble(1)), relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs((boosted.lsquare(1)*boosted.rsquare(-1))/mass+ldouble(1)), relativeTolerance);
      if(direction==0){
        // Raising/lowering must use the same stable boost ratio as the spinors.
        const cmatrix angleActual = boosted.lorentz_jm_lu()*boosted.rangle_matrix();
        const cmatrix angleExpected = boosted.rangle_matrix()*boosted.spin_jm_lu();
        const cmatrix squareActual = boosted.lorentz_jp_lu_dot()*boosted.lsquare_matrix();
        const cmatrix squareExpected = boosted.lsquare_matrix()*boosted.spin_jp_lu();
        for(int row=0; row<2; ++row)
          for(int column=0; column<2; ++column){
            const ldouble angleScale = std::max(ldouble(1),std::abs(angleExpected.get(row,column)));
            const ldouble squareScale = std::max(ldouble(1),std::abs(squareExpected.get(row,column)));
            SPINAS_CHECK_SMALL(std::abs(angleActual.get(row,column)-angleExpected.get(row,column))/angleScale,
              relativeTolerance);
            SPINAS_CHECK_SMALL(std::abs(squareActual.get(row,column)-squareExpected.get(row,column))/squareScale,
              relativeTolerance);
          }
      }
    }
  }
}

SPINAS_TEST_CASE(stable_spinor_zero_mass_and_mass_sign) {
  const ldouble zero[4] = {0,0,0,0};
  particle zeroParticle(zero,ldouble(0));
  for(int component=0; component<2; ++component){
    SPINAS_CHECK_EQUAL(zeroParticle.rangle().get(component),cdouble(0,0));
    SPINAS_CHECK_EQUAL(zeroParticle.lsquare().get(component),cdouble(0,0));
  }
  const ldouble photon[4] = {50,0,0,50};
  particle massless(photon,ldouble(0));
  for(int component=0; component<2; ++component){
    SPINAS_CHECK(std::isfinite(std::real(massless.rangle().get(component))));
    SPINAS_CHECK(std::isfinite(std::imag(massless.lsquare().get(component))));
  }
  const ldouble mass = ldouble(0.0005L);
  const ldouble massive[4] = {50,0,0,static_cast<ldouble>(49.9999999975)};
  particle positiveMass(massive,mass), negativeMass(massive,-mass);
  for(int spin=-1; spin<=1; spin+=2){
    SPINAS_CHECK_EQUAL(positiveMass.rangle(spin),negativeMass.rangle(spin));
    SPINAS_CHECK_EQUAL(positiveMass.lsquare(spin),negativeMass.lsquare(spin));
  }
  const ldouble offShell[4] = {2,0,0,0};
  SPINAS_CHECK_THROW(particle(offShell,ldouble(1)),std::runtime_error);
}

SPINAS_TEST_CASE(extreme_boost_shell_and_lightcone_matrices) {
  const ldouble mass = 1;
  const ldouble relativeTolerance = 128*std::numeric_limits<ldouble>::epsilon();
  for(const ldouble energy : {ldouble(1),ldouble(1000),ldouble(1e120L)}){
    const ldouble magnitude = std::sqrt(energy-mass)*std::sqrt(energy+mass);
    for(const int direction : {-1,1}){
      const ldouble momentum[4] = {energy,0,0,direction*magnitude};
      particle p(momentum,mass);
      const ldouble large = energy+magnitude;
      const ldouble small = mass*(mass/energy)/(1+magnitude/energy);
      const cmatrix upper = p.umat(), lower = p.lmat();
      const int smallUpper = direction>0 ? 0 : 1;
      SPINAS_CHECK_SMALL(std::abs(upper.get(smallUpper,smallUpper)/small-cdouble(1,0)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs(lower.get(1-smallUpper,1-smallUpper)/small-cdouble(1,0)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs(upper.get(1-smallUpper,1-smallUpper)/large-cdouble(1,0)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs(lower.get(smallUpper,smallUpper)/large-cdouble(1,0)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs((p.langle(-1)*p.rangle(1))/mass+ldouble(1)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs((p.lsquare(-1)*p.rsquare(1))/mass-ldouble(1)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs(upper.get_det()/mass/mass-ldouble(1)),relativeTolerance);
      SPINAS_CHECK_SMALL(std::abs(lower.get_det()/mass/mass-ldouble(1)),relativeTolerance);
      const ldouble materiallyOffShell[4] = {energy*ldouble(1.000001L),0,0,direction*magnitude};
      SPINAS_CHECK_THROW(particle(materiallyOffShell,mass),std::runtime_error);
    }
  }
}

SPINAS_TEST_CASE(stable_transverse_mass_and_exact_massless_axes) {
  const ldouble mass = ldouble(0.0005L), energy = 1000;
  const ldouble transverse = 2*mass;
  const ldouble longitudinal = std::sqrt(energy*energy-mass*mass-transverse*transverse);
  const ldouble tolerance = 128*std::numeric_limits<ldouble>::epsilon();
  for(const int direction : {-1,1}){
    const ldouble momentum[4] = {energy,transverse,0,direction*longitudinal};
    particle p(momentum,mass), signedMass(momentum,-mass);
    const ldouble small = (mass*mass+transverse*transverse)/(energy+longitudinal);
    const int index = direction>0 ? 0 : 1;
    SPINAS_CHECK_SMALL(std::abs(p.umat().get(index,index)/small-cdouble(1,0)),tolerance);
    SPINAS_CHECK_EQUAL(p.umat(),signedMass.umat());
    SPINAS_CHECK_EQUAL(p.lmat(),signedMass.lmat());
    const ldouble photonMomentum[4] = {energy,0,0,direction*energy};
    particle photon(photonMomentum,0);
    SPINAS_CHECK_EQUAL(photon.umat().get(index,index),cdouble(0,0));
    SPINAS_CHECK_EQUAL(photon.rangle().get(direction>0 ? 1 : 0),cdouble(0,0));
  }
}

SPINAS_TEST_CASE(finite_spinor_roots_when_energy_sum_overflows) {
  const ldouble energy = std::numeric_limits<ldouble>::max()/ldouble(1.5);
  const ldouble mass = 1;
  const ldouble tolerance = 128*std::numeric_limits<ldouble>::epsilon();
  SPINAS_CHECK(std::isinf(energy+energy));
  for(const int direction : {-1,1}){
    const ldouble momentum[4] = {energy,0,0,direction*energy};
    particle p(momentum,mass);
    for(const int spin : {-1,1})
      for(int component=0; component<2; component++){
        SPINAS_CHECK(std::isfinite(std::real(p.rangle(spin).get(component))));
        SPINAS_CHECK(std::isfinite(std::imag(p.lsquare(spin).get(component))));
      }
    SPINAS_CHECK_SMALL(std::abs((p.langle(-1)*p.rangle(1))/mass+ldouble(1)),tolerance);
    SPINAS_CHECK_SMALL(std::abs((p.lsquare(-1)*p.rsquare(1))/mass-ldouble(1)),tolerance);
  }
}

SPINAS_TEST_CASE(signed_energy_and_off_shell_matrix_conventions) {
  const ldouble energy = 1000, mass = 1;
  const ldouble magnitude = std::sqrt(energy-mass)*std::sqrt(energy+mass);
  const ldouble negative[4] = {-energy,0,0,magnitude};
  particle signedParticle(negative,mass);
  SPINAS_CHECK_EQUAL(signedParticle.umat().get(0,0),cdouble(negative[0]-negative[3],0));
  SPINAS_CHECK_EQUAL(signedParticle.umat().get(1,1),cdouble(negative[0]+negative[3],0));
  SPINAS_CHECK_EQUAL(signedParticle.lmat().get(0,0),cdouble(negative[0]+negative[3],0));
  SPINAS_CHECK(std::isnan(std::real(signedParticle.rangle(-1).get(0))));
  const ldouble virtualMomentum[4] = {2,0,0,1};
  const cmatrix genericUpper(virtualMomentum,true), genericLower(virtualMomentum,false);
  SPINAS_CHECK_EQUAL(genericUpper.get(0,0),cdouble(1,0));
  SPINAS_CHECK_EQUAL(genericLower.get(1,1),cdouble(1,0));
  SPINAS_CHECK_EQUAL(genericUpper.get_det(),ldouble(3));
  SPINAS_CHECK_THROW(particle(virtualMomentum,ldouble(1)),std::runtime_error);
}

SPINAS_TEST_CASE(stable_on_shell_scalar_products) {
  // Analytic edge cases remain valid even when E - |p| is not representable.
  const ldouble tolerance = 256*std::numeric_limits<ldouble>::epsilon();
  const ldouble huge = std::min(std::numeric_limits<ldouble>::max()/ldouble(2),ldouble(1e300L));
  const ldouble same[4] = {huge,0,0,huge};
  particle sameInstance(same,1);
  SPINAS_CHECK_SMALL(std::abs(sameInstance.dot(sameInstance)-ldouble(1)),tolerance);
  const ldouble maximumEnergy = std::numeric_limits<ldouble>::max()/ldouble(2);
  const ldouble maximumMomentum[4] = {maximumEnergy,0,0,maximumEnergy};
  particle maximumFirst(maximumMomentum,1), maximumSecond(maximumMomentum,1);
  SPINAS_CHECK_SMALL(std::abs(maximumFirst.dot(maximumSecond)-ldouble(1)),tolerance);
  const ldouble zero[4] = {0,0,0,0};
  particle nullParticle(zero,0);
  SPINAS_CHECK_EQUAL(nullParticle.dot(sameInstance),ldouble(0));
  particle masslessFirst(same,0), masslessSecond(same,0);
  SPINAS_CHECK_EQUAL(masslessFirst.dot(masslessSecond),ldouble(0));
  const ldouble temporarilyUnchecked[4] = {2,0,0,1};
  particle changedMass(temporarilyUnchecked,std::sqrt(ldouble(3)));
  changedMass.set_mass(ldouble(0.25));
  SPINAS_CHECK_EQUAL(changedMass.dot(changedMass),ldouble(3));
  const ldouble past[4] = {-5,1,2,3}, future[4] = {5,-1,-2,-3};
  particle negativeEnergy(past,std::sqrt(ldouble(11))), positiveEnergy(future,std::sqrt(ldouble(11)));
  SPINAS_CHECK_EQUAL(negativeEnergy.dot(positiveEnergy),ldouble(-11));
}

SPINAS_TEST_CASE(momentum_matrices_test) {
  SPINAS_TEST_MESSAGE("\t* Momentum Matrices");
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  int i = 0;
  cdouble zero = cdouble(0, 0), one = cdouble(1, 0);
  ldouble mom1[4] = {13, -4, 7, -9}, mom3[4] = {6, 1, -2, std::sqrt(8)};
  particle p1 = particle(mom3, std::sqrt(23));
  p1.set_momentum(mom1);
  SPINAS_CHECK_EQUAL(p1.umat(), cmatrix(mom1, true));
  SPINAS_CHECK_EQUAL(p1.lmat(), cmatrix(mom1, false));
  
  // Determinant
  cdouble p1detm23 = p1.umat().get_det() - cdouble(23, 0);
  SPINAS_CHECK_SMALL(std::abs(p1detm23), epsilon);
  SPINAS_CHECK_EQUAL(p1.lmat(), cmatrix(mom1, false));
  
  p1detm23 = p1.lmat().get_det() - cdouble(23, 0);
  SPINAS_CHECK_SMALL(std::abs(p1detm23), epsilon);
  
  // Anticommutation: p1*p2 + p2*p1 = 2p1.p2
  ldouble mom2[4] = {5, 1, -2, 3};
  particle p2 = particle(std::sqrt(25 - 1 - 4 - 9));
  p2.set_momentum(mom2);
  cmatrix id = cmatrix(one, zero, zero, one);
  cmatrix anticom = p1.umat() * p2.lmat() + p2.umat() * p1.lmat();
  SPINAS_CHECK_EQUAL(anticom, 2 * p1.dot(p2) * id);
  anticom = p1.lmat() * p2.umat() + p2.lmat() * p1.umat();
  SPINAS_CHECK_EQUAL(anticom,  2 * p1.dot(p2) * id);
}

SPINAS_TEST_CASE(helicity_spinors_test) {
  SPINAS_TEST_MESSAGE("\t* Helicity Spinors");
  ldouble mom1[4] = {std::sqrt(16 + 49 + 81), -4, 7, -9};
  particle p1 = particle(mom1, 0);

  for(int i=0;i<100;i++){
    
    choose_random_massless_momentum(mom1,-50,50);
    p1.set_momentum(mom1);
    
    //|1> == [1|*
    SPINAS_CHECK_EQUAL(p1.rangle(), p1.lsquare().get_conjugate());
    
    //<1| == |1]*
    SPINAS_CHECK_EQUAL(p1.langle(), p1.rsquare().get_conjugate());
    
    // p1|1> == 0
    SPINAS_CHECK(p1.umat() * p1.rangle() == cvector(cdouble(0, 0), cdouble(0, 0)));
    
    // [1|p1 == 0
    SPINAS_CHECK(p1.lsquare() * p1.umat() == cvector(cdouble(0, 0), cdouble(0, 0)));
    
    // <11> == 0
    SPINAS_CHECK(p1.langle() * p1.rangle() == cdouble(0, 0));
    
    // [11] == 0
    SPINAS_CHECK(p1.lsquare() * p1.rsquare() == cdouble(0, 0));
    
    // |1>[1| = p1
    SPINAS_CHECK(outer(p1.rangle(), p1.lsquare()) == p1.lmat());
    
    // |1]<1| = p1
    SPINAS_CHECK(outer(p1.rsquare(), p1.langle()) == p1.umat());
    
    // p1|1] == 0
    SPINAS_CHECK(p1.lmat() * p1.rsquare() == cvector(cdouble(0, 0), cdouble(0, 0)));
    
    // <1|p1 == 0
    SPINAS_CHECK(p1.langle() * p1.lmat() == cvector(cdouble(0, 0), cdouble(0, 0)));
    
    //Lorentz and Helicity Algebra
    constexpr ldouble half=0.5;
    cvector zero_vector=cvector(0,0);
    
    // LJ3 |1> == -1/2 |1>
    SPINAS_CHECK(p1.lorentz_j3_lu()*p1.rangle() == -half*p1.rangle());
    // LJ+ |1> == 0
    SPINAS_CHECK(p1.lorentz_jp_lu()*p1.rangle() == zero_vector);
    // LJ- |1> == 0
    SPINAS_CHECK(p1.lorentz_jm_lu()*p1.rangle() == zero_vector);
    // LJ3 <1| == -1/2 <1|
    SPINAS_CHECK(p1.lorentz_j3_ul()*p1.langle() == -half*p1.langle());
    // LJ+ <1| == 0
    SPINAS_CHECK(p1.lorentz_jp_ul()*p1.langle() == zero_vector);
    // LJ- <1| == 0
    SPINAS_CHECK(p1.lorentz_jm_ul()*p1.langle() == zero_vector);
    
    // LJ3 [1| == 1/2 [1|
    SPINAS_CHECK(p1.lorentz_j3_lu_dot()*p1.lsquare() == half*p1.lsquare());
    // LJ+ [1| == 0
    SPINAS_CHECK(p1.lorentz_jp_lu_dot()*p1.lsquare() == zero_vector);
    // LJ- [1| == 0
    SPINAS_CHECK(p1.lorentz_jm_lu_dot()*p1.lsquare() == zero_vector);
    // LJ3 |1] == 1/2 |1]
    SPINAS_CHECK(p1.lorentz_j3_ul_dot()*p1.rsquare() == half*p1.rsquare());
    // LJ+ |1] == 0
    SPINAS_CHECK(p1.lorentz_jp_ul_dot()*p1.rsquare() == zero_vector);
    // LJ- |1] == 0
    SPINAS_CHECK(p1.lorentz_jm_ul_dot()*p1.rsquare() == zero_vector);
    
  }
}


SPINAS_TEST_CASE(spin_spinors) {
  SPINAS_TEST_MESSAGE("\t* Spin Spinors");
  ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
  ldouble mom1[4] = {std::sqrt(100 + 16 + 49 + 81), -4, 7, -9};
  particle p1 = particle(mom1, 10);
  ldouble mass = 10;
  cmatrix epsU = cmatrix(cdouble(0, 0), cdouble(1, 0), cdouble(-1, 0), cdouble(0, 0));
  cmatrix epsL = cmatrix(cdouble(0,0),cdouble(-1,0),cdouble(1,0),cdouble(0,0));
  cmatrix delta = cmatrix(cdouble(1, 0), cdouble(0, 0), cdouble(0, 0), cdouble(1, 0));

  for(int i=0;i<100;i++){
    
    mass = choose_random_momentum(mom1,-50,50);
    p1.set_mass(mass);
    p1.set_momentum(mom1);
    
    // Conjugation
    for (int j = -1; j < 2; j = j + 2){
      //|1>^I == ([1|_I)*
      SPINAS_CHECK_EQUAL(p1.rangle(j), p1.lsquare(j, false).get_conjugate());
      
      //<1|^I == (|1]_I)*
      SPINAS_CHECK_EQUAL(p1.langle(j), p1.rsquare(j, false).get_conjugate());
      
      //|1>_I == -([1|^I)*
      SPINAS_CHECK_EQUAL(p1.rangle(j, false), -p1.lsquare(j).get_conjugate());
      
      //<1|_I == -(|1]^I)*
      SPINAS_CHECK_EQUAL(p1.langle(j, false), -p1.rsquare(j).get_conjugate());
    }
    
    
    
    // Spinor Products
    for(int i=-1;i<2;i=i+2)
      for(int j=-1;j<2;j=j+2){
	int epsi = 0, epsj=0;
	if(i==1) epsi=1;
	if(j==1) epsj=1;
	
	//<11>^{IJ} == -m eps^{IJ}
	SPINAS_CHECK_SMALL(std::abs(p1.langle(i)*p1.rangle(j) - (-mass*epsU.get(epsi,epsj))), epsilon);
	
	//<11>^I_J == -m delta^I_J
	SPINAS_CHECK_SMALL(std::abs(p1.langle(i)*p1.rangle(j,false) - (-mass*delta.get(epsi,epsj))), epsilon);
	
	//<11>_I^J == m delta_I^J
	SPINAS_CHECK_SMALL(std::abs(p1.langle(i,false)*p1.rangle(j) - mass*delta.get(epsi,epsj)), epsilon);
	
	//<11>_{IJ} == m eps_{IJ}
	SPINAS_CHECK_SMALL(std::abs(p1.langle(i,false)*p1.rangle(j,false) - mass*epsL.get(epsi,epsj)), epsilon);
	
	//[11]_{IJ} == -m eps_{IJ}
	SPINAS_CHECK_SMALL(std::abs(p1.lsquare(i,false)*p1.rsquare(j,false) - (-mass*epsL.get(epsi,epsj))), epsilon);
	
	//[11]_I^J == -m delta_I^J
	SPINAS_CHECK_SMALL(std::abs(p1.lsquare(i,false)*p1.rsquare(j) - (-mass*delta.get(epsi,epsj))), epsilon);
	
	//[11]^I_J == m delta^I_J
	SPINAS_CHECK_SMALL(std::abs(p1.lsquare(i)*p1.rsquare(j,false) - mass*delta.get(epsi,epsj)), epsilon);
	
	//[11]^{IJ} == m eps^{IJ}
	SPINAS_CHECK_SMALL(std::abs(p1.lsquare(i)*p1.rsquare(j) - mass*epsU.get(epsi,epsj)), epsilon);
	
      }
    
    
    //Outer Products
    //|1>^I[1|_I = p1
    cmatrix outerProd = outer(p1.rangle(-1), p1.lsquare(-1, false)) +
      outer(p1.rangle(1), p1.lsquare(1, false));
    SPINAS_CHECK(outerProd == p1.lmat());
    
    //|1>_I[1|^I = -p1
    outerProd = outer(p1.rangle(-1, false), p1.lsquare(-1)) +
      outer(p1.rangle(1, false), p1.lsquare(1));
    SPINAS_CHECK(outerProd == -p1.lmat());
    
    //|1]_I<1|^I = p1
    outerProd = outer(p1.rsquare(-1, false), p1.langle(-1)) +
      outer(p1.rsquare(1, false), p1.langle(1));
    SPINAS_CHECK(outerProd == p1.umat());
    
    //|1]^I<1|_I = -p1
    outerProd = outer(p1.rsquare(-1), p1.langle(-1, false)) +
      outer(p1.rsquare(1), p1.langle(1, false));
    SPINAS_CHECK(outerProd == -p1.umat());
    
    // Momentum times Spinor Gives Mass times Spinor
    for(int j=-1;j<2;j=j+2){
      //p1|1>^I=-m|1]^I
      SPINAS_CHECK_EQUAL(p1.umat()*p1.rangle(j), -mass*p1.rsquare(j));
      
      //<1|^I*p1=[1|^I*m
      SPINAS_CHECK_EQUAL(p1.langle(j)*p1.lmat(), mass*p1.lsquare(j));
      
      //p1|1>_I=-m|1]_I
      SPINAS_CHECK_EQUAL(p1.umat()*p1.rangle(j,false), -mass*p1.rsquare(j,false));
      
      //<1|_I*p1=[1|_I*m
      SPINAS_CHECK_EQUAL(p1.langle(j,false)*p1.lmat(), mass*p1.lsquare(j,false));
      
      //p1|1]^I=-m|1>^I
      SPINAS_CHECK_EQUAL(p1.lmat()*p1.rsquare(j), -mass*p1.rangle(j));
      
      //[1|^I*p1=<1|^I*m
      SPINAS_CHECK_EQUAL(p1.lsquare(j)*p1.umat(), mass*p1.langle(j));
      
      //p1|1]_I=-m|1>_I
      SPINAS_CHECK_EQUAL(p1.lmat()*p1.rsquare(j,false), -mass*p1.rangle(j,false));
      
      //[1|_I*p1=<1|_I*m
      SPINAS_CHECK_EQUAL(p1.lsquare(j,false)*p1.umat(), mass*p1.langle(j,false));
    }
    
    
    
    //Lorentz and Spin Algebra
    
    // LJ3 |1>^I == |1>^J SJ3_J^I
    SPINAS_CHECK(p1.lorentz_j3_lu()*p1.rangle_matrix() == p1.rangle_matrix()*p1.spin_j3_lu());
    // LJ+ |1>^I == |1>^J SJ+_J^I
    SPINAS_CHECK(p1.lorentz_jp_lu()*p1.rangle_matrix() == p1.rangle_matrix()*p1.spin_jp_lu());
    // LJ- |1>^I == |1>^J SJ-_J^I
    SPINAS_CHECK(p1.lorentz_jm_lu()*p1.rangle_matrix() == p1.rangle_matrix()*p1.spin_jm_lu());
    
    // LJ3 <1|^I == <1|^J SJ3_J^I
    SPINAS_CHECK(p1.lorentz_j3_ul()*p1.langle_matrix() == p1.langle_matrix()*p1.spin_j3_lu());
    // LJ+ <1|^I == <1|^J SJ+_J^I
    SPINAS_CHECK(p1.lorentz_jp_ul()*p1.langle_matrix() == p1.langle_matrix()*p1.spin_jp_lu());
    // LJ- <1|^I == <1|^J SJ-_J^I
    SPINAS_CHECK(p1.lorentz_jm_ul()*p1.langle_matrix() == p1.langle_matrix()*p1.spin_jm_lu());
    
    // LJ3 |1>_I == |1>_J SJ3^J_I
    SPINAS_CHECK(p1.lorentz_j3_lu()*p1.rangle_matrix(LOWER) == p1.rangle_matrix(LOWER)*p1.spin_j3_ul());
    // LJ+ |1>_I == |1>_J SJ+^J_I
    SPINAS_CHECK(p1.lorentz_jp_lu()*p1.rangle_matrix(LOWER) == p1.rangle_matrix(LOWER)*p1.spin_jp_ul());
    // LJ- |1>_I == |1>_J SJ-^J_I
    SPINAS_CHECK(p1.lorentz_jm_lu()*p1.rangle_matrix(LOWER) == p1.rangle_matrix(LOWER)*p1.spin_jm_ul());
    
    // LJ3 <1|_I == <1|_J SJ3^J_I
    SPINAS_CHECK(p1.lorentz_j3_ul()*p1.langle_matrix(LOWER) == p1.langle_matrix(LOWER)*p1.spin_j3_ul());
    // LJ+ <1|_I == <1|_J SJ+^J_I
    SPINAS_CHECK(p1.lorentz_jp_ul()*p1.langle_matrix(LOWER) == p1.langle_matrix(LOWER)*p1.spin_jp_ul());
    // LJ- <1|_I == <1|_J SJ-^J_I
    SPINAS_CHECK(p1.lorentz_jm_ul()*p1.langle_matrix(LOWER) == p1.langle_matrix(LOWER)*p1.spin_jm_ul());
    
    // LJ3 [1|^I == [1|^J SJ3_J^I
    SPINAS_CHECK(p1.lorentz_j3_lu_dot()*p1.lsquare_matrix() == p1.lsquare_matrix()*p1.spin_j3_lu());
    // LJ+ [1|^I == [1|^J SJ+_J^I
    SPINAS_CHECK(p1.lorentz_jp_lu_dot()*p1.lsquare_matrix() == p1.lsquare_matrix()*p1.spin_jp_lu());
    // LJ- [1|^I == [1|^J SJ-_J^I
    SPINAS_CHECK(p1.lorentz_jm_lu_dot()*p1.lsquare_matrix() == p1.lsquare_matrix()*p1.spin_jm_lu());
    
    // LJ3 |1]^I == |1]^J SJ3_J^I
    SPINAS_CHECK(p1.lorentz_j3_ul_dot()*p1.rsquare_matrix() == p1.rsquare_matrix()*p1.spin_j3_lu());
    // LJ+ |1]^I == |1]^J SJ+_J^I
    SPINAS_CHECK(p1.lorentz_jp_ul_dot()*p1.rsquare_matrix() == p1.rsquare_matrix()*p1.spin_jp_lu());
    // LJ- |1]^I == |1]^J SJ-_J^I
    SPINAS_CHECK(p1.lorentz_jm_ul_dot()*p1.rsquare_matrix() == p1.rsquare_matrix()*p1.spin_jm_lu());
    
    // LJ3 [1|_I == [1|_J SJ3^J_I
    SPINAS_CHECK(p1.lorentz_j3_lu_dot()*p1.lsquare_matrix(LOWER) == p1.lsquare_matrix(LOWER)*p1.spin_j3_ul());
    // LJ+ [1|_I == [1|_J SJ+^J_I
    SPINAS_CHECK(p1.lorentz_jp_lu_dot()*p1.lsquare_matrix(LOWER) == p1.lsquare_matrix(LOWER)*p1.spin_jp_ul());
    // LJ- [1|_I == [1|_J SJ-^J_I
    SPINAS_CHECK(p1.lorentz_jm_lu_dot()*p1.lsquare_matrix(LOWER) == p1.lsquare_matrix(LOWER)*p1.spin_jm_ul());
    
    // LJ3 |1]_I == |1]_J SJ3^J_I
    SPINAS_CHECK(p1.lorentz_j3_ul_dot()*p1.rsquare_matrix(LOWER) == p1.rsquare_matrix(LOWER)*p1.spin_j3_ul());
    // LJ+ |1]_I == |1]_J SJ+^J_I
    SPINAS_CHECK(p1.lorentz_jp_ul_dot()*p1.rsquare_matrix(LOWER) == p1.rsquare_matrix(LOWER)*p1.spin_jp_ul());
    // LJ- |1]_I == |1]_J SJ-^J_I
    SPINAS_CHECK(p1.lorentz_jm_ul_dot()*p1.rsquare_matrix(LOWER) == p1.rsquare_matrix(LOWER)*p1.spin_jm_ul());
    
  }
  
  
}



SPINAS_TEST_SUITE_END()


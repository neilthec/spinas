
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

//File:  SPINAS/source/process.cpp

#include <iostream>
#include <stdio.h>
#include <vector>

#include "types.h"
//#include "aliases.h"
#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "sproduct.h"
#include "process.h"

namespace spinas {
  //Constructors
  process::process()
  {}

  //Prototype set_momenta
  void process::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    throw std::runtime_error("set_momenta(p1, p2, p3, p4) not implemented yet!");
  }

  //Functions to determine the number of loops, the normalization factor and the spin indices for massive spin 1.
  //The individual functions are for users writing code by hand and are available up to some number of particles.
  //Above that, it is expected that the code would be output by code and would use the master functions with general numbers of particles.
  //Number of loops
  int process::get_num_spin_loops(const int& dsA){
    int dsList[1] = {dsA};
    return get_num_spin_loops(dsList,1);
  }
  int process::get_num_spin_loops(const int& dsA, const int& dsB){
    int dsList[2] = {dsA,dsB};
    return get_num_spin_loops(dsList,2);
  }
  int process::get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC){
    int dsList[3] = {dsA,dsB,dsC};
    return get_num_spin_loops(dsList,3);
  }
  int process::get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC, const int& dsD){
    int dsList[4] = {dsA,dsB,dsC,dsD};
    return get_num_spin_loops(dsList,4);
  }
  int process::get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE){
    int dsList[5] = {dsA,dsB,dsC,dsD,dsE};
    return get_num_spin_loops(dsList,5);
  }
  int process::get_num_spin_loops(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE, const int& dsF){
    int dsList[6] = {dsA,dsB,dsC,dsD,dsE,dsF};
    return get_num_spin_loops(dsList,6);
  }
  int process::get_num_spin_loops(const int* dsList, const int length){
    int nLoops = 1;
    for(int i=0;i<length;i++)
      if(dsList[i]==0) nLoops *= 2;
    return nLoops;
  }
  //Normalization Factor
  ldouble process::get_spin_normalization(const int& dsA){
    int dsList[1] = {dsA};
    return get_spin_normalization(dsList,1);
  }
  ldouble process::get_spin_normalization(const int& dsA, const int& dsB){
    int dsList[2] = {dsA,dsB};
    return get_spin_normalization(dsList,2);
  }
  ldouble process::get_spin_normalization(const int& dsA, const int& dsB, const int& dsC){
    int dsList[3] = {dsA,dsB,dsC};
    return get_spin_normalization(dsList,3);
  }
  ldouble process::get_spin_normalization(const int& dsA, const int& dsB, const int& dsC, const int& dsD){
    int dsList[4] = {dsA,dsB,dsC,dsD};
    return get_spin_normalization(dsList,4);
  }
  ldouble process::get_spin_normalization(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE){
    int dsList[5] = {dsA,dsB,dsC,dsD,dsE};
    return get_spin_normalization(dsList,5);
  }
  ldouble process::get_spin_normalization(const int& dsA, const int& dsB, const int& dsC, const int& dsD, const int& dsE, const int& dsF){
    int dsList[6] = {dsA,dsB,dsC,dsD,dsE,dsF};
    return get_spin_normalization(dsList,6);
  }
  ldouble process::get_spin_normalization(const int* dsList, const int length){
    ldouble normFactor = 1;
    //constexpr ldouble sqrt2=std::sqrt(2.0);
    for(int i=0;i<length;i++)
      if(dsList[i]==0) normFactor /= std::sqrt(2.0);
    return normFactor;
  }
  //Spinor Spins
  void process::get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& iterator){
    const int* dsList[1] = {&dsA};
    int* dsLista[1] = {&dsAa};
    int* dsListb[1] = {&dsAb};
    return get_spinor_spins(dsList, dsLista, dsListb, 1, iterator);
  }
  void process::get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& iterator){
    const int* dsList[2] = {&dsA, &dsB};
    int* dsLista[2] = {&dsAa, &dsBa};
    int* dsListb[2] = {&dsAb, &dsBb};
    return get_spinor_spins(dsList, dsLista, dsListb, 2, iterator);
  }
  void process::get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& iterator){
    const int* dsList[3] = {&dsA, &dsB, &dsC};
    int* dsLista[3] = {&dsAa, &dsBa, &dsCa};
    int* dsListb[3] = {&dsAb, &dsBb, &dsCb};
    return get_spinor_spins(dsList, dsLista, dsListb, 3, iterator);
  }
  void process::get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& dsD, int& dsDa, int& dsDb, const int& iterator){
    const int* dsList[4] = {&dsA, &dsB, &dsC, &dsD};
    int* dsLista[4] = {&dsAa, &dsBa, &dsCa, &dsDa};
    int* dsListb[4] = {&dsAb, &dsBb, &dsCb, &dsDb};
    return get_spinor_spins(dsList, dsLista, dsListb, 4, iterator);
  }
  void process::get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& dsD, int& dsDa, int& dsDb, const int& dsE, int& dsEa, int& dsEb, const int& iterator){
    const int* dsList[5] = {&dsA, &dsB, &dsC, &dsD, &dsE};
    int* dsLista[5] = {&dsAa, &dsBa, &dsCa, &dsDa, &dsEa};
    int* dsListb[5] = {&dsAb, &dsBb, &dsCb, &dsDb, &dsEb};
    return get_spinor_spins(dsList, dsLista, dsListb, 5, iterator);
  }
  void process::get_spinor_spins(const int& dsA, int& dsAa, int& dsAb, const int& dsB, int& dsBa, int& dsBb, const int& dsC, int& dsCa, int& dsCb, const int& dsD, int& dsDa, int& dsDb, const int& dsE, int& dsEa, int& dsEb, const int& dsF, int& dsFa, int& dsFb, const int& iterator){
    const int* dsList[6] = {&dsA, &dsB, &dsC, &dsD, &dsE, &dsF};
    int* dsLista[6] = {&dsAa, &dsBa, &dsCa, &dsDa, &dsEa, &dsFa};
    int* dsListb[6] = {&dsAb, &dsBb, &dsCb, &dsDb, &dsEb, &dsFb};
    return get_spinor_spins(dsList, dsLista, dsListb, 6, iterator);
  }
  void process::get_spinor_spins(const int* dsList[], int* dsLista[], int* dsListb[], int length, const int& iterator){
  // First map iterator into an array of local loop counters
  std::vector<int> loop_counters(length, 0);
  int temp_iterator = iterator;
  int spin0_count = 0; // add a counter for spin 0 particles
  for(int i = 0; i < length; i++) {
    if(*dsList[i] == 0){ // only particles with spin 0 have a local loop
      loop_counters[spin0_count] = (iterator >> spin0_count) & 1; // bitwise shift and bitwise AND to enumerate combinations
      spin0_count++; // increment the counter for each spin 0 particle
    }
  }

  // next, use the local loop counters to assign the spin states
  spin0_count = 0; // reset the counter for spin 0 particles
  for(int i=0;i<length;i++){
    if(*dsList[i] == 0){ //only particles with spin 0 have a local loop
      if(loop_counters[spin0_count] == 0){
        *dsLista[i] = 1;
        *dsListb[i] = -1;
      }
      else {
        *dsLista[i] = -1;
        *dsListb[i] = 1;
      }
      spin0_count++; // increment the counter for each spin 0 particle
    }
    else if(*dsList[i] < 0){
      *dsLista[i] = -1;
      *dsListb[i] = -1;
    }
    else if(*dsList[i] > 0){
      *dsLista[i] = 1;
      *dsListb[i] = 1;
    }
  }
}

    

  //Parameter Scan
  ldouble process::score_2to2_amp2(amp2ScanFunc amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20], const ldouble params[10]){
    ldouble score=0;
    ldouble En1, En2, En3, En4, Pout;
    set_test_initializations(m1, m2, m3, m4, Pin, En1, En2, Pout, En3, En4);
    ldouble p1[4], p2[4], p3[4], p4[4];
    ldouble cost;
    ldouble ampSquared;
    for(int j=0;j<20;j++){
      cost = 0.95-0.1*j;
      //Lab Frame
      set_test_momenta(p1, p2, p3, p4, En1, En2, Pin, En3, En4, Pout, cost);
      //Calculate amp^2
      set_momenta(p1,p2,p3,p4);
      ampSquared = amp2_func(params);
      if(std::isinf(ampSquared) || std::isnan(ampSquared) ||
	 (std::abs(ampSquared-data[j])>1e-20 &&
	  std::abs(ampSquared-data[j])/std::abs(ampSquared+data[j])>1e-9)){
	score += std::abs(ampSquared-data[j])/std::abs(ampSquared+data[j]);
      }
    }
    return score;
  }


  
  //  Tests
  int process::test_2to2_amp2(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]){
    int i=0;
    ldouble En1, En2, En3, En4, Pout;
    set_test_initializations(m1, m2, m3, m4, Pin, En1, En2, Pout, En3, En4);
    ldouble p1[4], p2[4], p3[4], p4[4];
    ldouble cost;
    ldouble ampSquared;
    for(int j=0;j<20;j++){
      cost = 0.95-0.1*j;
      //Lab Frame
      set_test_momenta(p1, p2, p3, p4, En1, En2, Pin, En3, En4, Pout, cost);
      //Calculate amp^2
      set_momenta(p1,p2,p3,p4);
      ampSquared = amp2_func();
      if(std::isinf(ampSquared) || std::isnan(ampSquared) ||
	 (std::abs(ampSquared-data[j])>1e-15 &&
	  std::abs(ampSquared-data[j])/std::abs(ampSquared+data[j])>1e-9)){
	  print_test_message("Lab Frame", m1, m2, m3, m4, ampSquared, data[j], cost);
	i++;
      }
    }
    return i;
  }



  
  int process::test_2to2_amp2_rotations(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]){
    int i=0;
    ldouble En1, En2, En3, En4, Pout;
    set_test_initializations(m1, m2, m3, m4, Pin, En1, En2, Pout, En3, En4);
    ldouble p1[4], p2[4], p3[4], p4[4];
    ldouble cost;
    ldouble ampSquared;
    for(int j=0;j<20;j++){
      cost = 0.95-0.1*j;
      //Lab Frame
      set_test_momenta(p1, p2, p3, p4, En1, En2, Pin, En3, En4, Pout, cost);
      //Random Rotations
      for(int l=0;l<10;l++){
	rotate_random_test_momenta(p1, p2, p3, p4);
	//Calculate amp^2
	set_momenta(p1,p2,p3,p4);
	ampSquared = amp2_func();
	if(std::isinf(ampSquared) || std::isnan(ampSquared) ||
	   (std::abs(ampSquared-data[j])>1e-15 &&
	    std::abs(ampSquared-data[j])/std::abs(ampSquared+data[j])>1e-9)){
	  print_test_message("Random Rotation", m1, m2, m3, m4, ampSquared, data[j], cost);
	  i++;
	}
      }
    }
    return i;
  }


  
  int process::test_2to2_amp2_boosts(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]){
    int i=0;
    ldouble En1, En2, En3, En4, Pout;
    set_test_initializations(m1, m2, m3, m4, Pin, En1, En2, Pout, En3, En4);
    ldouble p1[4], p2[4], p3[4], p4[4];
    ldouble cost;
    ldouble ampSquared;
    for(int j=0;j<20;j++){
      cost = 0.95-0.1*j;
      for(int l=0;l<10;l++){
	//Lab Frame
	set_test_momenta(p1, p2, p3, p4, En1, En2, Pin, En3, En4, Pout, cost);
	//Random Boost
	boost_random_test_momenta(p1, p2, p3, p4);
	//Calculate amp^2
	set_momenta(p1,p2,p3,p4);
	ampSquared = amp2_func();
	if(std::isinf(ampSquared) || std::isnan(ampSquared) ||
	   (std::abs(ampSquared-data[j])>1e-15 &&
	    std::abs(ampSquared-data[j])/std::abs(ampSquared+data[j])>1e-8)){
	  print_test_message("Random Boost", m1, m2, m3, m4, ampSquared, data[j], cost);
	  i++;
	}
      }
    }
    return i;
  }


  
  int process::test_2to2_amp2_boosts_and_rotations(amp2Func amp2_func, const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble& Pin, const ldouble data[20]){
    int i=0;
    ldouble En1, En2, En3, En4, Pout;
    set_test_initializations(m1, m2, m3, m4, Pin, En1, En2, Pout, En3, En4);
    ldouble p1[4], p2[4], p3[4], p4[4];
    ldouble cost;
    ldouble ampSquared;
    for(int j=0;j<20;j++){
      cost = 0.95-0.1*j;
      for(int l=0;l<10;l++){
	//Lab Frame
	set_test_momenta(p1, p2, p3, p4, En1, En2, Pin, En3, En4, Pout, cost);
	//Random Rotation
	rotate_random_test_momenta(p1, p2, p3, p4);
	//Random Boost
	boost_random_test_momenta(p1, p2, p3, p4);
	//Random Rotation
	rotate_random_test_momenta(p1, p2, p3, p4);
	//Calculate amp^2
	set_momenta(p1,p2,p3,p4);
	ampSquared = amp2_func();
	if(std::isinf(ampSquared) || std::isnan(ampSquared) ||
	   (std::abs(ampSquared-data[j])>1e-15 &&
	    std::abs(ampSquared-data[j])/std::abs(ampSquared+data[j])>1e-8)){
	  print_test_message("Random Rotation-Boost-Rotation", m1, m2, m3, m4, ampSquared, data[j], cost);
	  i++;
	}
      }
    }
    return i;
  }

  void process::set_test_initializations(const ldouble& m1, const ldouble& m2, const ldouble& m3, const ldouble& m4, const ldouble &Pin,
					 ldouble &En1, ldouble &En2, ldouble &Pout, ldouble &En3, ldouble &En4) const {
    En1 = std::sqrt(Pin*Pin+m1*m1);
    En2 = std::sqrt(Pin*Pin+m2*m2);
    ldouble sqrtS = En1+En2;
    ldouble lambda12=2*sqrtS*Pin;
    ldouble ms = m3 + m4, md = m3 - m4;
    ldouble lambda34 = std::sqrt((sqrtS*sqrtS - ms*ms)*(sqrtS*sqrtS - md*md));
    Pout = lambda34/(2*sqrtS);
    En3 = std::sqrt(Pout*Pout+m3*m3);
    En4 = std::sqrt(Pout*Pout+m4*m4);
  }

  void process::set_test_momenta(ldouble p1[4], ldouble p2[4], ldouble p3[4], ldouble p4[4], ldouble En1, ldouble En2,
			const ldouble &Pin, const ldouble &En3, const ldouble &En4, const ldouble &Pout, const ldouble &cost) const {
    ldouble sint = std::sqrt(1.-cost*cost);
    p1[0] = En1;
    p1[1] = 0;
    p1[2] = 0;
    p1[3] = Pin;
    p2[0] = En2;
    p2[1] = 0;
    p2[2] = 0;
    p2[3] = -Pin;
    p3[0] = En3;
    p3[1] = Pout*sint;
    p3[2] = 0;
    p3[3] = Pout*cost;
    p4[0] = En4;
    p4[1] = -p3[1];
    p4[2] = 0;
    p4[3] = -p3[3];
  }

  void process::rotate_random_test_momenta(ldouble p1[4], ldouble p2[4], ldouble p3[4], ldouble p4[4]) const {
    ldouble u[3], um, angle;
    for(int k=0;k<3;k++) u[k] = choose_random_ldouble(0,1000);
    um = std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    for(int k=0;k<3;k++) u[k] /= um;
    angle = choose_random_ldouble(0,2.*3.141592653589793238462643383279502884);
    rotate_momentum(p1,u,angle);
    rotate_momentum(p2,u,angle);
    rotate_momentum(p3,u,angle);
    rotate_momentum(p4,u,angle);
  }

  void process::boost_random_test_momenta(ldouble p1[4], ldouble p2[4], ldouble p3[4], ldouble p4[4]) const {
    ldouble v[3], vm;    
    vm = 2;
    while(vm>=0.99){
      for(int k=0;k<3;k++) v[k] = choose_random_ldouble(0,1);
      vm = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }
    boost_momentum(p1,v);
    boost_momentum(p2,v);
    boost_momentum(p3,v);
    boost_momentum(p4,v);
  }

  void process::print_test_message(const char *frame_string, const ldouble &m1, const ldouble &m2, const ldouble &m3, const ldouble &m4,
				   const ldouble &ampSquared, const ldouble &amp2_data, const ldouble &cost) const {
    std::cout<<"-----"<<frame_string<<"-------"<<std::endl;
    std::cout<<"m1,m2,m3,m4="<<m1<<","<<m2<<","<<m3<<","<<m4<<std::endl;
    std::cout<<"amp2 = "<<ampSquared<<" amp2_data="<<amp2_data<<" @ cos(theta) = "<<cost<<std::endl;
    std::cout<<"    amp2/amp2_data = ";
    printf("%.10E",double(ampSquared/amp2_data));
    std::cout<<std::endl;
  }

 
}

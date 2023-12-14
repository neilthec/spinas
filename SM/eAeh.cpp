
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

//File:  SPINAS/SM/eAeh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eAeh.h"

namespace spinas {

  eAeh::eAeh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& vev):
    e(echarge), me(masse), mh(massh), v(vev), prop(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(me);
    p2=particle(0);
    p3=particle(me);
    p4=particle(mh);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    a243s = sproduct(ANGLE,&p2,&p4,&p3);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    a241s = sproduct(ANGLE,&p2,&p4,&p1);
    s2342s = sproduct(SQUARE,&p2,&p3,&p4,&p2);
    a2342a = sproduct(ANGLE,&p2,&p3,&p4,&p2);
  }
  void eAeh::set_masses(const ldouble& masse, const ldouble& massh){
    me=masse;
    mh=massh;
    p1.set_mass(me);
    p3.set_mass(me);
    p4.set_mass(mh);
    prop.set_mass(me);
  }
  void eAeh::set_v(const ldouble& vev){
    v = vev;
  }
  void eAeh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s243a.update();
    a243s.update();
    s241a.update();
    a241s.update();
    s2342s.update();
    a2342a.update();
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=prop.den(propSP);
    pDenU=prop.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eAeh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*mh^2[12][23] - me[12][243> + me[23][241> - <13>[2342]
      mh^2<12><23> - me<12><243] + me<23><241] - [13]<2342>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    if(ds2>0){
      //mh^2[12][23] - me[12][243> - me[23][241> + <13>[2342]
      return sqrt2*e*me/v*(mh*mh*s12s.v(ds1)*s23s.v(ds3) - me*s12s.v(ds1)*s243a.v(ds3) - me*s23s.v(ds3)*s241a.v(ds1) + a13a.v(ds1,ds3)*s2342s.v())/pDenS/pDenU;
    }
    else if(ds2<0){
      //-mh^2<12><23> + me<12><243] + me<23><241] - [13]<2342>
      return sqrt2*e*me/v*(- mh*mh*a12a.v(ds1)*a23a.v(ds3) + me*a12a.v(ds1)*a243s.v(ds3) + me*a23a.v(ds3)*a241s.v(ds1) - s13s.v(ds1,ds3)*a2342a.v())/pDenS/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eAeh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2^2=1/4
    return amp2/4.0;
  }
  
  //set_momenta(...) must be called before amp2_Aplus().
  //Positive Helicity Photon Only
  ldouble eAeh::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	M = amp(j1,2,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/2^1
    return amp2/2.0;
  }
  
  //Alternate version that is a bit more efficient but requires more care from the author.
  //No need to call set_momenta(...) before this.
  ldouble eAeh::amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    ldouble propDS=std::real(prop.den(propSP));
    ldouble propDU=std::real(prop.den(propUP));

    ldouble amp2 = 0;

    //First calculate the spinor products
    cdouble s12[2], a12[2], s23[2], a23[2], s243[2], a243[2], s241[2], a241[2], a13[2][2], s13[2][2], s2342, a2342, numP, numM;
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	a13[i][j] = p1.langle(2*i-1)*p3.rangle(2*j-1);
	s13[i][j] = p1.lsquare(2*i-1)*p3.rsquare(2*j-1);
      }
      s12[i] = p1.lsquare(2*i-1)*p2.rsquare();
      a12[i] = p1.langle(2*i-1)*p2.rangle();
      s23[i] = p2.lsquare()*p3.rsquare(2*i-1);
      a23[i] = p2.langle()*p3.rangle(2*i-1);
      s243[i] = p2.lsquare()*p4.umat()*p3.rangle(2*i-1);
      a243[i] = p2.langle()*p4.lmat()*p3.rsquare(2*i-1);
    }
    s2342 = p2.lsquare()*p3.umat()*p4.lmat()*p2.rsquare();
    a2342 = p2.langle()*p3.lmat()*p4.umat()*p2.rangle();

    //Now use these to calculate the amp.
    //Sum over spins
    for(int j1=0;j1<2;j1++)
      for(int j3=0;j3<2;j3++){
	/*mh^2[12][23] - me[12][243> + me[23][241> - <13>[2342]
	  mh^2<12><23> - me<12><243] + me<23><241] - [13]<2342>*/
	//Sign changes due to p3 and p4 being outgoing
	/* mh^2[12][23] - me[12][243> - me[23][241> + <13>[2342]
	  -mh^2<12><23> + me<12><243] + me<23><241] - [13]<2342>*/
	numP =   mh*mh*s12[j1]*s23[j3] - me*s12[j1]*s243[j3] - me*s23[j3]*s241[j1] + a13[j1][j3]*s2342;
	numM = - mh*mh*a12[j1]*a23[j3] + me*a12[j1]*a243[j3] + me*a23[j3]*a241[j1] - s13[j1][j3]*a2342;

	amp2 += std::pow(std::abs(numP),2) +
	  std::pow(std::abs(numM),2);
      }
    //Average over initial spins 1/2^2=1/4
    return 2.0*e*e*me*me*amp2/v/v/propDS/propDS/propDU/propDU/4.0;
  }


  



  //  Tests
  int test_eAeh(){
    int n=0;//Number of fails
    std::cout<<"\t* e , A  -> e , h       :";
    {//amp^2
      int i=0;
      // me=0.0005, mmu=0.105, pspatial=250
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble vev=2.*MW*SW/EE;
      eAeh eAehAmp = eAeh(EE,me,mh,vev);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.022591814703232E-15,4.234637805402466E-15,8.922360995604105E-15,1.653588451249027E-14,2.764149229025680E-14,4.296168452080931E-14,6.343303552698764E-14,9.029182341493998E-14,1.252043360668657E-13,1.704716377010380E-13,2.293636616802653E-13,3.066891335380851E-13,4.098214324259149E-13,5.506714015883387E-13,7.498123157205161E-13,1.046104804802111E-12,1.522821647360727E-12,2.396472067567082E-12,4.460923220706207E-12,1.486097950409267E-11};
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH);
      //Positive Helicity Photon Only  (CH produces the same M^2s since it is averaging over spins and the to helicity cases give the same result.)    
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH);
      //Close to Higgs threshold
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.553213675482786E-14,3.913586224264316E-14,4.493068691511204E-14,5.331499244216995E-14,6.478996866663906E-14,7.999506469832603E-14,9.975919604994217E-14,1.251765311982659E-13,1.577218189598221E-13,1.994316172240947E-13,2.531999819392926E-13,3.232828777351256E-13,4.162061170437989E-13,5.425113968711459E-13,7.204071029544367E-13,9.842868527971993E-13,1.407851255334925E-12,2.182721893797211E-12,4.011511240006536E-12,1.321773153972579E-11};
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH2);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH2);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH2);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH2);
      //Positive Helicity Photon Only
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH2);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH2);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH2);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH2);
      //me~mh and close to threshold
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eAehAmp.set_masses(me,mh);
      ldouble dataCH4[20] = {9.081920004320313E-03,9.601790086425491E-03,1.012320513038232E-02,1.064572472843538E-02,1.116886444288800E-02,1.169209179548107E-02,1.221482186619280E-02,1.273641245950888E-02,1.325615879121750E-02,1.377328764313031E-02,1.428695092672449E-02,1.479621858943007E-02,1.530007078902417E-02,1.579738925218839E-02,1.628694772256517E-02,1.676740139140841E-02,1.723727518992671E-02,1.769495080638546E-02,1.813865227263857E-02,1.856642994361950E-02};
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH4);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH4);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH4);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH4);
      //Positive Helicity Photon Only
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH4);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH4);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH4);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH4);
      //Invert me and mh and close to threshold
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eAehAmp.set_masses(me,mh);
      ldouble dataCH3[20] = {9.150876604500420E-04,2.911208440362996E-03,5.149879565003035E-03,7.660069239357572E-03,1.047529955053485E-02,1.363448745031877E-02,1.718293495424712E-02,2.117347028204207E-02,2.566770930915505E-02,3.073733261475564E-02,3.646510981067254E-02,4.294504418917361E-02,5.028021613800197E-02,5.857510030083107E-02,6.791488669145936E-02,7.831384523625640E-02,8.958715827834346E-02,1.010218495012650E-01,1.104729734183341E-01,1.116053988448723E-01};
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH3);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH3);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH3);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2(); }, me,0,me,mh,pspatial,dataCH3);
      //Positive Helicity Photon Only
      i += eAehAmp.test_2to2_amp2([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH3);
      i += eAehAmp.test_2to2_amp2_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH3);
      i += eAehAmp.test_2to2_amp2_boosts([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH3);
      i += eAehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAehAmp.amp2_Aplus(); }, me,0,me,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

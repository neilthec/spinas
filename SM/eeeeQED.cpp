
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

//File:  SPINAS/SM/eeeeQED.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeeeQED.h"

namespace spinas {

  eeeeQED::eeeeQED(const ldouble& echarge, const ldouble& masse):
    e(echarge), me(masse), prop(0,0){
    p1=particle(me);
    p2=particle(me);
    p3=particle(me);
    p4=particle(me);
    a13a = sproduct(ANGLE,&p1,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a14a = sproduct(ANGLE,&p1,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a23a = sproduct(ANGLE,&p2,&p3);
    s23s = sproduct(SQUARE,&p2,&p3);
    a24a = sproduct(ANGLE,&p2,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a12a = sproduct(ANGLE,&p1,&p2);
    s34s = sproduct(SQUARE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a34a = sproduct(ANGLE,&p3,&p4);
  }
  void eeeeQED::set_masses(const ldouble& masse){
    me=masse;
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(me);
    p4.set_mass(me);
  }
  void eeeeQED::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a13a.update();
    s13s.update();
    a14a.update();
    s14s.update();
    a23a.update();
    s23s.update();
    a24a.update();
    s24s.update();
    a12a.update();
    s34s.update();
    s12s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
    }
    pDenS=prop.den(propSP);
    pDenT=prop.den(propTP);
  }

  

  
  //Amplitude
  //set_momenta(...) must be called before amp(....).
  cdouble eeeeQED::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //Sign changes due to p3 and p4 being outgoing.
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)/pDenS
    //- (<12>[34] + <14>[23] + [12]<34> + [14]<23>)/pDenT
    return -2.0*e*e*(a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3))/pDenS
      -     2.0*e*e*(a12a.v(ds1,ds2)*s34s.v(ds3,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3))/pDenT;
  }
  //set_momenta(...) must be called before amp2().
  ldouble eeeeQED::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    return amp2/4.0;
  }
  //Alternate version that is a bit more efficient but requires more care from the author.
  //No need to call set_momenta(...) before this.
  ldouble eeeeQED::amp2(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Propagator Momentum
    ldouble propSP[4], propTP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
    }
    pDenS=prop.den(propSP);
    pDenT=prop.den(propTP);
    

    ldouble amp2 = 0;

    cdouble a13[2][2], s13[2][2], a14[2][2], s14[2][2], a23[2][2], s23[2][2], a24[2][2], s24[2][2], a12[2][2], s34[2][2], s12[2][2], a34[2][2];
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	a13[i][j] = p1.langle(2*i-1)*p3.rangle(2*j-1);
	s13[i][j] = p1.lsquare(2*i-1)*p3.rsquare(2*j-1);
	a14[i][j] = p1.langle(2*i-1)*p4.rangle(2*j-1);
	s14[i][j] = p1.lsquare(2*i-1)*p4.rsquare(2*j-1);
	a23[i][j] = p2.langle(2*i-1)*p3.rangle(2*j-1);
	s23[i][j] = p2.lsquare(2*i-1)*p3.rsquare(2*j-1);
	a24[i][j] = p2.langle(2*i-1)*p4.rangle(2*j-1);
	s24[i][j] = p2.lsquare(2*i-1)*p4.rsquare(2*j-1);
	a12[i][j] = p1.langle(2*i-1)*p2.rangle(2*j-1);
	s34[i][j] = p3.lsquare(2*i-1)*p4.rsquare(2*j-1);
	s12[i][j] = p1.lsquare(2*i-1)*p2.rsquare(2*j-1);
	a34[i][j] = p3.langle(2*i-1)*p4.rangle(2*j-1);
      }
    //Calculate the amp^2
    cdouble M;
    ldouble two = 2;
    //Sum over spins
    for(int j1=0;j1<2;j1++)
      for(int j2=0;j2<2;j2++)
	for(int j3=0;j3<2;j3++)
	  for(int j4=0;j4<2;j4++){
	    //Sign changes due to p3 and p4 being outgoing.
	    //- e^2*(<13>[24] + <14>[23] + [13]<24> + [14]<23>)/pDenS
	    //- e^2*(<12>[34] + <14>[23] + [12]<34> + [14]<23>)/pDenT
	    M = - two*(a13[j1][j3]*s24[j2][j4] + a14[j1][j4]*s23[j2][j3] + s13[j1][j3]*a24[j2][j4] + s14[j1][j4]*a23[j2][j3])/pDenS
	      -   two*(a12[j1][j2]*s34[j3][j4] + a14[j1][j4]*s23[j2][j3] + s12[j1][j2]*a34[j3][j4] + s14[j1][j4]*a23[j2][j3])/pDenT;

	    amp2 += std::pow(std::abs(M),2);
	  }

    //Average over initial spins 1/2*1/2=1/4
    return e*e*e*e*amp2/4.0;
  }


  



  //  Tests
  int test_eeeeQED(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> e , E  (QED):";
    {//amp^2
      int i=0;
      // me=0.0005, pspatial=250
      ldouble me=0.0005;
      eeeeQED eeee = eeeeQED(0.31333,me);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.871563061112942E+01,5.936012537439090E+00,1.957210979996527E+00,9.216345343259433E-01,5.191209307696146E-01,3.267840212188762E-01,2.224262785273987E-01,1.607080319744738E-01,1.218716475954900E-01,9.627792825246460E-02,7.881254444247270E-02,6.658016286402389E-02,5.785489151073496E-02,5.156384234713712E-02,4.701648818953046E-02,4.375525014799866E-02,4.146932248102294E-02,3.994308122428439E-02,3.902418760902603E-02,3.860330743648967E-02};
      i += eeee.test_2to2_amp2([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeee.test_2to2_amp2_rotations([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeee.test_2to2_amp2_boosts([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeee.test_2to2_amp2_boosts_and_rotations([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH);
      // me=0.0005, pspatial=0.0001
      pspatial = 0.0001;
      ldouble dataCH2[20] = {1.117150686913463E+04,1.225628210347236E+03,4.356317898014751E+02,2.194270939672959E+02,1.310376625237332E+02,8.658811679350988E+01,6.119088892168838E+01,4.536152527090215E+01,3.485262717835317E+01,2.753312103931522E+01,2.223919225224745E+01,1.829203933406196E+01,1.527435211414919E+01,1.291827968521231E+01,1.104560593976382E+01,9.534083797699072E+00,8.297634573673667E+00,7.274254310141318E+00,6.418382287808866E+00,5.695946032560230E+00};
      i += eeee.test_2to2_amp2([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeee.test_2to2_amp2_rotations([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeee.test_2to2_amp2_boosts([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeee.test_2to2_amp2_boosts_and_rotations([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH2);
      // me=0.10, pspatial=0.05
      me=0.10;
      pspatial = 0.05;
      eeee.set_masses(me);
      ldouble dataCH3[20] = {5.411094133789898E+02,5.707805268575833E+01,1.949015874365259E+01,9.423541497121073E+00,5.397421917447446E+00,3.417828689338026E+00,2.312677921429887E+00,1.640169708614789E+00,1.204621473587737E+00,9.089382006405856E-01,7.006885842916543E-01,5.496391857261927E-01,4.374176782978349E-01,3.523700447544143E-01,2.868346948015292E-01,2.356199523091129E-01,1.951153336937824E-01,1.627528729463165E-01,1.366698141322213E-01,1.154913054245666E-01};
      i += eeee.test_2to2_amp2([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH3);
      i += eeee.test_2to2_amp2_rotations([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH3);
      i += eeee.test_2to2_amp2_boosts([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH3);
      i += eeee.test_2to2_amp2_boosts_and_rotations([&]() { return eeee.amp2(); }, me,me,me,me,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 


   
  

}

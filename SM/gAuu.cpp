
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

//File:  SPINAS/SM/gAuu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gAuu.h"

namespace spinas {

  gAuu::gAuu(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu):
    e(echarge), Qu(2.0/3.0), gs(gscharge), mu(massu), prop(massu,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(mu);
    p4=particle(mu);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
  }
  void gAuu::set_masses(const ldouble& massu){
    mu=massu;
    p3.set_mass(mu);
    p4.set_mass(mu);
    prop.set_mass(mu);
  }
  void gAuu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s132a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT = prop.den(propTP);
    pDenU = prop.den(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gAuu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds2>0){
      //gAuu:   mu[12]^2<34>
      //34 out: mu[12]^2<34>
      return 2.0*e*Qu*gs*mu*s12s.v()*s12s.v()*a34a.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1<0&&ds2<0){
      //gAuu:   mu<12>^2[34]
      //34 out: mu<12>^2[34]
      return 2.0*e*Qu*gs*mu*a12a.v()*a12a.v()*s34s.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1>0&&ds2<0){
      //gAuu:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*e*Qu*gs*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v()/pDenT/pDenU;
    }
    else if(ds1<0&&ds2>0){
      //gAuu:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*e*Qu*gs*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble gAuu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 //The color factor is the same for both diagrams
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8
    return amp2/32.0;
  }

  //A+, A+ -> e, E
  ldouble gAuu::amp2_gplus_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,2,j3,j4);
	amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
      }
    //Average over initial colors 1/8
    return amp2/8.0;
  }

  //A+, A- -> u, U
  ldouble gAuu::amp2_gplus_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,-2,j3,j4);
	amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
      }
    //Average over initial colors 1/8
    return amp2/8.0;
  }



  //  Tests
  int test_gAuu(){
    int n=0;//Number of fails
    std::cout<<"\t* g , A  -> u , U       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  mu=0.0042, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333, gs=1.238;
      gAuu gAuuAmp = gAuu(EE,gs,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.609833328985836E+00,8.302119250171076E-01,4.776773984854614E-01,3.294526542972494E-01,2.497619315962512E-01,2.016727023736482E-01,1.710928848141727E-01,1.515829612101962E-01,1.399069455491202E-01,1.344200960842295E-01,1.344200960842295E-01,1.399069455491201E-01,1.515829612101962E-01,1.710928848141726E-01,2.016727023736483E-01,2.497619315962512E-01,3.294526542972494E-01,4.776773984854615E-01,8.302119250171084E-01,2.609833328985840E+00};
      i += gAuuAmp.test_2to2_amp2([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gAuuAmp.test_2to2_amp2_rotations([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gAuuAmp.test_2to2_amp2_boosts([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      i += gAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCHpp[20] = {1.588407925571510E-08,1.960854204192353E-09,7.888876658365939E-10,4.527592231183086E-10,3.103722256383079E-10,2.374159595116177E-10,1.960997535619606E-10,1.718021957810840E-10,1.580293689147378E-10,1.517558715342471E-10,1.517558694065003E-10,1.580293393837811E-10,1.718021970008115E-10,1.960997441260136E-10,2.374159406532761E-10,3.103722217826139E-10,4.527592136993022E-10,7.888876597769202E-10,1.960854205945711E-09,1.588407928852280E-08};
      i += gAuuAmp.test_2to2_amp2([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += gAuuAmp.test_2to2_amp2_rotations([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += gAuuAmp.test_2to2_amp2_boosts([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      i += gAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCHpp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCHpm[20] = {5.219666642087591E+00,1.660423848073361E+00,9.553547961820352E-01,6.589053081417395E-01,4.995238628821302E-01,4.033454045098805E-01,3.421857694322456E-01,3.031659222485902E-01,2.798138909402110E-01,2.688401920167031E-01,2.688401920167031E-01,2.798138909402110E-01,3.031659222485902E-01,3.421857694322455E-01,4.033454045098805E-01,4.995238628821302E-01,6.589053081417395E-01,9.553547961820353E-01,1.660423848073362E+00,5.219666642087602E+00};
      i += gAuuAmp.test_2to2_amp2([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += gAuuAmp.test_2to2_amp2_rotations([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += gAuuAmp.test_2to2_amp2_boosts([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      i += gAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  mu=0.0042, pspatial=0.006\n";
      pspatial = 0.006;
      ldouble dataCH2[20] = {3.842461341366394E-01,3.363397625956756E-01,2.988679463586453E-01,2.699432241644780E-01,2.477353001109418E-01,2.308554314363499E-01,2.183129220063971E-01,2.094217191870949E-01,2.037247207071231E-01,2.009424965177974E-01,2.009424965177974E-01,2.037247207071231E-01,2.094217191870948E-01,2.183129220063971E-01,2.308554314363499E-01,2.477353001109417E-01,2.699432241644778E-01,2.988679463586451E-01,3.363397625956753E-01,3.842461341366390E-01};
      i += gAuuAmp.test_2to2_amp2([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += gAuuAmp.test_2to2_amp2_rotations([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += gAuuAmp.test_2to2_amp2_boosts([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      i += gAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAuuAmp.amp2(); }, 0,0,mu,mu,pspatial,dataCH2);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCH2pp[20] = {6.794393213797689E-01,4.962661488649108E-01,3.891924652770362E-01,3.215749214342284E-01,2.767179998449365E-01,2.461371637330606E-01,2.251801128504456E-01,2.111703195453086E-01,2.025444898891490E-01,1.984284339548514E-01,1.984284339548514E-01,2.025444898891489E-01,2.111703195453085E-01,2.251801128504456E-01,2.461371637330606E-01,2.767179998449364E-01,3.215749214342283E-01,3.891924652770360E-01,4.962661488649107E-01,6.794393213797684E-01};
      i += gAuuAmp.test_2to2_amp2([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += gAuuAmp.test_2to2_amp2_rotations([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += gAuuAmp.test_2to2_amp2_boosts([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      i += gAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAuuAmp.amp2_gplus_Aplus(); }, 0,0,mu,mu,pspatial,dataCH2pp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCH2pm[20] = {8.905294689350994E-02,1.764133763264402E-01,2.085434274402544E-01,2.183115268947275E-01,2.187526003769471E-01,2.155736991396393E-01,2.114457311623486E-01,2.076731188288812E-01,2.049049515250972E-01,2.034565590807434E-01,2.034565590807434E-01,2.049049515250972E-01,2.076731188288811E-01,2.114457311623486E-01,2.155736991396392E-01,2.187526003769470E-01,2.183115268947274E-01,2.085434274402543E-01,1.764133763264399E-01,8.905294689350955E-02};
      i += gAuuAmp.test_2to2_amp2([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += gAuuAmp.test_2to2_amp2_rotations([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += gAuuAmp.test_2to2_amp2_boosts([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      i += gAuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAuuAmp.amp2_gplus_Aminus(); }, 0,0,mu,mu,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

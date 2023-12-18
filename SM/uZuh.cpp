
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

//File:  SPINAS/SM/uZuh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uZuh.h"

namespace spinas {

  uZuh::uZuh(const ldouble& echarge, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mu(massu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propu(massu,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(mu);
    p2=particle(MZ);
    p3=particle(mu);
    p4=particle(mh);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[312>,[213>,[343>,[143>,[341>
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s242a = sproduct(SQUARE,&p2,&p4,&p2);
    s142a = sproduct(SQUARE,&p1,&p4,&p2);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    //Couplings
    preTU = sqrt2*e*e*mu/(4.0*MW*MW*SW*SW);
    gL=-4.0/3.0*SW*SW+1.0;
    gR=-4.0/3.0*SW*SW;
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void uZuh::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mu);
    p2.set_mass(MZ);
    p3.set_mass(mu);
    p4.set_mass(mh);
    propu.set_mass(mu);
    propZ.set_mass(MZ);
    //Couplings
    preTU = sqrt2*e*e*mu/(4.0*MW*MW*SW*SW);
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void uZuh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    //[312>,[213>,[343>,[143>,[341>
    s213a.update();
    s312a.update();
    s242a.update();
    s142a.update();
    s241a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propu.denominator(propSP);
    pDenT=propZ.denominator(propTP);
    pDenU=propu.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uZuh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;

    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);
      
      //S-Channel Z
      //preZ = e*e/(4.0*MW*MW*SW*SW);
      //uUZh all ingoing: 
      //+ preZ (Mu (gLe-gRe) (<12>-[12]) [343> + 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //uZUh: 2<->3
      //+ preZ (Mu (gLe-gRe) (<13>-[13]) [242> - 2 MZ^2(gLe[23] <12> + gRe[12] <23>)))/(t-MZ^2)
      //34 out:
      //+ preZ (Mu (gLe-gRe) (<13>+[13]) [242> - 2 MZ^2(gLe[23] <12> - gRe[12] <23>)))/(t-MZ^2)
      amplitude += normFactor*preZ*( mu*(gL-gR)*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*s242a.v(ds2a,ds2b) - 2.0*MZ*MZ*(gL*s23s.v(ds2a,ds3)*a12a.v(ds1,ds2b) - gR*a23a.v(ds2a,ds3)*s12s.v(ds1,ds2b)))/pDenT;
      
      //T-Channel e
      //preTU = e*e*mu/(4.0*MW*MW*SW*SW);
      //uUZh all ingoing:
      //+ preh*(gLe <13> (Mu [23]-[312>+MZ <23>)+gRe [13] (MZ [23]-[213>+Mu <23>)))/(t-Mu^2)
      //uZUh: 2<->3
      //+ preh*(gLe <12> (-Mu [23]-[213>-MZ <23>)+gRe [12] (-MZ [23]-[312>-Mu <23>)))/(s-Mu^2)
      //34 out:
      //+ preh*(gLe <12> (-Mu [23]+[213>+MZ <23>)+gRe [12] (-MZ [23]-[312>+Mu <23>)))/(s-Mu^2)
      amplitude += normFactor*preTU*( gL*a12a.v(ds1,ds2a)*(-mu*s23s.v(ds2b,ds3)+s213a.v(ds2b,ds3)+MZ*a23a.v(ds2b,ds3)) + gR*s12s.v(ds1,ds2a)*(mu*a23a.v(ds2b,ds3)-s312a.v(ds3,ds2b)-MZ*s23s.v(ds2b,ds3)) )/pDenS;

      //U-Channel e
      //preTU = e*e*mu/(4.0*MW*MW*SW*SW);
      //uUZh all ingoing:
      //+ preh (gLe [23] ([143>+2 Mu <13>)+gRe <23> (2 Mu [13]+[341>) )/(u-Mu^2)
      //uZUh: 2<->3
      //- preh (gLe [23] ([142>+2 Mu <12>)+gRe <23> (2 Mu [12]+[241>) )/(u-Mu^2)
      //34 out:
      //- preh (gLe [23] (-[142>+2 Mu <12>)-gRe <23> (2 Mu [12]-[241>) )/(u-Mu^2)
      amplitude += -normFactor*preTU*( gL*s23s.v(ds2a,ds3)*(-s142a.v(ds1,ds2b)+2.0*mu*a12a.v(ds1,ds2b)) - gR*a23a.v(ds2a,ds3)*(-s241a.v(ds2b,ds1)+2.0*mu*s12s.v(ds1,ds2b)) )/pDenU;
    }

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uZuh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/3
    return amp2/18.0;
  }
  



  //  Tests
  int test_uZuh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , Z  -> u , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      uZuh uZuhAmp = uZuh(EE,mu,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.338606532238737E+00,6.553997723514051E-01,2.929974084560649E-01,1.607208670088717E-01,9.884403582762943E-02,6.534370525406891E-02,4.536499125458444E-02,3.260205307008907E-02,2.401891463198660E-02,1.801242118607368E-02,1.367385516141607E-02,1.045836943308995E-02,8.024027511215488E-03,6.147957655445044E-03,4.680184535010175E-03,3.516997792183091E-03,2.584970601956794E-03,1.831037002868180E-03,1.216145610541731E-03,7.110936858635595E-04};
      i += uZuhAmp.test_2to2_amp2([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH);
      //i += uZuhAmp.test_2to2_amp2_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH);
      //i += uZuhAmp.test_2to2_amp2_boosts([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH);
      //i += uZuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.491704724962019E-01,2.020203965578147E-01,1.292846403340397E-01,8.838390713503516E-02,6.328212531857208E-02,4.686969286536376E-02,3.561163266845574E-02,2.759369926346900E-02,2.170838547325242E-02,1.728013963704075E-02,1.387867649168208E-02,1.121986671577469E-02,9.110274489061295E-03,7.414713316626544E-03,6.036534495484997E-03,4.905252858216252E-03,3.968553369898092E-03,3.186995938725167E-03,2.530428416299414E-03,1.975519986453076E-03};
      //i += uZuhAmp.test_2to2_amp2([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH2);
      //i += uZuhAmp.test_2to2_amp2_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH2);
      //i += uZuhAmp.test_2to2_amp2_boosts([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH2);
      //i += uZuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, mh=125, MW=80.385, pspatial=95\n";
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uZuhAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {3.187054511007267E-01,2.483609584329007E-01,2.018885012961885E-01,1.699670029214243E-01,1.474255634669530E-01,1.312204705642079E-01,1.194752577594514E-01,1.109911512675071E-01,1.049825151823831E-01,1.009271044255139E-01,9.847826512559818E-02,9.741233122836643E-02,9.759712151309063E-02,9.897395316268086E-02,1.015491802078270E-01,1.053934678490681E-01,1.106486387555627E-01,1.175435137064815E-01,1.264222451243715E-01,1.377919888468090E-01};
      //i += uZuhAmp.test_2to2_amp2([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH4);
      //i += uZuhAmp.test_2to2_amp2_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH4);
      //i += uZuhAmp.test_2to2_amp2_boosts([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH4);
      //i += uZuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uZuhAmp.set_masses(mu,mh,MW);
      ldouble dataCH3[20] = {6.200168402908511E-01,3.390311550331094E-01,2.170020130974172E-01,1.543763852462991E-01,1.187005827609352E-01,9.694768986396513E-02,8.311586342593842E-02,7.416075589425289E-02,6.843265216328426E-02,6.500353711954587E-02,6.335311408921382E-02,6.321603806627263E-02,6.451037951358415E-02,6.731481391788434E-02,7.188544286683653E-02,7.872091105091214E-02,8.870968325075045E-02,1.034488133064262E-01,1.259758454333256E-01,1.626494200647483E-01};
      //i += uZuhAmp.test_2to2_amp2([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH3);
      //i += uZuhAmp.test_2to2_amp2_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH3);
      //i += uZuhAmp.test_2to2_amp2_boosts([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH3);
      //i += uZuhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZuhAmp.amp2(); }, mu,MZ,mu,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}


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

//File:  SPINAS/SM/ddZh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddZh.h"

namespace spinas {

  ddZh::ddZh(const ldouble& echarge, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), md(massd), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propd(massd,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(md);
    p2=particle(md);
    p3=particle(MZ);
    p4=particle(mh);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[312>,[213>,[343>,[143>,[341>
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    s343a = sproduct(SQUARE,&p3,&p4,&p3);
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    s342a = sproduct(SQUARE,&p3,&p4,&p2);
    //Couplings
    preTU = sqrt2*e*e*md/(4.0*MW*MW*SW*SW);
    gL=2.0/3.0*SW*SW-1.0;
    gR=2.0/3.0*SW*SW;
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void ddZh::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(md);
    p2.set_mass(md);
    p3.set_mass(MZ);
    p4.set_mass(mh);
    propd.set_mass(md);
    propZ.set_mass(MZ);
    //Couplings
    preTU = sqrt2*e*e*md/(4.0*MW*MW*SW*SW);
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void ddZh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s312a.update();
    s213a.update();
    s343a.update();
    s143a.update();
    s341a.update();
    s243a.update();
    s342a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propZ.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propd.denominator(propTP);
    pDenU=propd.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddZh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;

    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //S-Channel Z
      //preZ = e*e/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preZ (-Md (gLe-gRe) (<12>-[12]) (MZ [33]-[343>) + 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //34 outgoing:
      //preZ (-Md (gLe-gRe) (<12>-[12]) (MZ [33]-[343>) - 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //[33]=0 since it is symmetrized
      //preZ ( Md (gLe-gRe) (<12>-[12]) [343> - 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      amplitude += normFactor*preZ*( md*(gL-gR)*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*s343a.v(ds3a,ds3b) - 2.0*MZ*MZ*(gL*s23s.v(ds2,ds3a)*a13a.v(ds1,ds3b) + gR*a23a.v(ds2,ds3a)*s13s.v(ds1,ds3b)))/pDenS;
      
      //T-Channel e
      //preTU = e*e*md/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //preh*(gLe <13> (Md [23]-[312>+MZ <23>)+gRe [13] (MZ [23]-[213>+Md <23>)))/(t-Md^2)
      //34 outgoing:
      //- preh*( gLe <13> (Md [23]-[312>-MZ <23>) + gRe [13] (Md <23>-[213>-MZ [23])) )/(t-Md^2)
      //amplitude += -normFactor*preTU*( gL*a13a.v(ds1,ds3a)*(md*s23s.v(ds2,ds3b)-s312a.v(ds3b,ds2)-MZ*a23a.v(ds2,ds3b)) + gR*s13s.v(ds1,ds3a)*(md*a23a.v(ds2,ds3b)-s213a.v(ds2,ds3b)-MZ*s23s.v(ds2,ds3b)) )/pDenT;
      amplitude += normFactor*preTU*(
				     gR*s13s.v(ds1,ds3a)*(-2.0*md*a23a.v(ds2,ds3b)+s243a.v(ds2,ds3b))
				     +gL*a13a.v(ds1,ds3a)*(-2.0*md*s23s.v(ds2,ds3b)+s342a.v(ds3b,ds2))
				     )/pDenT;

      //U-Channel e
      //preTU = e*e*md/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //preh (gLe [23] ([143>+2 Md <13>)+gRe <23> (2 Md [13]+[341>) )/(u-Md^2)
      //34 outgoing:
      //preh (gLe [23] ([143>-2 Md <13>)+gRe <23> ([341> - 2 Md [13]) )/(u-Md^2)
      amplitude += normFactor*preTU*(
				     gL*s23s.v(ds2,ds3a)*(s143a.v(ds1,ds3b)-2.0*md*a13a.v(ds1,ds3b))
				     + gR*a23a.v(ds2,ds3a)*(s341a.v(ds3b,ds1)-2.0*md*s13s.v(ds1,ds3b))
				     )/pDenU;
    }

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ddZh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }
  



  //  Tests
  int test_ddZh(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> Z , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, mh=125, MW=80.385, pspatial=250\n";
      ldouble md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      ddZh ddZhAmp = ddZh(EE,md,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.560074908456033E-04,6.481945945245533E-04,8.190277953187908E-04,9.685068800375551E-04,1.096631820312737E-03,1.203402608031152E-03,1.288819240018146E-03,1.352881714800917E-03,1.395590031632140E-03,1.416944190129631E-03,1.416944190129631E-03,1.395590031632140E-03,1.352881714800917E-03,1.288819240018147E-03,1.203402608031152E-03,1.096631820312737E-03,9.685068800375551E-04,8.190277953187909E-04,6.481945945245532E-04,4.560074908456034E-04};
      i += ddZhAmp.test_2to2_amp2([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH);
      i += ddZhAmp.test_2to2_amp2_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH);
      i += ddZhAmp.test_2to2_amp2_boosts([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH);
      i += ddZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.787275890981122E-03,1.860548916678241E-03,1.925680482336043E-03,1.982670598469260E-03,2.031519268106150E-03,2.072226492299318E-03,2.104792271458588E-03,2.129216605754586E-03,2.145499495260637E-03,2.153640940007921E-03,2.153640940007922E-03,2.145499495260637E-03,2.129216605754586E-03,2.104792271458588E-03,2.072226492299318E-03,2.031519268106150E-03,1.982670598469261E-03,1.925680482336043E-03,1.860548916678241E-03,1.787275890981122E-03};
      i += ddZhAmp.test_2to2_amp2([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH2);
      i += ddZhAmp.test_2to2_amp2_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH2);
      i += ddZhAmp.test_2to2_amp2_boosts([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH2);
      i += ddZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, MW=80.385, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      ddZhAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {6.033282925044231E-02,5.583624162747997E-02,5.215252345259221E-02,4.915974346744843E-02,4.675536545091499E-02,4.485830640947455E-02,4.340664431323751E-02,4.235464667811276E-02,4.167020465281561E-02,4.133293086347801E-02,4.133293086347801E-02,4.167020465281560E-02,4.235464667811276E-02,4.340664431323751E-02,4.485830640947456E-02,4.675536545091499E-02,4.915974346744844E-02,5.215252345259223E-02,5.583624162747999E-02,6.033282925044232E-02};
      i += ddZhAmp.test_2to2_amp2([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH4);
      i += ddZhAmp.test_2to2_amp2_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH4);
      i += ddZhAmp.test_2to2_amp2_boosts([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH4);
      i += ddZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      ddZhAmp.set_masses(md,mh,MW);
      ldouble dataCH3[20] = {9.679175944831529E-02,8.272984360106940E-02,7.261780105161846E-02,6.517920168678068E-02,5.964078518328260E-02,5.551781750404314E-02,5.249826084750503E-02,5.037904541053353E-02,4.903008088063862E-02,4.837375977592163E-02,4.837375977592163E-02,4.903008088063862E-02,5.037904541053353E-02,5.249826084750503E-02,5.551781750404314E-02,5.964078518328261E-02,6.517920168678068E-02,7.261780105161847E-02,8.272984360106941E-02,9.679175944831529E-02};
      i += ddZhAmp.test_2to2_amp2([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH3);
      i += ddZhAmp.test_2to2_amp2_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH3);
      i += ddZhAmp.test_2to2_amp2_boosts([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH3);
      i += ddZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZhAmp.amp2(); }, md,md,MZ,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

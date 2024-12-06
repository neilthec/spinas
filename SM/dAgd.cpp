
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

//File:  SPINAS/SM/dAgd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dAgd.h"

namespace spinas {

  dAgd::dAgd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd):
    e(echarge), Qd(-1.0/3.0), gs(gscharge), md(massd), prop(massd,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(md);
    p2=particle(0);
    p3=particle(0);
    p4=particle(md);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    s342a = sproduct(SQUARE,&p3,&p4,&p2,2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3,2);
  }
  void dAgd::set_masses(const ldouble& massd){
    md=massd;
    p1.set_mass(md);
    p4.set_mass(md);
    prop.set_mass(md);
  }
  void dAgd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s342a.update();
    s243a.update();
    //Propagator Momentum
    ldouble propTP[4], propSP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propSP[j] = mom1[j]+mom2[j];
    }
    pDenT = prop.denominator(propTP);
    pDenS = prop.denominator(propSP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dAgd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds3>0&&ds2>0){
      //gADd:             + md[12]^2<34>/tu
      //dAgD: 4->1->3->4: - md[23]^2<14>/st
      //34 out:           + md[23]^2<14>/st
      return + 2.0*e*Qd*gs*md*s23s.v()*s23s.v()*a14a.v(ds1,ds4)/pDenS/pDenT;
    }
    else if(ds3<0&&ds2<0){
      //gADd:             + md<12>^2[34]/tu
      //dAgD: 4->1->3->4: - md<23>^2[14]/st
      //34 out:           - md<23>^2[14]/st
      return - 2.0*e*Qd*gs*md*a23a.v()*a23a.v()*s14s.v(ds1,ds4)/pDenS/pDenT;
    }
    else if(ds3>0&&ds2<0){
      //gADd:             + ([31]<42>+[41]<32>)[132>/tu
      //dAgD: 4->1->3->4: - ([34]<12>+[13]<24>)[342>/st
      //34 out:           + ([34]<12>-[13]<24>)[342>/st
      return + 2.0*e*Qd*gs*(s34s.v(ds4)*a12a.v(ds1)-s13s.v(ds1)*a24a.v(ds4))*s342a.v()/pDenS/pDenT;
    }
    else if(ds3<0&&ds2>0){
      //gADd:             + (<31>[42]+<41>[32])[231>
      //dAgD: 4->1->3->4: - (<34>[12]+<13>[24])[243>
      //34 out:           - (<34>[12]-<13>[24])[243>
      return -2.0*e*Qd*gs*(a34a.v(ds4)*s12s.v(ds1)-a13a.v(ds1)*s24s.v(ds4))*s243a.v()/pDenS/pDenT;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) mdst be called before amp2().
  ldouble dAgd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 //The color factor is the same for both diagrams
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/3
    return amp2/12.0;
  }

  //d, A+ -> g, d
  ldouble dAgd::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3
    return amp2/6.0;
  }

  //d, A- -> g, d
  ldouble dAgd::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,-2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3
    return amp2/6.0;
  }



  //  Tests
  int test_dAgd(){
    int n=0;//Number of fails
    std::cout<<"\t* d , A  -> g , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  md=0.0075, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333, gs=1.238;
      dAgd dAgdAmp = dAgd(EE,gs,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.784443517922574E+00,5.977867244978341E-01,3.622386930550193E-01,2.625633429434969E-01,2.081788866996537E-01,1.743812004158773E-01,1.516686979633304E-01,1.356073058121814E-01,1.238495732690622E-01,1.150364499320808E-01,1.083266188855583E-01,1.031713951843725E-01,9.919767306577673E-02,9.614289662261635E-02,9.381693916348097E-02,9.207873894956799E-02,9.082143102722975E-02,8.996257671793180E-02,8.943756192227789E-02,8.919502671510389E-02};
      i += dAgdAmp.test_2to2_amp2([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dAgdAmp.test_2to2_amp2_rotations([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dAgdAmp.test_2to2_amp2_boosts([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dAgdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      //std::cout<<"\n#  d , A+  -> g , d\n";
      ldouble dataCHpp[20] = {1.784443517922574E+00,5.977867244978341E-01,3.622386930550193E-01,2.625633429434969E-01,2.081788866996537E-01,1.743812004158773E-01,1.516686979633304E-01,1.356073058121814E-01,1.238495732690622E-01,1.150364499320808E-01,1.083266188855583E-01,1.031713951843725E-01,9.919767306577673E-02,9.614289662261635E-02,9.381693916348097E-02,9.207873894956799E-02,9.082143102722975E-02,8.996257671793180E-02,8.943756192227789E-02,8.919502671510389E-02};
      i += dAgdAmp.test_2to2_amp2([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dAgdAmp.test_2to2_amp2_rotations([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dAgdAmp.test_2to2_amp2_boosts([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dAgdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      //std::cout<<"\n#  d , A-  -> g , d\n";
      ldouble dataCHpm[20] = {1.784443517922574E+00,5.977867244978341E-01,3.622386930550193E-01,2.625633429434969E-01,2.081788866996537E-01,1.743812004158773E-01,1.516686979633304E-01,1.356073058121814E-01,1.238495732690622E-01,1.150364499320808E-01,1.083266188855583E-01,1.031713951843725E-01,9.919767306577673E-02,9.614289662261635E-02,9.381693916348097E-02,9.207873894956799E-02,9.082143102722975E-02,8.996257671793180E-02,8.943756192227789E-02,8.919502671510389E-02};
      i += dAgdAmp.test_2to2_amp2([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dAgdAmp.test_2to2_amp2_rotations([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dAgdAmp.test_2to2_amp2_boosts([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dAgdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  md=0.0075, pspatial=0.008\n";
      pspatial = 0.008;
      ldouble dataCH2[20] = {2.375459886244632E-01,1.726659894514116E-01,1.374556648555788E-01,1.166401365590591E-01,1.036409339131127E-01,9.524744360836557E-02,8.974186950698752E-02,8.613788404960215E-02,8.383651400984633E-02,8.245704547007839E-02,8.174815475902523E-02,8.153854655706765E-02,8.170823486305655E-02,8.217116154520408E-02,8.286431270426310E-02,8.374069875682848E-02,8.476470707084442E-02,8.590895339701259E-02,8.715210428401018E-02,8.847734287831657E-02};
      i += dAgdAmp.test_2to2_amp2([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dAgdAmp.test_2to2_amp2_rotations([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dAgdAmp.test_2to2_amp2_boosts([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dAgdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAgdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      //std::cout<<"\n#  d , A+  -> g , d\n";
      ldouble dataCH2pp[20] = {2.375459886244632E-01,1.726659894514116E-01,1.374556648555788E-01,1.166401365590591E-01,1.036409339131127E-01,9.524744360836557E-02,8.974186950698752E-02,8.613788404960215E-02,8.383651400984633E-02,8.245704547007839E-02,8.174815475902523E-02,8.153854655706765E-02,8.170823486305655E-02,8.217116154520408E-02,8.286431270426310E-02,8.374069875682848E-02,8.476470707084442E-02,8.590895339701259E-02,8.715210428401018E-02,8.847734287831657E-02};
      i += dAgdAmp.test_2to2_amp2([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dAgdAmp.test_2to2_amp2_rotations([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dAgdAmp.test_2to2_amp2_boosts([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dAgdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAgdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      //std::cout<<"\n#  d , A-  -> g , d\n";
      ldouble dataCH2pm[20] = {2.375459886244632E-01,1.726659894514116E-01,1.374556648555788E-01,1.166401365590591E-01,1.036409339131127E-01,9.524744360836557E-02,8.974186950698752E-02,8.613788404960215E-02,8.383651400984633E-02,8.245704547007839E-02,8.174815475902523E-02,8.153854655706765E-02,8.170823486305655E-02,8.217116154520408E-02,8.286431270426310E-02,8.374069875682848E-02,8.476470707084442E-02,8.590895339701259E-02,8.715210428401018E-02,8.847734287831657E-02};
      i += dAgdAmp.test_2to2_amp2([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dAgdAmp.test_2to2_amp2_rotations([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dAgdAmp.test_2to2_amp2_boosts([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dAgdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAgdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

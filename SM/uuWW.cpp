
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

//File:  SPINAS/SM/uuWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuWW.h"

namespace spinas {

  uuWW::uuWW(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(2.0/3.0), mu(massu), md(massd), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propd(md,0), propA(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(mu);
    p2=particle(mu);
    p3=particle(MW);
    p4=particle(MW);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[413>,[314>,[231>,[132>
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*mu/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW+1.0;
    gR=-2.0*Qf*SW*SW;
  }
  void uuWW::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW){
    mu=massu;
    md=massd;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(MW);
    p4.set_mass(MW);
    propd.set_mass(md);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*mu/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
  }
  void uuWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[413>,[314>,[231>,[132>
    s413a.update();
    s314a.update();
    s231a.update();
    s132a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenhS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDendT=propd.denominator(propTP);
    pDendU=propd.denominator(propUP);
    pDenZS=propZ.denominator(propSP);
    pDenAS=propA.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uuWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //S-Channel h
      //preh = e*e*mu/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //34 outgoing:
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      amplitude += normFactor*preh*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenhS;
      
      //T-Channel d
      //pred = e*e/(2.0*MW*MW*SW*SW);
      //all ingoing:
      // pred [24] <13> ([314>+MW <34>))/t
      //34 outgoing:
      // pred [24] <13> (-MW <34> + [314>)/t
      amplitude += - normFactor*pred*s24s.v(ds2,ds4a)*a13a.v(ds1,ds3a)*(-MW*a34a.v(ds3b,ds4b)+s314a.v(ds3b,ds4b))/pDendT;

      //S-Channel A
      //Guess based on ee->mm (and agreed with Z Diagram):
      //2e^2/MW*( <13>[24] + [13]<24> + <14>[23] + [14]<23> )( <34> + [34] )/s
      //Added based on Z Diagram
      //2e^2/MW/MW*[34]<34>([231>+[132>)/s
      amplitude +=
	normFactor*e*e*Qf*two/MW*(
				a13a.v(ds1,ds3a)*s24s.v(ds2,ds4a) + s13s.v(ds1,ds3a)*a24a.v(ds2,ds4a)
				+ a14a.v(ds1,ds4a)*s23s.v(ds2,ds3a) + s14s.v(ds1,ds4a)*a23a.v(ds2,ds3a)
				)*(s34s.v(ds3b,ds4b)+a34a.v(ds3b,ds4b))/pDenAS
	+ normFactor*e*e*Qf*two/MW/MW*(
				s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s231a.v(ds2,ds1)+s132a.v(ds1,ds2))
				)/pDenAS;

      //S-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=pred
      //all ingoing
      //+ preZ gLe ( Mu[34]<34>( [12]-<12>) + 2[34]<34>[231> + 2MW([23]<14>+[24]<13>)([34]+<34>) )/(s-MZ^2)
      //+ preZ gRe ( Mu[34]<34>(-[12]+<12>) + 2[34]<34>[132> + 2MW([13]<24>+[14]<23>)([34]+<34>) )/(s-MZ^2)
      // = preZ ( (gL-gR)Mu[34]<34>([12]-<12>)
      //          + 2[34]<34>(gL[231>+gR[132>)
      //          + 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //        )/(s-MZ^2)
      //34 outgoing
      //preZ ( (gL-gR)Mu[34]<34>([12]-<12>)
      //        - 2[34]<34>(gL[231>+gR[132>)
      //        - 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //      )/(s-MZ^2)
      amplitude += - normFactor*preZ*( (gL-gR)*mu*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)-a12a.v(ds1,ds2))
				       -two*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(gL*s231a.v(ds2,ds1)+gR*s132a.v(ds1,ds2))
				       -two*MW*(s34s.v(ds3a,ds4a)+a34a.v(ds3a,ds4a))*( gL*(s23s.v(ds2,ds3b)*a14a.v(ds1,ds4b)+s24s.v(ds2,ds4b)*a13a.v(ds1,ds3b))
										       + gR*(s13s.v(ds1,ds3b)*a24a.v(ds2,ds4b)+s14s.v(ds1,ds4b)*a23a.v(ds2,ds3b)) )
				       )/pDenZS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }

  
  

  //  Tests
  int test_uuWW(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      uuWW uuWWAmp = uuWW(EE,mu,md,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.247646142163802E-01,1.888910323223591E-01,1.000296574600400E-01,6.302685285395956E-02,4.341266153845872E-02,3.165659973774527E-02,2.406750799999728E-02,1.891540968772569E-02,1.528172335703240E-02,1.263348657269628E-02,1.064022315279356E-02,9.086122420708290E-03,7.824276865829918E-03,6.751215949257134E-03,5.791957066298060E-03,4.890832663798485E-03,4.005647228980510E-03,3.103834356577405E-03,2.159858266105236E-03,1.153413692141902E-03};
      i += uuWWAmp.test_2to2_amp2([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH);
      i += uuWWAmp.test_2to2_amp2_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH);
      i += uuWWAmp.test_2to2_amp2_boosts([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH);
      i += uuWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=2500\n";
      pspatial = 2500;
      ldouble dataCH1[20] = {5.910834351891432E-01,1.696022622016166E-01,8.773075151023706E-02,5.421579452457936E-02,3.669578091189068E-02,2.634480550550595E-02,1.976828715782862E-02,1.538208735285421E-02,1.234754429202850E-02,1.017911683728050E-02,8.576231851749161E-03,7.342863354176461E-03,6.345713660626196E-03,5.490980691114860E-03,4.710738799315221E-03,3.954591535327784E-03,3.184359763042151E-03,2.370591007300014E-03,1.490200899876400E-03,5.248403185917215E-04};
      i += uuWWAmp.test_2to2_amp2([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH1);
      i += uuWWAmp.test_2to2_amp2_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH1);
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {4.128163421739068E-01,2.064408960921143E-01,1.284017493243048E-01,8.902135251088936E-02,6.577439332779816E-02,5.070650136209449E-02,4.031554377861337E-02,3.282270294714974E-02,2.722948174280134E-02,2.293238055139696E-02,1.954525136554401E-02,1.680957055818435E-02,1.454575494036501E-02,1.262517667802449E-02,1.095329206612887E-02,9.459055727330645E-03,8.088048579622710E-03,6.797883781027288E-03,5.555055945661364E-03,4.332731036682001E-03};
      i += uuWWAmp.test_2to2_amp2([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH2);
      i += uuWWAmp.test_2to2_amp2_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH2);
      i += uuWWAmp.test_2to2_amp2_boosts([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH2);
      i += uuWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, mh=125, MW=80.385, pspatial=95\n";
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uuWWAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH4[20] = {1.009360920291309E+00,6.864734575694987E-01,4.857472334587417E-01,3.656288862148899E-01,2.872353942318390E-01,2.325185592648959E-01,1.923956855039418E-01,1.618623839841352E-01,1.379498121410491E-01,1.187898877992239E-01,1.031500761604619E-01,9.018498116811144E-02,7.929571883044502E-02,7.004631924847049E-02,6.211193405395155E-02,5.524559514545877E-02,4.925621904375153E-02,4.399366007715938E-02,3.933831309699978E-02,3.519372945454213E-02};
      i += uuWWAmp.test_2to2_amp2([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH4);
      i += uuWWAmp.test_2to2_amp2_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH4);
      i += uuWWAmp.test_2to2_amp2_boosts([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH4);
      i += uuWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uuWWAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH3[20] = {1.661158904064196E+00,8.438165492399010E-01,5.379786573029638E-01,3.835540200824558E-01,2.912174738790127E-01,2.301515056358585E-01,1.869877465538987E-01,1.550044314852050E-01,1.304560151391403E-01,1.110912095282891E-01,9.547580334526272E-02,8.265297221130229E-02,7.195998123766496E-02,6.292329229085364E-02,5.519556424784462E-02,4.851624914684415E-02,4.268608171453567E-02,3.755006544832401E-02,3.298582790892773E-02,2.889546748882251E-02};
      i += uuWWAmp.test_2to2_amp2([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH3);
      i += uuWWAmp.test_2to2_amp2_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH3);
      i += uuWWAmp.test_2to2_amp2_boosts([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH3);
      i += uuWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuWWAmp.amp2(); }, mu,mu,MW,MW,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}

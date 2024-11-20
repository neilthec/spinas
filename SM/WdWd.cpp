
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

//File:  SPINAS/SM/WdWd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/WdWd.h"

namespace spinas {

  WdWd::WdWd(const ldouble& echarge, const ldouble& massd, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(-1.0/3.0), md(massd), mu(massu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propu(mu,0), propA(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(MW);
    p2=particle(md);
    p3=particle(MW);
    p4=particle(md);
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
    //[143>,[341>,[234>,[432>
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    s234a = sproduct(SQUARE,&p2,&p3,&p4);
    s432a = sproduct(SQUARE,&p4,&p3,&p2);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*md/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW-1.0;
    gR=-2.0*Qf*SW*SW;
  }
  void WdWd::set_masses(const ldouble& massd, const ldouble& massu, const ldouble& massh, const ldouble& massW){
    md=massd;
    mu=massu;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(MW);
    p2.set_mass(md);
    p3.set_mass(MW);
    p4.set_mass(md);
    propu.set_mass(mu);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*md/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
  }
  void WdWd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[143>,[341>,[234>,[432>
    s143a.update();
    s341a.update();
    s234a.update();
    s432a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenhT=proph.denominator(propTP);
    pDendS=propu.denominator(propSP);
    pDendU=propu.denominator(propUP);
    pDenZT=propZ.denominator(propTP);
    pDenAT=propA.denominator(propTP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble WdWd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds1a, ds1b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds1);
    ldouble normFactor=get_spin_normalization(ds3,ds1);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds1,ds1a,ds1b, i);
      
      //S-Channel h
      //preh = e*e*md/(4.0*MW*MW*SW*SW);
      //DdW-W+ all ingoing: 
      //+ preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //W+dW-D: 1<->4
      //- preh [13] <13> ([24]+<24>)/(t-Mh^2)
      //34 out:
      //+ preh [13] <13> ([24]-<24>)/(t-Mh^2)
      amplitude += normFactor*preh*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenhT;
      
      //T-Channel u
      //pred = e*e/(2.0*MW*MW*SW*SW);
      //DdW-W+ all ingoing:
      //- pred [24] <13> ([314>+MW <34>))/t
      //W+dW-D: 1<->4
      //- pred [12] <34> ([341>-MW <13>))/s
      //34 out:
      //+ pred [12] <34> ([341>-MW <13>))/s
      amplitude += + normFactor*pred*s12s.v(ds1a,ds2)*a34a.v(ds3a,ds4)*(-MW*a13a.v(ds1b,ds3b)+s341a.v(ds3b,ds1b))/pDendS;

      //S-Channel A
      //DdW-W+ all in:
      //+ 2e^2/MW*( <13>[24] + [13]<24> + <14>[23] + [14]<23> )( <34> + [34] )/s
      //+ 2e^2/MW/MW*[34]<34>([231>+[132>)/s
      //W+dW-D: 1<->4
      //- 2e^2/MW*( <34>[12] + [34]<12> - <14>[23] - [14]<23> )( <13> + [13] )/t
      //+ 2e^2/MW/MW*[13]<13>([234>+[432>)/t
      //34 out:
      //- 2e^2/MW*( <34>[12] + [34]<12> + <14>[23] + [14]<23> )(-<13> + [13] )/t
      //- 2e^2/MW/MW*[13]<13>([234>-[432>)/t
      amplitude +=
	- normFactor*e*e*Qf*two/MW*(
				a34a.v(ds3a,ds4)*s12s.v(ds1a,ds2) + s34s.v(ds3a,ds4)*a12a.v(ds1a,ds2)
				+ a14a.v(ds1a,ds4)*s23s.v(ds2,ds3a) + s14s.v(ds1a,ds4)*a23a.v(ds2,ds3a)
				)*(s13s.v(ds1b,ds3b)-a13a.v(ds1b,ds3b))/pDenAT
	- normFactor*e*e*Qf*two/MW/MW*(
				s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*(s234a.v(ds2,ds4)-s432a.v(ds4,ds2))
				)/pDenAT;
      
      //S-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=pred
      //DdW-W+ all ingoing:
      //+ preZ ( (gL-gR)md[34]<34>([12]-<12>)
      //          + 2[34]<34>(gL[231>+gR[132>)
      //          + 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //        )/(s-MZ^2)
      //W+dW-D: 1<->4
      //+ preZ (-(gL-gR)md[13]<13>([24]-<24>)
      //          + 2[13]<13>(gL[234>+gR[432>)
      //          - 2MW([13]+<13>)( gL(-[23]<14>+[12]<34>) + gR([34]<12>-[14]<23>) )
      //        )/(t-MZ^2)
      //34 out:
      //+ preZ (+(gL-gR)md[13]<13>([24]+<24>)
      //          - 2[13]<13>(gL[234>-gR[432>)
      //          - 2MW([13]-<13>)( gL([23]<14>+[12]<34>) + gR([34]<12>+[14]<23>) )
      //        )/(t-MZ^2)
      amplitude += + normFactor*preZ*(
				      + (gL-gR)*md*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*(s24s.v(ds2,ds4)+a24a.v(ds2,ds4))
				      - two*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*(gL*s234a.v(ds2,ds4)-gR*s432a.v(ds4,ds2))
				      - two*MW*(s13s.v(ds1a,ds3a)-a13a.v(ds1a,ds3a))*(
										      + gL*(s23s.v(ds2,ds3b)*a14a.v(ds1b,ds4)+s12s.v(ds1b,ds2)*a34a.v(ds3b,ds4))
										      + gR*(s34s.v(ds3b,ds4)*a12a.v(ds1b,ds2)+s14s.v(ds1b,ds4)*a23a.v(ds2,ds3b))
										      )
				       )/pDenZT;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble WdWd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/3
    return amp2/18.0;
  }

  
  

  //  Tests
  int test_WdWd(){
    int n=0;//Number of fails
    std::cout<<"\t* W+, d  -> W+, d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, mu=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble md=0.0075, mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      WdWd WdWdAmp = WdWd(EE,md,mu,mh,0,MW,SW,0);
      ldouble pspatial=250;
      //2ND ROW: h only
      ldouble dataCH[20] = {3.199002895303244E+01,6.232099908788318E+00,2.482487444996788E+00,1.260402394232308E+00,7.269202242222387E-01,4.534497642775587E-01,2.981121533100848E-01,2.033561287231379E-01,1.424509293353464E-01,1.017095809719861E-01,7.358737292874543E-02,5.367779604097386E-02,3.928457871040673E-02,2.869360391385163E-02,2.077979571634165E-02,1.478449361835410E-02,1.018418120865246E-02,6.610354491066208E-03,3.799116338770002E-03,1.558625488632791E-03};
      i += WdWdAmp.test_2to2_amp2([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH);
      i += WdWdAmp.test_2to2_amp2_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH);
      i += WdWdAmp.test_2to2_amp2_boosts([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH);
      i += WdWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH);
      //std::cout<<"\n# md=0.0075, mu=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.430569480128815E+01,2.897159367406667E+00,1.347273744201156E+00,7.685360521273549E-01,4.827335485104632E-01,3.211166702309883E-01,2.219063303264571E-01,1.575066043139582E-01,1.139668511505201E-01,8.360346780819858E-02,6.190698236344274E-02,4.609705894435375E-02,3.439035781982458E-02,2.560528928405755E-02,1.893742531356968E-02,1.382635669230064E-02,9.873888109781563E-03,6.792212764354121E-03,4.370195964937085E-03,2.450923582219008E-03};
      i += WdWdAmp.test_2to2_amp2([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH2);
      i += WdWdAmp.test_2to2_amp2_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH2);
      i += WdWdAmp.test_2to2_amp2_boosts([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH2);
      i += WdWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, MW=80.385, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      WdWdAmp.set_masses(md,mu,mh,MW);
      ldouble dataCH4[20] = {2.839726153097277E+01,5.966840321491049E+00,3.107310449346867E+00,2.010312567413830E+00,1.432576281403560E+00,1.079583086827000E+00,8.444719679936906E-01,6.786855145782578E-01,5.568786472551179E-01,4.645237449628551E-01,3.927124894174743E-01,3.356963396116794E-01,2.896141305459427E-01,2.517870246163494E-01,2.203059696275169E-01,1.937793770576304E-01,1.711730135071380E-01,1.517051986075421E-01,1.347763525924100E-01,1.199205269931044E-01};
      i += WdWdAmp.test_2to2_amp2([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH4);
      i += WdWdAmp.test_2to2_amp2_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH4);
      i += WdWdAmp.test_2to2_amp2_boosts([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH4);
      i += WdWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      WdWdAmp.set_masses(md,mu,mh,MW);
      ldouble dataCH3[20] = {8.098478655605966E+01,1.250610820453296E+01,5.380084480959685E+00,3.051850771476602E+00,1.975403648491651E+00,1.383875467712544E+00,1.022795421219608E+00,7.859983233036016E-01,6.222984851047857E-01,5.044090005293752E-01,4.166617379809900E-01,3.495288386688777E-01,2.969449399442449E-01,2.549014623193210E-01,2.206649652814502E-01,1.923213868531052E-01,1.684994822789103E-01,1.481971899601037E-01,1.306693617336142E-01,1.153532790961160E-01};
      i += WdWdAmp.test_2to2_amp2([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH3);
      i += WdWdAmp.test_2to2_amp2_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH3);
      i += WdWdAmp.test_2to2_amp2_boosts([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH3);
      i += WdWdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return WdWdAmp.amp2(); }, MW,md,MW,md,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}

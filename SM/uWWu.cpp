
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

//File:  SPINAS/SM/uWWu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uWWu.h"

namespace spinas {

  uWWu::uWWu(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(2.0/3.0), mu(massu), md(massd), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propd(md,0), propA(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(mu);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(mu);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    //[213>,[312>,[431>,[134>
    s213a = sproduct(SQUARE,&p2,&p1,&p3,2);
    s312a = sproduct(SQUARE,&p3,&p1,&p2,2);
    s431a = sproduct(SQUARE,&p4,&p3,&p1,2);
    s134a = sproduct(SQUARE,&p1,&p3,&p4,2);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*mu/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW+1.0;
    gR=-2.0*Qf*SW*SW;
  }
  void uWWu::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW){
    mu=massu;
    md=massd;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mu);
    p2.set_mass(MW);
    p3.set_mass(MW);
    p4.set_mass(mu);
    propd.set_mass(md);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    pred = 2.0*e*e/(2.0*MW*MW*SW*SW);
    preh = 2.0*e*e*mu/(4.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
  }
  void uWWu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[213>,[312>,[431>,[134>
    s213a.update();
    s312a.update();
    s431a.update();
    s134a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenhU=proph.denominator(propUP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDendT=propd.denominator(propTP);
    pDendS=propd.denominator(propSP);
    pDenZU=propZ.denominator(propUP);
    pDenAU=propA.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uWWu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds2a, ds2b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);
      
      //S-Channel h
      //preh = e*e*mu/(4.0*MW*MW*SW*SW);
      //uUW-W+ all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //uW+W-U: 2<->4
      //preh [23] <23> ([14]+<14>)/(u-Mh^2)
      //34 out:
      //- preh [23] <23> ([14]-<14>)/(u-Mh^2)
      amplitude += - normFactor*preh*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))/pDenhU;
      
      //T-Channel d
      //pred = e*e/(2.0*MW*MW*SW*SW);
      //uUW-W+ all ingoing:
      //- pred [24] <13> ([314>+MW <34>))/t
      //uW+W-U: 2<->4
      //+ pred [24] <13> ([312>-MW <23>))/t
      //34 out:
      //- pred [24] <13> ([312>+MW <23>))/t
      amplitude += - normFactor*pred*s24s.v(ds2a,ds4)*a13a.v(ds1,ds3a)*(MW*a23a.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDendT;

      //S-Channel A
      //uUW-W+ all in:
      //- 2e^2/MW*( <13>[24] + [13]<24> + <14>[23] + [14]<23> )( <34> + [34] )/s
      //- 2e^2/MW/MW*[34]<34>([231>+[132>)/s
      //uW+W-U: 2<->4
      //- 2e^2/MW*( <13>[24] + [13]<24> + <12>[34] + [12]<34> )( <23> + [23] )/u
      //- 2e^2/MW/MW*[23]<23>([431>+[134>)/u
      //34 out:
      //- 2e^2/MW*(-<13>[24] - [13]<24> + <12>[34] + [12]<34> )(-<23> + [23] )/u
      //- 2e^2/MW/MW*[23]<23>([431>-[134>)/u
      amplitude +=
	- normFactor*e*e*Qf*two/MW*(
				- a13a.v(ds1,ds3a)*s24s.v(ds2a,ds4) - s13s.v(ds1,ds3a)*a24a.v(ds2a,ds4)
				+ a12a.v(ds1,ds2a)*s34s.v(ds3a,ds4) + s12s.v(ds1,ds2a)*a34a.v(ds3a,ds4)
				)*(s23s.v(ds2b,ds3b)-a23a.v(ds2b,ds3b))/pDenAU
	- normFactor*e*e*Qf*two/MW/MW*(
				s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s431a.v(ds4,ds1)-s134a.v(ds1,ds4))
				)/pDenAU;

      //S-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=pred
      //uUW-W+ all ingoing:
      // - preZ ( (gL-gR)Mu[34]<34>([12]-<12>)
      //          + 2[34]<34>(gL[231>+gR[132>)
      //          + 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //        )/(s-MZ^2)
      //uW+W-U: 2<->4
      // - preZ ( (gL-gR)Mu[23]<23>([14]-<14>)
      //          + 2[23]<23>(gL[431>+gR[134>)
      //          + 2MW([23]+<23>)( gL([34]<12>+[24]<13>) + gR([13]<24>+[12]<34>) )
      //        )/(u-MZ^2)
      //34 out:
      // - preZ (-(gL-gR)Mu[23]<23>([14]+<14>)
      //          + 2[23]<23>(gL[431>-gR[134>)
      //          + 2MW([23]-<23>)( gL([34]<12>-[24]<13>) - gR([13]<24>-[12]<34>) )
      //        )/(u-MZ^2)
      amplitude += - normFactor*preZ*(
				      - (gL-gR)*mu*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)+a14a.v(ds1,ds4))
				      + two*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(gL*s431a.v(ds4,ds1)-gR*s134a.v(ds1,ds4))
				      + two*MW*(s23s.v(ds2a,ds3a)-a23a.v(ds2a,ds3a))*(
										      + gL*(s34s.v(ds3b,ds4)*a12a.v(ds1,ds2b)-s24s.v(ds2b,ds4)*a13a.v(ds1,ds3b))
										      - gR*(s13s.v(ds1,ds3b)*a24a.v(ds2b,ds4)-s12s.v(ds1,ds2b)*a34a.v(ds3b,ds4)) )
				      )/pDenZU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uWWu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
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
  int test_uWWu(){
    int n=0;//Number of fails
    std::cout<<"\t* u , W+ -> W+, u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, mh=125, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      uWWu uWWuAmp = uWWu(EE,mu,md,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.345360740048891E+00,4.859624800731109E-01,3.292564155990730E-01,2.699862387031694E-01,2.444656921360539E-01,2.360527095869248E-01,2.390129088625200E-01,2.514942860445614E-01,2.736277891206455E-01,3.070722979304301E-01,3.552275614372647E-01,4.240669421645826E-01,5.239788231782706E-01,6.737420854304733E-01,9.097165462691197E-01,1.309649422945971E+00,2.065286598476736E+00,3.765483234736999E+00,9.086265267832191E+00,5.418453935907692E+01};
      i += uWWuAmp.test_2to2_amp2([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH);
      //i += uWWuAmp.test_2to2_amp2_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH);
      //i += uWWuAmp.test_2to2_amp2_boosts([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH);
      //i += uWWuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.033040686995564E+00,5.572889256942820E-01,3.402337385817036E-01,2.690811144363999E-01,2.398385347395948E-01,2.297222172836701E-01,2.313583128119445E-01,2.422080012745883E-01,2.619027289964149E-01,2.915102628560960E-01,3.335299190193242E-01,3.924028267386383E-01,4.757275894744924E-01,5.968682921242393E-01,7.808529498196600E-01,1.079306297828339E+00,1.615296416818430E+00,2.759009301100405E+00,6.234961143406906E+00,4.043807000875533E+01};
      //i += uWWuAmp.test_2to2_amp2([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH2);
      //i += uWWuAmp.test_2to2_amp2_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH2);
      //i += uWWuAmp.test_2to2_amp2_boosts([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH2);
      //i += uWWuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, mh=125, MW=80.385, pspatial=95\n";
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uWWuAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH4[20] = {4.333295843011305E+02,1.122897480553641E+01,4.203061910197297E+00,2.561233826644194E+00,1.860809663938117E+00,1.481858624286367E+00,1.250460763494224E+00,1.100098123034463E+00,1.000650809304884E+00,9.374044467491455E-01,9.035908648806030E-01,8.977265228156570E-01,9.233755409606215E-01,9.912286190120589E-01,1.125512720534533E+00,1.382257703305132E+00,1.909295712000470E+00,3.204132618722058E+00,7.916876599578562E+00,6.764952076595783E+01};
      //i += uWWuAmp.test_2to2_amp2([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH4);
      //i += uWWuAmp.test_2to2_amp2_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH4);
      //i += uWWuAmp.test_2to2_amp2_boosts([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH4);
      //i += uWWuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uWWuAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH3[20] = {4.058983591753044E+01,4.700667129258085E+00,2.481955372081903E+00,1.730077929994683E+00,1.362013236269201E+00,1.150890096886114E+00,1.020691536631484E+00,9.393828514234707E-01,8.918072676739532E-01,8.707315488622759E-01,8.734870076903062E-01,9.008332462861118E-01,9.570882196037531E-01,1.051535932304384E+00,1.202039796276925E+00,1.443860502978666E+00,1.853758972778155E+00,2.632590855376165E+00,4.540750673417542E+00,1.655533456854371E+01};
      //i += uWWuAmp.test_2to2_amp2([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH3);
      //i += uWWuAmp.test_2to2_amp2_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH3);
      //i += uWWuAmp.test_2to2_amp2_boosts([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH3);
      //i += uWWuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uWWuAmp.amp2(); }, mu,MW,MW,mu,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}

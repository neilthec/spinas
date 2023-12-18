
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

//File:  SPINAS/SM/ZWWh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ZWWh.h"

namespace spinas {

  ZWWh::ZWWh(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW, const ldouble& widthZ):
    e(echarge), mh(massh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW), WZ(widthZ) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;
    propW = propagator(MW,WW);
    propZ = propagator(MZ,WZ);
    p1=particle(MZ);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(mh);
    //<23>,[23],<12>,[12],<13>,[13]
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[131>,[121>,[231>,[242>,[321>,[343>,[143>,[142>,[341>,[241>
    s131a = sproduct(SQUARE,&p1,&p3,&p1);
    s121a = sproduct(SQUARE,&p1,&p2,&p1);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    s242a = sproduct(SQUARE,&p2,&p4,&p2);
    s321a = sproduct(SQUARE,&p3,&p2,&p1);
    s343a = sproduct(SQUARE,&p3,&p4,&p3);
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s142a = sproduct(SQUARE,&p1,&p4,&p2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    //Couplings
    pre = sqrt2*e*e/(2.0*MW*MW*SW*SW);
    preS = pre;
    preTU = pre/(MZ*MZ);
  }
  void ZWWh::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;
    p1.set_mass(MZ);
    p2.set_mass(MW);
    p3.set_mass(MW);
    p4.set_mass(mh);
    propW.set_mass(MW);
    //Couplings
    pre = sqrt2*e*e/(2.0*MW*MW*SW*SW);
    preS = pre;
    preTU = pre/(MZ*MZ);
  }
  void ZWWh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<23>,[23],<12>,[12],<13>,[13]
    s23s.update();
    a23a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    //[131>,[121>,[231>,[242>,[321>,[343>,[143>,[142>,[341>,[241>
    s131a.update();
    s121a.update();
    s231a.update();
    s242a.update();
    s321a.update();
    s343a.update();
    s143a.update();
    s142a.update();
    s341a.update();
    s241a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);
    pDenT=propW.denominator(propTP);
    pDenU=propZ.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ZWWh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds1a, ds1b, ds3a, ds3b, ds2a, ds2b;

    //Symmetrize the Z- & W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds1,ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);

      //pre = e^2 / 2 MW^2 SW^2

      //S-Channel Z
      //preS = pre;
      //ZhW-W+ All ingoing:
      //- ( [34]<34>([131>-[141>) + 2MW([13]<14>+[14]<13>)([34]+<34>) )/(s-MZ^2)
      //ZW+W-h: 2<->4
      //- ( [23]<23>([131>-[121>) - 2MW([13]<12>+[12]<13>)([23]+<23>) )/(u-MZ^2)
      //34 out:
      //+ (-[23]<23>([131>+[121>) + 2MW([13]<12>-[12]<13>)([23]-<23>) )/(u-MZ^2)
      amplitude += normFactor*preS*(
				    - s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s131a.v(ds1a,ds1b)+s121a.v(ds1a,ds1b))
				    + 2*MW*(s13s.v(ds1a,ds3a)*a12a.v(ds1b,ds2a)-s12s.v(ds1a,ds2a)*a13a.v(ds1b,ds3a))*(s23s.v(ds2b,ds3b)-a23a.v(ds2b,ds3b))
				     )/pDenU;
      
      //T-Channel W
      //preTU = pre / MZ^2
      //ZhW-W+ All ingoing:
      // - ( 2MW^2[13]<34>[431> - MZ^2[13]<13>[424> + 2MW^3[14]<13><34> + 2MW^2MZ[13][14]<34> + 2MW^2MZ[34]<13><14> )/(t-MW^2)
      //ZW+W-h: 2<->4:
      // + ( 2MW^2[13]<23>[231> + MZ^2[13]<13>[242> + 2MW^3[12]<13><23> + 2MW^2MZ[13][12]<23> + 2MW^2MZ[23]<13><12> )/(t-MW^2)
      //34 out:
      // + ( 2MW^2[13]<23>[231> + MZ^2[13]<13>[242> + 2MW^3[12]<13><23> - 2MW^2MZ[13][12]<23> - 2MW^2MZ[23]<13><12> )/(t-MW^2)
      amplitude += normFactor*preTU*(
				     2.0*MW*MW*s13s.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)*s231a.v(ds2b,ds1b)
				     + MZ*MZ*s13s.v(ds1a,ds3a)*a13a.v(ds1b,ds3b)*s242a.v(ds2a,ds2b)
				     + 2.0*MW*MW*MW*s12s.v(ds1a,ds2a)*a13a.v(ds1b,ds3a)*a23a.v(ds2b,ds3b)
				     - 2.0*MW*MW*MZ*s13s.v(ds1a,ds3a)*s12s.v(ds1b,ds2a)*a23a.v(ds2b,ds3b)
				     - 2.0*MW*MW*MZ*a13a.v(ds1a,ds3a)*a12a.v(ds1b,ds2a)*s23s.v(ds2b,ds3b)
				     )/pDenT;

      //U-Channel W
      //ZhW-W+ All ingoing:
      // - ( 2MW^2[14]<34>[341> + MZ^2[14]<14>[323> + 2MW^3[13]<14><34> + 2MW^2MZ[13][14]<34> + 2MW^2MZ[34]<13><14> )/(u-MW^2)
      //ZW+W-h: 2<->4:
      // - (-2MW^2[12]<23>[321> + MZ^2[12]<12>[343> - 2MW^3[13]<12><23> - 2MW^2MZ[13][12]<23> - 2MW^2MZ[23]<13><12> )/(s-MW^2)
      //34 out:
      // - ( 2MW^2[12]<23>[321> + MZ^2[12]<12>[343> + 2MW^3[13]<12><23> + 2MW^2MZ[13][12]<23> + 2MW^2MZ[23]<13><12> )/(s-MW^2)
      amplitude += - normFactor*preTU*(
				       + 2.0*MW*MW*s12s.v(ds1a,ds2a)*a23a.v(ds2b,ds3a)*s321a.v(ds3b,ds1b)
				       + MZ*MZ*s12s.v(ds1a,ds2a)*a12a.v(ds1b,ds2b)*s343a.v(ds3a,ds3b)
				       + 2.0*MW*MW*MW*s13s.v(ds1a,ds3a)*a12a.v(ds1b,ds2a)*a23a.v(ds2b,ds3b)
				       + 2.0*MW*MW*MZ*s12s.v(ds1a,ds2a)*s13s.v(ds1b,ds3a)*a23a.v(ds2b,ds3b)
				       + 2.0*MW*MW*MZ*a12a.v(ds1a,ds2a)*a13a.v(ds1b,ds3a)*s23s.v(ds2b,ds3b)
				       )/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ZWWh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3^2=1/9
    return amp2/9.0;
  }


  



  //  Tests
  int test_ZWWh(){
    int n=0;//Number of fails
    std::cout<<"\t* Z , W+ -> W+, h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW = std::sqrt(1.0-SW*SW);
      ldouble MZ = MW/CW;
      ZWWh ZWWhAmp = ZWWh(EE,mh,MW,SW,0,0);
      ldouble pspatial=250;
      //2nd Z only
      //3rd W only
      ldouble dataCH[20] = {1.650679054460746E+01,4.149265763443914E+00,1.825527626116234E+00,1.016918302864870E+00,6.489834392601030E-01,4.553436722959721E-01,3.450246203115704E-01,2.805578812003502E-01,2.447392226529418E-01,2.294593510257980E-01,2.315082474207477E-01,2.511150817884552E-01,2.919970641230435E-01,3.629295788527520E-01,4.820818661899239E-01,6.883016199829451E-01,1.073830570753030E+00,1.899704122118246E+00,4.155735000676124E+00,1.457523899289568E+01};
      i += ZWWhAmp.test_2to2_amp2([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH);
      i += ZWWhAmp.test_2to2_amp2_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH);
      i += ZWWhAmp.test_2to2_amp2_boosts([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH);
      i += ZWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {3.442745250725873E+00,1.885184503820596E+00,1.198606336479058E+00,8.394023818468438E-01,6.309820125108367E-01,5.022221606134775E-01,4.201918490811433E-01,3.681738315228095E-01,3.372382908357218E-01,3.226455072925095E-01,3.222299310405782E-01,3.357288181528834E-01,3.646803947977510E-01,4.128312388294069E-01,4.872525102186058E-01,6.007873429364929E-01,7.774630377081975E-01,1.065427259176863E+00,1.571912561185825E+00,2.576061569565989E+00};
      i += ZWWhAmp.test_2to2_amp2([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH2);
      i += ZWWhAmp.test_2to2_amp2_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH2);
      i += ZWWhAmp.test_2to2_amp2_boosts([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH2);
      i += ZWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=60\n";
      pspatial = 60;
      ldouble dataCH4[20] = {4.901908090684942E-01,4.743061788570282E-01,4.603452792176039E-01,4.481284858350041E-01,4.375043011919574E-01,4.283450909885342E-01,4.205436571728521E-01,4.140104799180147E-01,4.086715002423873E-01,4.044663447920651E-01,4.013469171706001E-01,3.992762979747758E-01,3.982279097375769E-01,3.981849143008509E-01,3.991398195057738E-01,4.010942800946172E-01,4.040590848432681E-01,4.080543285984926E-01,4.131097744460565E-01,4.192654180448243E-01};
      i += ZWWhAmp.test_2to2_amp2([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH4);
      i += ZWWhAmp.test_2to2_amp2_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH4);
      i += ZWWhAmp.test_2to2_amp2_boosts([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH4);
      i += ZWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH4);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=95\n";
      mh = 0.0005;
      pspatial = 95;
      ZWWhAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {1.996856089730392E+00,1.255516386421362E+00,8.738199821179560E-01,6.531768458036358E-01,5.159149330819059E-01,4.266432383351487E-01,3.674615652003468E-01,3.286606534423537E-01,3.047786915063659E-01,2.927841401629772E-01,2.912087073771666E-01,2.997650454868821E-01,3.192686231523025E-01,3.518275058914491E-01,4.013813275487333E-01,4.748455605290416E-01,5.844957519798349E-01,7.532118965885626E-01,1.027150022905959E+00,1.510723061306611E+00};
      i += ZWWhAmp.test_2to2_amp2([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH3);
      i += ZWWhAmp.test_2to2_amp2_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH3);
      i += ZWWhAmp.test_2to2_amp2_boosts([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH3);
      i += ZWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZWWhAmp.amp2(); }, MZ,MW,MW,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }
  
  

}

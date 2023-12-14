
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

//File:  SPINAS/SM/neneZh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/neneZh.h"

namespace spinas {

  neneZh::neneZh(const ldouble& echarge, const ldouble& massh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mh(massh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(0);
    p2=particle(0);
    p3=particle(MZ);
    p4=particle(mh);
    //[23], <13>
    s23s = sproduct(SQUARE,&p2,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //Couplings
    preZ = sqrt2*two*e*e/(4.0*CW*CW*SW*SW);
  }
  void neneZh::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p3.set_mass(MZ);
    p4.set_mass(mh);
    propZ.set_mass(MZ);
  }
  void neneZh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);    
    //[23], <13>
    s23s.update();
    a13a.update();
    //Propagator Momentum
    ldouble propSP[4];
    for(int j=0;j<4;j++)
      propSP[j] = mom1[j]+mom2[j];
    pDenS=propZ.den(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble neneZh::amp(const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    
    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //S-Channel Z
      //preZ = 2.0*e*e*MZ*MZ/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preZ [23] <13> /(s-MZ^2)
      //34 outgoing:
      //- preZ [23] <13> /(s-MZ^2)
      amplitude += - normFactor*preZ*s23s.v(ds3a)*a13a.v(ds3b)/pDenS;

    }

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble neneZh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=2){
      M = amp(j3);
      amp2 += std::pow(std::abs(M),2);
    }
    //Nothing to average over.
    return amp2;
  }
  



  //  Tests
  int test_neneZh(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Ne -> Z , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      neneZh neneZhAmp = neneZh(EE,mh,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.342098465167624E-03,1.043647864158071E-02,1.318703879839235E-02,1.559377893560253E-02,1.765669905321126E-02,1.937579915121854E-02,2.075107922962435E-02,2.178253928842872E-02,2.247017932763163E-02,2.281399934723308E-02,2.281399934723308E-02,2.247017932763163E-02,2.178253928842872E-02,2.075107922962436E-02,1.937579915121854E-02,1.765669905321127E-02,1.559377893560254E-02,1.318703879839236E-02,1.043647864158072E-02,7.342098465167634E-03};
      i += neneZhAmp.test_2to2_amp2([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH);
      i += neneZhAmp.test_2to2_amp2_rotations([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH);
      i += neneZhAmp.test_2to2_amp2_boosts([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH);
      i += neneZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {2.877665381853260E-02,2.995641110461577E-02,3.100508424780082E-02,3.192267324808774E-02,3.270917810547653E-02,3.336459881996718E-02,3.388893539155970E-02,3.428218782025410E-02,3.454435610605036E-02,3.467544024894850E-02,3.467544024894850E-02,3.454435610605036E-02,3.428218782025410E-02,3.388893539155970E-02,3.336459881996718E-02,3.270917810547653E-02,3.192267324808774E-02,3.100508424780082E-02,2.995641110461577E-02,2.877665381853260E-02};
      i += neneZhAmp.test_2to2_amp2([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH2);
      i += neneZhAmp.test_2to2_amp2_rotations([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH2);
      i += neneZhAmp.test_2to2_amp2_boosts([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH2);
      i += neneZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH2);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=125.1\n";
      mh = 0.0005;
      pspatial = 125.1;
      neneZhAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {3.007320575205575E-02,3.364659737695089E-02,3.682294548796879E-02,3.960225008510946E-02,4.198451116837289E-02,4.396972873775908E-02,4.555790279326803E-02,4.674903333489974E-02,4.754312036265422E-02,4.794016387653146E-02,4.794016387653146E-02,4.754312036265422E-02,4.674903333489974E-02,4.555790279326803E-02,4.396972873775908E-02,4.198451116837289E-02,3.960225008510946E-02,3.682294548796879E-02,3.364659737695089E-02,3.007320575205575E-02};
      i += neneZhAmp.test_2to2_amp2([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH3);
      i += neneZhAmp.test_2to2_amp2_rotations([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH3);
      i += neneZhAmp.test_2to2_amp2_boosts([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH3);
      i += neneZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neneZhAmp.amp2(); }, 0,0,MZ,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

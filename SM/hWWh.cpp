
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

//File:  SPINAS/SM/hWWh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hWWh.h"

namespace spinas {

  hWWh::hWWh(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WW(widthW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    CW = std::sqrt(1.0-sinW*sinW);
    propW = propagator(MW,WW);
    proph = propagator(mh,wh);  
    p1=particle(mh);
    p2=particle(MW);
    p3=particle(MW);
    p4=particle(mh);
    //<23>,[23]
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    //[312>,[213>
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    //Couplings
    pre = 2.0*e*e/(4.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hWWh::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    p1.set_mass(mh);
    p2.set_mass(MW);
    p3.set_mass(MW);
    p4.set_mass(mh);
    propW.set_mass(MW);
    proph.set_mass(mh);
    //Couplings
    pre = 2.0*e*e/(4.0*MW*MW*SW*SW);
    preS = 3.0*pre*mh*mh;
  }
  void hWWh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<23>,[23]
    s23s.update();
    a23a.update();
    //[314>,[413>
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propW.denominator(propTP);
    pDenU=proph.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hWWh::amp(const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1;
    int ds3a, ds3b, ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);


      //pre = sqrt2*e*e/(4.0*MW*MW*SW*SW);
      
      //4-Point
      //hhW-W+ all in:  - [34]<34>
      //hW+W-h: 2<->4:  - [23]<23>
      //34 out:         + [23]<23>
      amplitude += + normFactor*pre*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b);
          
      //U-Channel h
      //preS = 3.0*pre*mh*mh;
      //hhW-W+ all ingoing: 
      //+ preS [34] <34> /(s-Mh^2)
      //hW+W-h: 2<->4:
      //+ preS [23] <23> /(u-Mh^2)
      //34 out:
      //- preS [23] <23> /(u-Mh^2)
      amplitude += - normFactor*preS*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenU;
      
      //T-Channel Z
      //hhW-W+ all ingoing:
      //+ pre ( 2MW^2[34]<34> + MW([34][314>+<34>[413>) + [314>[413> )/(t-MW^2)
      //hW+W-h: 2<->4:
      //+ pre ( 2MW^2[23]<23> - MW([23][312>+<23>[213>) + [312>[213> )/(t-MW^2)
      //34 out:
      //- pre ( 2MW^2[23]<23> + MW([23][312>+<23>[213>) + [312>[213> )/(t-MW^2)
      amplitude +=
	- normFactor*pre*2.0*MW*MW*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenT
      	- normFactor*pre*MW*(s23s.v(ds2a,ds3a)*s312a.v(ds3b,ds2b)+a23a.v(ds2a,ds3a)*s213a.v(ds2b,ds3b))/pDenT
      	- normFactor*pre*s312a.v(ds3a,ds2a)*s213a.v(ds2b,ds3b)/pDenT;
      
      //S-Channel Z
      //hhW-W+ all ingoing:
      //+ pre ( 2MW^2[34]<34> - MW([34][413>+<34>[314>) + [314>[413> )/(u-MW^2)
      //hW+W-h: 2<->4:
      //+ pre ( 2MW^2[32]<23> + MW([23][213>+<23>[312>) + [312>[213> )/(s-MW^2)
      //34 out:
      //- pre ( 2MW^2[32]<23> + MW([23][213>+<23>[312>) + [312>[213> )/(s-MW^2)
      amplitude +=
	- normFactor*pre*2.0*MW*MW*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)/pDenS
      	- normFactor*pre*MW*(s23s.v(ds2a,ds3a)*s213a.v(ds2b,ds3b)+a23a.v(ds2a,ds3a)*s312a.v(ds3b,ds2b))/pDenS
      	- normFactor*pre*s312a.v(ds3a,ds2a)*s213a.v(ds2b,ds3b)/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hWWh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2){
	M = amp(j2,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/3
    return amp2/3.0;
  }

  



  //  Tests
  int test_hWWh(){
    int n=0;//Number of fails
    std::cout<<"\t* h , W+ -> W+, h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      hWWh hWWhAmp = hWWh(EE,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.681773240591237E+00,2.001674794658378E+00,9.449006628722678E-01,5.723711252288808E-01,3.988506896480206E-01,3.042138277868801E-01,2.471352106707728E-01,2.102854155919365E-01,1.853623880672921E-01,1.679946779096416E-01,1.557220427080687E-01,1.471048248250286E-01,1.413051187625992E-01,1.378949228304824E-01,1.368053479311971E-01,1.384247312571526E-01,1.439945077739038E-01,1.569173431219364E-01,1.879176738819908E-01,2.852586582594986E-01};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {2.226551174616142E+00,1.154791532709749E+00,7.331990841642416E-01,5.251102334540041E-01,4.077165050090676E-01,3.355633366348908E-01,2.886644519605909E-01,2.571345423622791E-01,2.356630532944946E-01,2.212368443228950E-01,2.121054987446865E-01,2.072812138834188E-01,2.062959223105020E-01,2.091018898003663E-01,2.160753816727608E-01,2.281287446371894E-01,2.469875560412330E-01,2.757869241452156E-01,3.203833817004011E-01,3.924858090529695E-01};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH2);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH2);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH2);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=95\n";
      pspatial = 1;
      ldouble dataCH4[20] = {9.787896963465105E-01,9.787570644633118E-01,9.787244380486615E-01,9.786918171016848E-01,9.786592016215065E-01,9.786265916072522E-01,9.785939870580471E-01,9.785613879730169E-01,9.785287943512877E-01,9.784962061919852E-01,9.784636234942359E-01,9.784310462571663E-01,9.783984744799029E-01,9.783659081615728E-01,9.783333473013026E-01,9.783007918982198E-01,9.782682419514520E-01,9.782356974601267E-01,9.782031584233716E-01,9.781706248403148E-01};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH4);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH4);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH4);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH4);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=95\n";
      mh = 5;
      pspatial = 95;
      hWWhAmp.set_masses(mh,MW);
      ldouble dataCH6[20] = {5.008563155448694E-01,2.711420503346680E-01,1.651139227224619E-01,1.103999072730748E-01,8.010493028199926E-02,6.257895023012340E-02,5.220126895261586E-02,4.603463554929923E-02,4.244694873063239E-02,4.048627363969884E-02,3.957513831406459E-02,3.935145782510715E-02,3.958123996002667E-02,4.010813632843805E-02,4.082218003681976E-02,4.163725805366606E-02,4.246724084611152E-02,4.317599381739845E-02,4.336240423451299E-02,3.979473745411363E-02};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH6);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH6);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH6);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH6);
      //std::cout<<"\n# mh=5, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH7[20] = {6.016473550376888E+00,1.383359251416995E+00,5.759146399771949E-01,3.108265115603657E-01,1.960071796077275E-01,1.377732797535955E-01,1.050671343064102E-01,8.534867348423211E-02,7.283331175479468E-02,6.458264980134022E-02,5.899086193173564E-02,5.512855049407270E-02,5.243231522847404E-02,5.054720785842814E-02,4.924192365885802E-02,4.836052375020735E-02,4.779265249169617E-02,4.745017493473591E-02,4.722444646113974E-02,4.656529099795843E-02};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH7);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH7);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH7);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH7);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=125.1\n";
      mh = 0.0005;
      pspatial = 95;
      hWWhAmp.set_masses(mh,MW);
      ldouble dataCH3[20] = {5.003211533494921E-01,2.705205089694293E-01,1.645904131775197E-01,1.099971195750472E-01,7.981346490697251E-02,6.238507263606435E-02,5.209325112196023E-02,4.600396650021147E-02,4.248834298115188E-02,4.059753905181460E-02,3.975723467592578E-02,3.960891908529168E-02,3.992318366807508E-02,4.055031208496723E-02,4.139109861551038E-02,4.237907684118517E-02,4.346940672218055E-02,4.463176363455228E-02,4.584569592797596E-02,4.709754036847440E-02};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH3);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH3);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH3);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH3);
      //std::cout<<"\n# mh=0.0005, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH5[20] = {6.013955837855818E+00,1.382291243996661E+00,5.753151220558532E-01,3.104371622114474E-01,1.957334626986550E-01,1.375717665380038E-01,1.049148245940374E-01,8.523221224376645E-02,7.274449352280909E-02,6.451631683873246E-02,5.894385571143065E-02,5.509919685655410E-02,5.242023898036451E-02,5.055343155593812E-02,4.926934841040964E-02,4.841516006369568E-02,4.788678194581249E-02,4.761198162169108E-02,4.753982310050738E-02,4.763389704789561E-02};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH5);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH5);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH5);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH5);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=95\n";
      mh = 0;
      pspatial = 95;
      hWWhAmp.set_masses(mh,MW);
      ldouble dataCH8[20] = {5.003211533441732E-01,2.705205089632275E-01,1.645904131722942E-01,1.099971195710288E-01,7.981346490406853E-02,6.238507263413788E-02,5.209325112089429E-02,4.600396649992132E-02,4.248834298158500E-02,4.059753905294972E-02,3.975723467777334E-02,3.960891908789822E-02,3.992318367153370E-02,4.055031208943870E-02,4.139109862126522E-02,4.237907684869590E-02,4.346940673234566E-02,4.463176364937227E-02,4.584569595348351E-02,4.709754044697265E-02};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH8);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH8);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH8);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH8);
      //std::cout<<"\n# mh=0, MW=80.385, pspatial=250\n";
      pspatial = 250;
      ldouble dataCH9[20] = {6.013955837830647E+00,1.382291243985982E+00,5.753151220498599E-01,3.104371622075562E-01,1.957334626959204E-01,1.375717665359915E-01,1.049148245925168E-01,8.523221224260497E-02,7.274449352192355E-02,6.451631683807302E-02,5.894385571096334E-02,5.509919685626414E-02,5.242023898024720E-02,5.055343155600450E-02,4.926934841068815E-02,4.841516006424718E-02,4.788678194676133E-02,4.761198162332052E-02,4.753982310368522E-02,4.763389705872388E-02};
      i += hWWhAmp.test_2to2_amp2([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH9);
      i += hWWhAmp.test_2to2_amp2_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH9);
      i += hWWhAmp.test_2to2_amp2_boosts([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH9);
      i += hWWhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hWWhAmp.amp2(); }, mh,MW,MW,mh,pspatial,dataCH9);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  


}

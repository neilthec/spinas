
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

//File:  SPINAS/SM/enAW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/enAW.h"

namespace spinas {

  enAW::enAW(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), me(masse), MW(massW), SW(sinW), WW(widthW) {
    //constexpr ldouble sqrt2 = std::sqrt(2);
    propW = propagator(MW,WW);
    prope = propagator(me,0);
    p1=particle(me);
    p2=particle(0);
    p3=particle(0);
    p4=particle(MW);
    //[34],<12>,<24>,[13],[314>
    s34s = sproduct(SQUARE,&p3,&p4);
    a12a = sproduct(ANGLE,&p1,&p2);
    a24a = sproduct(ANGLE,&p2,&p4);
    s13s = sproduct(SQUARE,&p1,&p3);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    //<23>,[14],<34>,[413>
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    //prefactor
    pre = std::sqrt(2.0)*e*e/(MW*SW);
  }
  void enAW::set_masses(const ldouble& masse, const ldouble& massW){
    //constexpr ldouble sqrt2 = std::sqrt(2);
    me=masse;
    MW=massW;
    p1.set_mass(me);
    p4.set_mass(MW);
    propW.set_mass(MW);
    prope.set_mass(me);
    pre = std::sqrt(2.0)*e*e/(MW*SW);
  }
  void enAW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //[34],<12>,<24>,[13],[314>
    s34s.update();
    a12a.update();
    a24a.update();
    s13s.update();
    s314a.update();
    //<23>,[14],<34>,[413>
    a23a.update();
    s14s.update();
    a34a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
    }
    pDenS=propW.denominator(propSP);
    pDenT=prope.denominator(propTP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble enAW::amp(const int& ds1, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds4);
    ldouble normFactor=get_spin_normalization(ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds4,ds4a,ds4b, i);

      if(ds3>0){
      
	//ST Diagram
	//pre = sqrt(2)*e*e/(MW*SW);
	//EnAW- all ingoing:
	// -(me[34]<12>+MW<24>[13])(MW[34]-[314>))/((s-MW^2)(t-me^2))
	//34 out:
	// -(me[34]<12>-MW<24>[13])(MW[34]+[314>))/((s-MW^2)(t-me^2))
	amplitude += - normFactor*pre*(me*s34s.v(ds4a)*a12a.v(ds1)-MW*a24a.v(ds4a)*s13s.v(ds1))*(MW*s34s.v(ds4b)+s314a.v(ds4b))/pDenS/pDenT;
	
	
      }
      else if(ds3<0){

	//ST Diagram
	//pre = sqrt(2)*e*e/(MW*SW);
	//EnAW- all ingoing:
	// (<23>[14](-me^2 <34>+MW[413>))/((s-MW^2)(t-Me^2))
	//34 out:
	// <23>[14](me^2 <34>+MW[413>)/((s-MW^2)(t-Me^2))
	amplitude += normFactor*pre*a23a.v()*s14s.v(ds1,ds4a)*(me*me*a34a.v(ds4b)+MW*s413a.v(ds4b))/pDenS/pDenT;

	
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble enAW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-2;j4<=2;j4+=2){
	  M = amp(j1,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    return amp2/2.0;
  }



  



  //  Tests
  int test_enAW(){
    int n=0;//Number of fails
    std::cout<<"\t* E , ne -> A , W+      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,me=0.0005,MW=80.385, SW=0.474;
      enAW enAWAmp = enAW(EE,me,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.682652587080813E+00,4.845025681572092E-01,2.509642793518078E-01,1.548599571283782E-01,1.042809610905685E-01,7.416129840275455E-02,5.486020258964682E-02,4.187864929609518E-02,3.283480209878003E-02,2.634970306840871E-02,2.156971158208489E-02,1.793810549444025E-02,1.507631374662470E-02,1.271793448241011E-02,1.067005601280325E-02,8.789550310486848E-03,6.967987051735442E-03,5.121719986818159E-03,3.185188541621299E-03,1.106280464926496E-03};
      i += enAWAmp.test_2to2_amp2([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH);
      i += enAWAmp.test_2to2_amp2_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH);
      i += enAWAmp.test_2to2_amp2_boosts([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH);
      i += enAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH);
      //std::cout<<"\n# MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {3.041535675541980E+00,9.142348272538774E-01,4.948663723929549E-01,3.191205504127525E-01,2.242962767274148E-01,1.660205282620704E-01,1.272267576790584E-01,9.995057622046494E-02,7.997554277903628E-02,6.486055987157335E-02,5.309447191321392E-02,4.369174265294812E-02,3.598220744008225E-02,2.949426892782770E-02,2.388642558919712E-02,1.890530532352087E-02,1.435895956968463E-02,1.009931372339942E-02,6.010309311032901E-03,1.999694725338483E-03};
      i += enAWAmp.test_2to2_amp2([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH2);
      i += enAWAmp.test_2to2_amp2_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH2);
      i += enAWAmp.test_2to2_amp2_boosts([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH2);
      i += enAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH2);
      me=80;
      enAWAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH3[20] = {1.174079874385891E+00,5.238170737774492E-01,3.118603307454951E-01,2.091549809524941E-01,1.499021649908431E-01,1.121706066451834E-01,8.656159821724277E-02,6.837592061319130E-02,5.500294979561043E-02,4.487633529875522E-02,3.699932052381426E-02,3.070740742376865E-02,2.554056249749350E-02,2.117019254100855E-02,1.735533054497665E-02,1.391525211622676E-02,1.071175682211316E-02,7.637355828698096E-03,4.607189803500653E-03,1.553371475233226E-03};
      i += enAWAmp.test_2to2_amp2([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH3);
      i += enAWAmp.test_2to2_amp2_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH3);
      i += enAWAmp.test_2to2_amp2_boosts([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH3);
      i += enAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH3);
      //std::cout<<"\n# MW=80.385, pspatial=81\n";
      pspatial = 1;
      ldouble dataCH4[20] = {1.613125514505544E-03,1.683401747981417E-03,1.736414181823207E-03,1.772290664485233E-03,1.791158081147125E-03,1.793142361813453E-03,1.778368489285631E-03,1.746960507053331E-03,1.699041527147541E-03,1.634733737901329E-03,1.554158411632046E-03,1.457435912266635E-03,1.344685702880264E-03,1.216026353170366E-03,1.071575546882984E-03,9.114500891055071E-04,7.357659135921240E-04,5.446380899263780E-04,3.381808306528724E-04,1.165074983843148E-04};
      i += enAWAmp.test_2to2_amp2([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH4);
      i += enAWAmp.test_2to2_amp2_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH4);
      i += enAWAmp.test_2to2_amp2_boosts([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH4);
      i += enAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH4);
      MW=1;
      enAWAmp.set_masses(me,MW);
      pspatial=250;
      ldouble dataCH5[20] = {2.613654310557742E+03,1.224008270773222E+03,7.685659156661995E+02,5.423035168376197E+02,4.070140779100245E+02,3.170168709204182E+02,2.528294957193332E+02,2.047419359255139E+02,1.673720558580072E+02,1.374961166289887E+02,1.130653961950543E+02,9.271554927912557E+01,7.550286381568523E+01,6.075375797590848E+01,4.797458175510469E+01,3.679529743814910E+01,2.693306774932553E+01,1.816799867313745E+01,1.032655286553813E+01,3.269945484727895E+00};
      i += enAWAmp.test_2to2_amp2([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH5);
      i += enAWAmp.test_2to2_amp2_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH5);
      i += enAWAmp.test_2to2_amp2_boosts([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH5);
      i += enAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH5);
      //std::cout<<"\n# MW=80.385, pspatial=81\n";
      pspatial = 1;
      ldouble dataCH6[20] = {3.382927805493410E+00,3.197487137211937E+00,3.013386532341626E+00,2.830617630158225E+00,2.649172128586750E+00,2.469041783730283E+00,2.290218409402476E+00,2.112693876661131E+00,1.936460113350010E+00,1.761509103643154E+00,1.587832887594687E+00,1.415423560689925E+00,1.244273273403821E+00,1.074374230757980E+00,9.057186918940632E-01,7.382989696311696E-01,5.721074300508850E-01,4.071364920638724E-01,2.433786269936593E-01,8.082635816493844E-02};
      i += enAWAmp.test_2to2_amp2([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH6);
      i += enAWAmp.test_2to2_amp2_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH6);
      i += enAWAmp.test_2to2_amp2_boosts([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH6);
      i += enAWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return enAWAmp.amp2(); }, me,0,0,MW,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}

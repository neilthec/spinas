
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

//File:  SPINAS/SM/ddnn.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddnn.h"

namespace spinas {
  //Constructors
  ddnn::ddnn(const ldouble& echarge, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), md(massd), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    p1(particle(md)), p2(particle(md)),
    p3(particle(0)), p4(particle(0)),
    s23s(sproduct(SQUARE,&p2,&p3)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    gLn=1.0;
    gRn=0;
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void ddnn::set_masses(const ldouble& massd, const ldouble& massW){
    md=massd;
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p2.set_mass(md);
  }
  void ddnn::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s23s.update();
    a14a.update();
    s13s.update();
    a24a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddnn::amp(const int& ds1, const int& ds2){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //gLn=1.0;
    //gRn=0;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //all in:
    //- preZ ( gLd gLn [23] <14> + gLn gRd [13] <24>)/(s-MZ^2)
    //34 out:
    //+ preZ ( gLd gLn [23] <14> + gLn gRd [13] <24>)/(s-MZ^2)
    amplitude += 
      + two*preZ*(
	        gLd*gLn*s23s.v(ds2)*a14a.v(ds1)
	      + gLn*gRd*s13s.v(ds1)*a24a.v(ds2)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ddnn::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += three*std::pow(std::abs(M),2);// Color factor 3
      }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3*1/3 = 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_ddnn(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> ne, Ne      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, MW=80.385, pspatial=250\n";
      ldouble md=0.0042, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      ddnn ddnnAmp = ddnn(0.31333,md,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.866461852521740E-03,4.380939939085495E-03,3.921808017818836E-03,3.489066088721762E-03,3.082714151794274E-03,2.702752207036373E-03,2.349180254448057E-03,2.021998294029326E-03,1.721206325780182E-03,1.446804349700623E-03,1.198792365790650E-03,9.771703740502626E-04,7.819383744794613E-04,6.130963670782457E-04,4.706443518466159E-04,3.545823287845721E-04,2.649102978921141E-04,2.016282591692417E-04,1.647362126159552E-04,1.542341582322546E-04};
      i += ddnnAmp.test_2to2_amp2([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      i += ddnnAmp.test_2to2_amp2_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      i += ddnnAmp.test_2to2_amp2_boosts([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      i += ddnnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      //std::cout<<"########### md=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {1.123590614097453E-19,1.071099877110247E-19,1.019667959116505E-19,9.692948601162282E-20,9.199805801094152E-20,8.717251190960667E-20,8.245284770761825E-20,7.783906540497625E-20,7.333116500168069E-20,6.892914649773155E-20,6.463300989312885E-20,6.044275518787258E-20,5.635838238196274E-20,5.237989147539933E-20,4.850728246818234E-20,4.474055536031180E-20,4.107971015178768E-20,3.752474684260999E-20,3.407566543277874E-20,3.073246592229392E-20};
      i += ddnnAmp.test_2to2_amp2([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      i += ddnnAmp.test_2to2_amp2_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      i += ddnnAmp.test_2to2_amp2_boosts([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      i += ddnnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      //std::cout<<"########### md=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      ddnnAmp.set_masses(md,MW);
      ldouble dataCH3[20] = {1.678332448832793E-04,1.599925860130963E-04,1.523100853058488E-04,1.447857427615367E-04,1.374195583801601E-04,1.302115321617189E-04,1.231616641062132E-04,1.162699542136429E-04,1.095364024840080E-04,1.029610089173086E-04,9.654377351354469E-05,9.028469627271618E-05,8.418377719482312E-05,7.824101627986550E-05,7.245641352784333E-05,6.682996893875661E-05,6.136168251260534E-05,5.605155424938951E-05,5.089958414910912E-05,4.590577221176419E-05};
      i += ddnnAmp.test_2to2_amp2([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH3);
      i += ddnnAmp.test_2to2_amp2_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH3);
      i += ddnnAmp.test_2to2_amp2_boosts([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH3);
      i += ddnnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH3);
      //std::cout<<"########### md=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      ddnnAmp.set_masses(md,MW);
      ldouble dataCH4[20] = {1.421529744539256E-03,1.355120197321707E-03,1.290050232932999E-03,1.226319851373133E-03,1.163929052642109E-03,1.102877836739926E-03,1.043166203666586E-03,9.847941534220864E-04,9.277616860064290E-04,8.720688014196132E-04,8.177154996616391E-04,7.647017807325068E-04,7.130276446322161E-04,6.626930913607669E-04,6.136981209181596E-04,5.660427333043940E-04,5.197269285194699E-04,4.747507065633877E-04,4.311140674361471E-04,3.888170111377482E-04};
      i += ddnnAmp.test_2to2_amp2([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH4);
      i += ddnnAmp.test_2to2_amp2_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH4);
      i += ddnnAmp.test_2to2_amp2_boosts([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH4);
      i += ddnnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddnnAmp.amp2(); }, md,md,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

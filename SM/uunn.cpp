
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

//File:  SPINAS/SM/uunn.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uunn.h"

namespace spinas {
  //Constructors
  uunn::uunn(const ldouble& echarge, const ldouble& massu, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mu(massu), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    p1(particle(mu)), p2(particle(mu)),
    p3(particle(0)), p4(particle(0)),
    a13a(sproduct(ANGLE,&p1,&p3,2)),
    s13s(sproduct(SQUARE,&p1,&p3,2)),
    a14a(sproduct(ANGLE,&p1,&p4,2)),
    s14s(sproduct(SQUARE,&p1,&p4,2)),
    a23a(sproduct(ANGLE,&p2,&p3,2)),
    s23s(sproduct(SQUARE,&p2,&p3,2)),
    a24a(sproduct(ANGLE,&p2,&p4,2)),
    s24s(sproduct(SQUARE,&p2,&p4,2)),
    s12s(sproduct(SQUARE,&p1,&p2,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    s34s(sproduct(SQUARE,&p3,&p4,2)),
    a34a(sproduct(ANGLE,&p3,&p4,2))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLn=1.0;
    gRn=0;
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void uunn::set_masses(const ldouble& massu, const ldouble& massW){
    mu=massu;
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mu);
  }
  void uunn::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a13a.update();
    s13s.update();
    a14a.update();
    s14s.update();
    a23a.update();
    s23s.update();
    a24a.update();
    s24s.update();
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uunn::amp(const int& ds1, const int& ds2){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLn=1.0;
    //gRn=0;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //all in:
    //- preZ 2( gLu gLe [23] <14> + gLe gRu [13] <24> )/(s-MZ^2)
    //34 out:
    //+ preZ 2( gLu gLe [23] <14> + gLe gRu [13] <24> )/(s-MZ^2)
    amplitude += 
      + two*preZ*(
	        gLu*gLn*s23s.v(ds2)*a14a.v(ds1)
	      + gLn*gRu*s13s.v(ds1)*a24a.v(ds2)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uunn::amp2(){
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
  int test_uunn(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> ne, Ne      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      uunn uunnAmp = uunn(0.31333,mu,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.303163710740519E-03,2.976281791494586E-03,2.669949012801782E-03,2.384165374662106E-03,2.118930877075559E-03,1.874245520042140E-03,1.650109303561849E-03,1.446522227634687E-03,1.263484292260653E-03,1.100995497439748E-03,9.590558431719712E-04,8.376653294573229E-04,7.368239562958030E-04,6.565317236874117E-04,5.967886316321488E-04,5.575946801300144E-04,5.389498691810084E-04,5.408541987851309E-04,5.633076689423818E-04,6.063102796527612E-04};
      i += uunnAmp.test_2to2_amp2([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      i += uunnAmp.test_2to2_amp2_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      i += uunnAmp.test_2to2_amp2_boosts([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      i += uunnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      //std::cout<<"########### mu=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {5.133564814199820E-20,4.813490249644133E-20,4.501660409140353E-20,4.198075292688479E-20,3.902734900288510E-20,3.615639231940449E-20,3.336788287644293E-20,3.066182067400043E-20,2.803820571207700E-20,2.549703799067263E-20,2.303831750978732E-20,2.066204426942108E-20,1.836821826957389E-20,1.615683951024577E-20,1.402790799143671E-20,1.198142371314672E-20,1.001738667537578E-20,8.135796878123910E-21,6.336654321391101E-21,4.619959005177353E-21};
      i += uunnAmp.test_2to2_amp2([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      i += uunnAmp.test_2to2_amp2_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      i += uunnAmp.test_2to2_amp2_boosts([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      i += uunnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      //std::cout<<"########### mu=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      uunnAmp.set_masses(mu,MW);
      ldouble dataCH3[20] = {7.668120664018439E-05,7.190018122932916E-05,6.724230910699710E-05,6.270759027318819E-05,5.829602472790246E-05,5.400761247113989E-05,4.984235350290048E-05,4.580024782318424E-05,4.188129543199117E-05,3.808549632932126E-05,3.441285051517451E-05,3.086335798955093E-05,2.743701875245051E-05,2.413383280387326E-05,2.095380014381918E-05,1.789692077228826E-05,1.496319468928050E-05,1.215262189479592E-05,9.465202388834494E-06,6.900936171396236E-06};
      i += uunnAmp.test_2to2_amp2([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH3);
      i += uunnAmp.test_2to2_amp2_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH3);
      i += uunnAmp.test_2to2_amp2_boosts([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH3);
      i += uunnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH3);
      //std::cout<<"########### mu=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      uunnAmp.set_masses(mu,MW);
      ldouble dataCH4[20] = {6.494816695106573E-04,6.089868924737520E-04,5.695352107002151E-04,5.311266241900466E-04,4.937611329432466E-04,4.574387369598150E-04,4.221594362397519E-04,3.879232307830572E-04,3.547301205897310E-04,3.225801056597731E-04,2.914731859931837E-04,2.614093615899628E-04,2.323886324501103E-04,2.044109985736263E-04,1.774764599605107E-04,1.515850166107635E-04,1.267366685243848E-04,1.029314157013745E-04,8.016925814173268E-05,5.845019584545929E-05};
      i += uunnAmp.test_2to2_amp2([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH4);
      i += uunnAmp.test_2to2_amp2_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH4);
      i += uunnAmp.test_2to2_amp2_boosts([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH4);
      i += uunnAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uunnAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

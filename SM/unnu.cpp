
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

//File:  SPINAS/SM/unnu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/unnu.h"

namespace spinas {
  //Constructors
  unnu::unnu(const ldouble& echarge, const ldouble& massu, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mu(massu), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propZ(MZ,WZ),
    p1(particle(mu)), p2(particle(0)),
    p3(particle(0)), p4(particle(mu)),
    a13a(sproduct(ANGLE,&p1,&p3)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3)),
    s23s(sproduct(SQUARE,&p2,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4)),
    s24s(sproduct(SQUARE,&p2,&p4)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s34s(sproduct(SQUARE,&p3,&p4)),
    a34a(sproduct(ANGLE,&p3,&p4))
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
  void unnu::set_masses(const ldouble& massu, const ldouble& massW){
    mu=massu;
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p4.set_mass(mu);
  }
  void unnu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
      propP[j] = mom1[j]-mom4[j];
    pDenUZ = propZ.den(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble unnu::amp(const int& ds1, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLn=1.0;
    //gRn=0;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //uUNn:
    //all in:
    //- preZ ( gLu gLn [23] <14> + gLn gRu [13] <24> )/(s-MZ^2)
    //unNU: 2<->4
    //+ preZ ( gLu gLn [34] <12> + gLn gRu [13] <24> )/(u-MZ^2)
    //34 out:
    //+ preZ ( gLu gLn [34] <12> - gLn gRu [13] <24> )/(u-MZ^2)
    amplitude += 
      + two*preZ*(
		  gLu*gLn*s34s.v(ds4)*a12a.v(ds1)
		  - gLn*gRu*s13s.v(ds1)*a24a.v(ds4)
		  )/pDenUZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble unnu::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += three*std::pow(std::abs(M),2);// Color factor 3
      }
    //Average over initial spins 1/2
    //Average over colors 1/3
    return amp2/6.0;
  }

  



  //  Tests
  int test_unnu(){
    int n=0;//Number of fails
    std::cout<<"\t* u , nm -> nm, u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      unnu unnuAmp = unnu(0.31333,mu,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.916056388688479E-02,2.123147452350355E-02,2.367639213386586E-02,2.658768088540233E-02,3.008765248320469E-02,3.434093203296972E-02,3.957319161077887E-02,4.610028351130355E-02,5.437491641230738E-02,6.506404058329010E-02,7.918234316775721E-02,9.833354326233971E-02,1.251715333540637E-01,1.643434798063317E-01,2.245884634071902E-01,3.239419658336019E-01,5.046570408439137E-01,8.857843894166476E-01,1.919400960133950E+00,6.719051978929588E+00};
      i += unnuAmp.test_2to2_amp2([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      i += unnuAmp.test_2to2_amp2_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      i += unnuAmp.test_2to2_amp2_boosts([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      i += unnuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      //std::cout<<"########### mu=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {5.054688591009471E-20,4.984082892004365E-20,4.913887572176173E-20,4.844102631524953E-20,4.774728070050764E-20,4.705763887753667E-20,4.637210084633720E-20,4.569066660690982E-20,4.501333615925513E-20,4.434010950337370E-20,4.367098663926614E-20,4.300596756693304E-20,4.234505228637498E-20,4.168824079759256E-20,4.103553310058636E-20,4.038692919535699E-20,3.974242908190502E-20,3.910203276023106E-20,3.846574023033567E-20,3.783355149221948E-20};
      i += unnuAmp.test_2to2_amp2([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      i += unnuAmp.test_2to2_amp2_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      i += unnuAmp.test_2to2_amp2_boosts([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      i += unnuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      //std::cout<<"########### mu=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      unnuAmp.set_masses(mu,MW);
      ldouble dataCH3[20] = {4.059142056530517E-05,4.007891690566802E-05,3.956828365087382E-05,3.905952987867516E-05,3.855266470592558E-05,3.804769728876196E-05,3.754463682278779E-05,3.704349254325747E-05,3.654427372526151E-05,3.604698968391277E-05,3.555164977453354E-05,3.505826339284376E-05,3.456683997515005E-05,3.407738899853590E-05,3.358991998105272E-05,3.310444248191192E-05,3.262096610167806E-05,3.213950048246292E-05,3.166005530812063E-05,3.118264030444385E-05};
      i += unnuAmp.test_2to2_amp2([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH3);
      i += unnuAmp.test_2to2_amp2_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH3);
      i += unnuAmp.test_2to2_amp2_boosts([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH3);
      i += unnuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH3);
      //std::cout<<"########### mu=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      unnuAmp.set_masses(mu,MW);
      ldouble dataCH4[20] = {1.843293284951736E-01,1.996320455427036E-01,2.171805743239515E-01,2.374506761923247E-01,2.610528185877115E-01,2.887810077880042E-01,3.216836857911163E-01,3.611688992433957E-01,4.091641237186144E-01,4.683658801398981E-01,5.426419463780157E-01,6.377031750062917E-01,7.622736645375712E-01,9.302324230471924E-01,1.164773217027691E+00,1.507095140994865E+00,2.036300184351735E+00,2.920712325462931E+00,4.573566188367653E+00,8.251952851891856E+00};
      i += unnuAmp.test_2to2_amp2([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH4);
      i += unnuAmp.test_2to2_amp2_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH4);
      i += unnuAmp.test_2to2_amp2_boosts([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH4);
      i += unnuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return unnuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

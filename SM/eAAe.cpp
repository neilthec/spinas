
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

//File:  SPINAS/SM/eAAe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eAAe.h"

namespace spinas {

  eAAe::eAAe(const ldouble& echarge, const ldouble& masse):
    e(echarge), me(masse), prop(masse,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(me);
    p2=particle(0);
    p3=particle(0);
    p4=particle(me);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s342a = sproduct(SQUARE,&p3,&p4,&p2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
  }
  void eAAe::set_masses(const ldouble& masse){
    me=masse;
    p1.set_mass(me);
    p4.set_mass(me);
    prop.set_mass(me);
  }
  void eAAe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s342a.update();
    s243a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
    }
    pDenS = prop.den(propSP);
    pDenT = prop.den(propTP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eAAe::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds3>0&&ds2>0){
      //AAEe all in:
      //+ me[12]^2<34>
      //eAAE: 4->1->3->4
      //- me[23]^2<14>
      //34 out:
      //+ me[23]^2<14>
      return + 2.0*e*e*me*s23s.v()*s23s.v()*a14a.v(ds1,ds4)/pDenS/pDenT;
    }
    else if(ds3<0&&ds2<0){
      //AAEe all in:
      //+ me<12>^2[34]
      //eAAE: 4->1->3->4
      //- me<23>^2[14]
      //34 out:
      //- me<23>^2[14]
      return - 2.0*e*e*me*a23a.v()*a23a.v()*s14s.v(ds1,ds4)/pDenS/pDenT;
    }
    else if(ds3>0&&ds2<0){
      //AAEe all in:
      //+ ([13]<24>+[14]<23>)[132>
      //eAAE: 4->1->3->4
      //- ([34]<12>+[13]<24>)[342>
      //34 out:
      //+ ([34]<12>-[13]<24>)[342>
      return + 2.0*e*e*(s34s.v(ds4)*a12a.v(ds1)-s13s.v(ds1)*a24a.v(ds4))*s342a.v()/pDenS/pDenT;
    }
    else if(ds3<0&&ds2>0){
      //AAEe all in:
      //+ (<31>[42]+<41>[32])[231>
      //eAAE: 4->1->3->4
      //- (<34>[12]+<13>[24])[243>
      //34 out:
      //- (<34>[12]-<13>[24])[243>
      return -2.0*e*e*(a34a.v(ds4)*s12s.v(ds1)-a13a.v(ds1)*s24s.v(ds4))*s243a.v()/pDenS/pDenT;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble eAAe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    return amp2/4.0;
  }

  //A+, A+ -> e, E
  ldouble eAAe::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spin 1/2
    return amp2/2.0;
  }



  //  Tests
  int test_eAAe(){
    int n=0;//Number of fails
    std::cout<<"\t* e , A  -> A , e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  me=0.0005, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333;
      eAAe eAAeAmp = eAAe(EE,me);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.715591945406541E-01,2.584715269726006E-01,1.566250707269783E-01,1.135273589580032E-01,9.001256201475140E-02,7.539908999232167E-02,6.557863908111648E-02,5.863400083665116E-02,5.355018255301988E-02,4.973955687417102E-02,4.683835448445317E-02,4.460933452583255E-02,4.289117321484610E-02,4.157034640757527E-02,4.056464696532040E-02,3.981308249485267E-02,3.926944663766538E-02,3.889809448862478E-02,3.867108814911252E-02,3.856622057150888E-02};
      i += eAAeAmp.test_2to2_amp2([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      i += eAAeAmp.test_2to2_amp2_rotations([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      i += eAAeAmp.test_2to2_amp2_boosts([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      i += eAAeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH);
      //std::cout<<"\n#  e , A+ -> A , e\n";
      ldouble dataCHpp[20] = {7.715591945406541E-01,2.584715269726006E-01,1.566250707269783E-01,1.135273589580032E-01,9.001256201475140E-02,7.539908999232167E-02,6.557863908111648E-02,5.863400083665116E-02,5.355018255301988E-02,4.973955687417102E-02,4.683835448445317E-02,4.460933452583255E-02,4.289117321484610E-02,4.157034640757527E-02,4.056464696532040E-02,3.981308249485267E-02,3.926944663766538E-02,3.889809448862478E-02,3.867108814911252E-02,3.856622057150888E-02};
      i += eAAeAmp.test_2to2_amp2([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCHpp);
      i += eAAeAmp.test_2to2_amp2_rotations([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCHpp);
      i += eAAeAmp.test_2to2_amp2_boosts([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCHpp);
      i += eAAeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCHpp);
      //Close to threshold
      //std::cout<<"\n#  me=0.0005, pspatial=0.0001\n";
      pspatial = 0.0001;
      ldouble dataCH2[20] = {3.872099767650941E-02,3.369235237512078E-02,2.963215132430944E-02,2.642745798456125E-02,2.397936285598197E-02,2.220104378555605E-02,2.101612465676533E-02,2.035728184916531E-02,2.016505723083364E-02,2.038684394214673E-02,2.097601724990601E-02,2.189118760822890E-02,2.309555699830748E-02,2.455636282108254E-02,2.624439623212314E-02,2.813358395211684E-02,3.020062435062220E-02,3.242467005748808E-02,3.478705056324546E-02,3.727102927292252E-02};
      i += eAAeAmp.test_2to2_amp2([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      i += eAAeAmp.test_2to2_amp2_rotations([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      i += eAAeAmp.test_2to2_amp2_boosts([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      i += eAAeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAAeAmp.amp2(); }, me,0,0,me,pspatial,dataCH2);
      //std::cout<<"\n#  e , A+ -> A , e\n";
      ldouble dataCH2pp[20] = {3.872099767650941E-02,3.369235237512078E-02,2.963215132430944E-02,2.642745798456125E-02,2.397936285598197E-02,2.220104378555605E-02,2.101612465676533E-02,2.035728184916531E-02,2.016505723083364E-02,2.038684394214673E-02,2.097601724990601E-02,2.189118760822890E-02,2.309555699830748E-02,2.455636282108254E-02,2.624439623212314E-02,2.813358395211684E-02,3.020062435062220E-02,3.242467005748808E-02,3.478705056324546E-02,3.727102927292252E-02};
      i += eAAeAmp.test_2to2_amp2([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCH2pp);
      i += eAAeAmp.test_2to2_amp2_rotations([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCH2pp);
      i += eAAeAmp.test_2to2_amp2_boosts([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCH2pp);
      i += eAAeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eAAeAmp.amp2_Aplus(); }, me,0,0,me,pspatial,dataCH2pp);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}


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

//File:  SPINAS/SM/ueeu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ueeu.h"

namespace spinas {
  //Constructors
  ueeu::ueeu(const ldouble& echarge, const ldouble& massu, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mu(massu), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propA(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(mu)), p2(particle(me)),
    p3(particle(me)), p4(particle(mu)),
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
    preh = e*e*mu*me/(4.0*MW*MW*SW*SW);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLe=-1.0+2.0*SW*SW;
    gRe=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLe-gRe)*mu*me/MZ/MZ;//=-preh!
  }
  void ueeu::set_masses(const ldouble& massu, const ldouble& masse, const ldouble& massh, const ldouble& massW){
    mu=massu;
    me=masse;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(me);
    p3.set_mass(me);
    p4.set_mass(mu);
    preh = e*e*mu*me/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLe-gRe)*mu*me/MZ/MZ;//=-preh
  }
  void ueeu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenUA = propA.denominator(propP);
    pDenUh = proph.denominator(propP);
    pDenUZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ueeu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //uUEe:
    //all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //ueEU: 2<->4
    //+ (+ <13>[24] + <12>[34] + [13]<24> + [12]<34>)
    //34 out:
    //+ (- <13>[24] + <12>[34] - [13]<24> + [12]<34>)
    amplitude = - two*e*e*2.0/3.0*(
				   - a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				   - s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				   )/pDenUA;
    
    //Higgs
    //preh = e*e*mu*me/(4*MW*MW*SW*SW);
    //uUEe:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //ueeu: 2<->4
    //- preh ([14]+<14>) ([23]+<23>)/(u-Mh^2)
    //34 out:
    //- preh ([14]-<14>) ([23]-<23>)/(u-Mh^2)
    amplitude += - preh*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))*(s23s.v(ds2,ds3)-a23a.v(ds2,ds3))/pDenUh;
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLe=-1.0+2.0*SW*SW;
    //gRe=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLu-gRu)*(gLe-gRe)*mu*me/MZ/MZ; // = -preh
    //uUEe:
    //all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ ( gLu gLe [23] <14> + gLe gRu [13] <24> + gLu gRe [24] <13> + gRu gRe [14] <23> )/(s-MZ^2)
    //ueeu: 2<->4
    //+ preZ0 (<14>-[14]) (<23>-[23]) / (u-MZ^2)
    //+ preZ ( gLu gLe [34] <12> + gLe gRu [13] <24> + gLu gRe [24] <13> + gRu gRe [12] <34> )/(u-MZ^2)
    //34 out:
    //+ preZ0 (<14>+[14]) (<23>+[23]) / (u-MZ^2)
    //+ preZ ( gLu gLe [34] <12> - gLe gRu [13] <24> - gLu gRe [24] <13> + gRu gRe [12] <34> )/(u-MZ^2)
    amplitude += 
      + preZ0*(a14a.v(ds1,ds4)+s14s.v(ds1,ds4))*(a23a.v(ds2,ds3)+s23s.v(ds2,ds3))/pDenUZ
      + two*preZ*(
	        gLu*gLe*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      - gLe*gRu*s13s.v(ds1,ds3)*a24a.v(ds2,ds4)
	      - gLu*gRe*s24s.v(ds2,ds4)*a13a.v(ds1,ds3)
	      + gRu*gRe*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenUZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ueeu::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += three*std::pow(std::abs(M),2);// Color factor 3
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3
    return amp2/12.0;
  }

  



  //  Tests
  int test_ueeu(){
    int n=0;//Number of fails
    std::cout<<"\t* u , e  -> u , e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrame\n";
      //std::cout<<"########### mu=0.0042, me=0.0005, pspatial=250\n";
      ldouble mu=0.0042, me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrame.
      ueeu ueeuAmp = ueeu(0.31333,mu,me,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.203337029588701E-02,2.446129945255240E-02,2.732779584667764E-02,3.074310518098513E-02,3.485415109108403E-02,3.986015911289195E-02,4.603661867274324E-02,5.377307936226636E-02,6.363471441051986E-02,7.646638806104436E-02,9.357637497348949E-02,1.170778331197952E-01,1.505641858992534E-01,2.005513653720402E-01,2.798721708492081E-01,4.167560877996676E-01,6.838430979324290E-01,1.319435075871543E+00,3.540716635724392E+00,2.922597308766511E+01};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, me=1.23, pspatial=1.25\n";
      me=1.23;
      pspatial = 1.25;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH2[20] = {9.228238439802364E-03,1.074629127352355E-02,1.261690953536694E-02,1.493873544589525E-02,1.784535573962299E-02,2.152084902616552E-02,2.622391269187930E-02,3.232638366008015E-02,4.037677560622735E-02,5.120910167238656E-02,6.613766308940777E-02,8.732446322743519E-02,1.185175213119366E-01,1.666551574525751E-01,2.457158329256048E-01,3.872485529460365E-01,6.749240937717426E-01,1.393762687021175E+00,4.076166722808565E+00,3.859558856009919E+01};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH2);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH2);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH2);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, me=0.0042, pspatial=0.005\n";
      mu=1.23;
      me=0.0042;
      pspatial = 0.005;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH3[20] = {1.993164495189151E+02,2.367560567318787E+02,2.816961895126527E+02,3.361218468130789E+02,4.027009560304600E+02,4.850827762788198E+02,5.883589567833865E+02,7.197965002865220E+02,8.900416151275052E+02,1.115172981100251E+03,1.420361837859717E+03,1.846747511449453E+03,2.465197807510277E+03,3.406089059823199E+03,4.930484827023270E+03,7.624053648998741E+03,1.303077746137484E+04,2.637871656740895E+04,7.560308903360760E+04,7.013877569590977E+05};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH3);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH3);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH3);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, me=1.23, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      pspatial = 0.3;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH4[20] = {1.228422298435706E+00,1.381287435839678E+00,1.562125117762099E+00,1.778046487751821E+00,2.038554180887490E+00,2.356572054376795E+00,2.750030473682624E+00,3.244377473668823E+00,3.876688497873196E+00,4.702651017617347E+00,5.808971112700570E+00,7.336597351490647E+00,9.527036835181494E+00,1.282223983520084E+01,1.810249890693422E+01,2.733228142154390E+01,4.566317989440661E+01,9.044581708962342E+01,2.538737369831091E+02,2.308650432107055E+03};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH4);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH4);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH4);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, me=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=2.11;
      pspatial = 0.3;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH5[20] = {1.245019009655291E+00,1.398173536910559E+00,1.579326473259488E+00,1.795593580521372E+00,2.056483302578716E+00,2.374926906481467E+00,2.768864361524316E+00,3.263756377959540E+00,3.896695481697555E+00,4.723392731945036E+00,5.830587725812539E+00,7.359278319198388E+00,9.551047184472786E+00,1.284796666601203E+01,1.813054016818550E+01,2.736363332711973E+01,4.569969207612730E+01,9.049155805598335E+01,2.539409111068414E+02,2.308824466487422E+03};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH5);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH5);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH5);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, me=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH6[20] = {2.006053378475435E+07,2.006053425567244E+07,2.006053479910064E+07,2.006053543139952E+07,2.006053617402702E+07,2.006053705561879E+07,2.006053811514705E+07,2.006053940685516E+07,2.006054100822027E+07,2.006054303329402E+07,2.006054565605802E+07,2.006054915350404E+07,2.006055399027343E+07,2.006056099842555E+07,2.006057179897114E+07,2.006058992801319E+07,2.006062444311221E+07,2.006070510776233E+07,2.006098589104002E+07,2.006432729618960E+07};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH6);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH6);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH6);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, me=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH7[20] = {1.294997352847269E+00,1.451064936267191E+00,1.635464838807675E+00,1.855373538242107E+00,2.120375863072413E+00,2.443500536621567E+00,2.842813860427261E+00,3.343943303342997E+00,3.984206152916398E+00,4.819623853520953E+00,5.937377075979554E+00,7.479112314063769E+00,9.687406818285654E+00,1.300593955933138E+01,1.831798902347429E+01,2.759366266833309E+01,4.599663890199915E+01,9.090896379988746E+01,2.546394001295249E+02,2.310928409497080E+03};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH7);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH7);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH7);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, me=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      ueeuAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH8[20] = {2.773579092586076E+07,2.771927952436107E+07,2.770276873709342E+07,2.768626612972908E+07,2.766978125292176E+07,2.765332632716639E+07,2.763691723219014E+07,2.762057497150636E+07,2.760432790153063E+07,2.758821523529590E+07,2.757229276061343E+07,2.755664259779099E+07,2.754139077055937E+07,2.752674101317383E+07,2.751304550327032E+07,2.750097004660498E+07,2.749194298090001E+07,2.748968434824418E+07,2.750787911276266E+07,2.767115590031227E+07};
      i += ueeuAmp.test_2to2_amp2([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH8);
      i += ueeuAmp.test_2to2_amp2_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH8);
      i += ueeuAmp.test_2to2_amp2_boosts([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH8);
      i += ueeuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ueeuAmp.amp2(); }, mu,me,me,mu,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

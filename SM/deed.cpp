
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

//File:  SPINAS/SM/deed.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/deed.h"

namespace spinas {
  //Constructors
  deed::deed(const ldouble& echarge, const ldouble& massd, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), md(massd), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propA(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(md)), p2(particle(me)),
    p3(particle(me)), p4(particle(md)),
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
    preh = e*e*md*me/(4.0*MW*MW*SW*SW);
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    gLe=-1.0+2.0*SW*SW;
    gRe=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLd-gRd)*(gLe-gRe)*md*me/MZ/MZ;//=-preh!
  }
  void deed::set_masses(const ldouble& massd, const ldouble& masse, const ldouble& massh, const ldouble& massW){
    md=massd;
    me=masse;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p2.set_mass(me);
    p3.set_mass(me);
    p4.set_mass(md);
    preh = e*e*md*me/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLd-gRd)*(gLe-gRe)*md*me/MZ/MZ;//=-preh
  }
  void deed::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenUA = propA.den(propP);
    pDenUh = proph.den(propP);
    pDenUZ = propZ.den(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble deed::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //dDEe all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //deEd: 2<->4
    //+ (<13>[24] + <12>[34] + [13]<24> + [12]<34>)
    //34 out:
    //+ (- <13>[24] + <12>[34] - [13]<24> + [12]<34>)
    amplitude = two*e*e*1.0/3.0*(
				 - a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				 - s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				 )/pDenUA;
    
    //Higgs
    //preh = e*e*md*me/(4*MW*MW*SW*SW);
    //dDEe all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //deEd: 2<->4
    //- preh ([14]+<14>) ([23]+<23>)/(u-Mh^2)
    //34 out:
    //- preh ([14]-<14>) ([23]-<23>)/(u-Mh^2)
    amplitude += - preh*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))*(s23s.v(ds2,ds3)-a23a.v(ds2,ds3))/pDenUh;
    
    //Z Boson
    //Defined above:
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //gLe=-1.0+2.0*SW*SW;
    //gRe=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLd-gRd)*(gLe-gRe)*md*me/MZ/MZ; // = -preh
    //dDEe all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ ( gLd gLe [23] <14> + gLe gRd [13] <24> + gLd gRe [24] <13> + gRd gRe [14] <23> )/(s-MZ^2)
    //deED: 2<->4
    //+ preZ0 (<14>-[14]) (<23>-[23]) / (u-MZ^2)
    //+ preZ ( gLd gLe [34] <12> + gLe gRd [13] <24> + gLd gRe [24] <13> + gRd gRe [12] <34> )/(u-MZ^2)
    //34 out:
    //+ preZ0 (<14>+[14]) (<23>+[23]) / (u-MZ^2)
    //+ preZ ( gLd gLe [34] <12> - gLe gRd [13] <24> - gLd gRe [24] <13> + gRd gRe [12] <34> )/(u-MZ^2)
    amplitude += 
      + preZ0*(a14a.v(ds1,ds4)+s14s.v(ds1,ds4))*(a23a.v(ds2,ds3)+s23s.v(ds2,ds3))/pDenUZ
      + two*preZ*(
	        gLd*gLe*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      - gLe*gRd*s13s.v(ds1,ds3)*a24a.v(ds2,ds4)
	      - gLd*gRe*s24s.v(ds2,ds4)*a13a.v(ds1,ds3)
	      + gRd*gRe*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenUZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble deed::amp2(){
    ldouble amp2 = 0, two=2, three = 3;
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
  int test_deed(){
    int n=0;//Number of fails
    std::cout<<"\t* d , e  -> e , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrame\n";
      //std::cout<<"########### md=0.0042, me=0.0005, pspatial=250\n";
      ldouble md=0.0042, me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrame.
      deed deedAmp = deed(0.31333,md,me,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.164453261807580E-02,1.291503274143451E-02,1.441024897907178E-02,1.618564937205634E-02,1.831483160359411E-02,2.089710794996745E-02,2.406906305603918E-02,2.802265229767543E-02,3.303443834137202E-02,3.951454719114548E-02,4.809215116262735E-02,5.977231652783146E-02,7.624148083200348E-02,1.005076524183719E-01,1.383722347323524E-01,2.022499689144707E-01,3.228890399456751E-01,5.955794144069932E-01,1.462702525085183E+00,9.430590219434052E+00};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH);
      //std::cout<<"########### md=0.0042, me=1.23, pspatial=1.25\n";
      me=1.23;
      pspatial = 1.25;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH2[20] = {2.309713958665334E-03,2.689332907311996E-03,3.157096429170806E-03,3.737665661857699E-03,4.464438021792189E-03,5.383434096317623E-03,6.559329443981080E-03,8.085084895395346E-03,1.009783102565485E-02,1.280607409792561E-02,1.653839365131599E-02,2.183529675100201E-02,2.963379774473075E-02,4.166849207291529E-02,6.143402167706753E-02,9.681768814778503E-02,1.687372933070583E-01,3.484481752290549E-01,1.019051958730694E+00,9.648920979148528E+00};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH2);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH2);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH2);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH2);
      //std::cout<<"########### md=1.23, me=0.0042, pspatial=0.005\n";
      md=1.23;
      me=0.0042;
      pspatial = 0.005;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH3[20] = {4.982911245217164E+01,5.918901426348092E+01,7.042404746766475E+01,8.403046180285122E+01,1.006752391185769E+02,1.212706941936168E+02,1.470897393346267E+02,1.799491252276575E+02,2.225104039581496E+02,2.787932454754385E+02,3.550904596944814E+02,4.616868781279563E+02,6.162994521888141E+02,8.515222653267424E+02,1.232621207208166E+03,1.906013412819620E+03,3.257694366098422E+03,6.594679142939546E+03,1.890077226026527E+04,1.753469392455181E+05};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH3);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH3);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH3);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH3);
      //std::cout<<"########### md=1.2, me=1.23, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      pspatial = 0.3;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH4[20] = {3.071087049960479E-01,3.453250875546069E-01,3.905346169547107E-01,4.445150810311884E-01,5.096421410039287E-01,5.891467643010051E-01,6.875115463405598E-01,8.110985011877466E-01,9.691764969367987E-01,1.175667411374483E+00,1.452247778599143E+00,1.834154761509371E+00,2.381765167307171E+00,3.205566615385937E+00,4.525632333609583E+00,6.833079332912893E+00,1.141580610244909E+01,2.261146927004695E+01,6.346845826463134E+01,5.771626771337659E+02};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH4);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH4);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH4);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH4);
      //std::cout<<"########### md=1.2, me=1.23, MW=2.11, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=2.11;
      pspatial = 0.3;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH5[20] = {3.185971915532846E-01,3.570471255593702E-01,4.025124478296030E-01,4.567749285069662E-01,5.222152377806313E-01,6.020707327096595E-01,7.008322903385753E-01,8.248728551276465E-01,9.834760259459647E-01,1.190584021408694E+00,1.467902279121034E+00,1.850710464996784E+00,2.399450555054329E+00,3.224715438481942E+00,4.546760169286276E+00,6.857045203419450E+00,1.144420497792540E+01,2.264780957698951E+01,6.352326442231602E+01,5.773096174540921E+02};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH5);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH5);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH5);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH5);
      //std::cout<<"########### md=1.2, me=1.23, MW=0.006, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH6[20] = {2.006052662008340E+07,2.006052662122529E+07,2.006052663347065E+07,2.006052666052690E+07,2.006052670752465E+07,2.006052678168033E+07,2.006052689333538E+07,2.006052705763249E+07,2.006052729730782E+07,2.006052764752010E+07,2.006052816457577E+07,2.006052894253801E+07,2.006053014690112E+07,2.006053208840088E+07,2.006053540164941E+07,2.006054153782672E+07,2.006055439446001E+07,2.006058741102691E+07,2.006071369881926E+07,2.006237278434942E+07};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH6);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH6);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH6);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH6);
      //std::cout<<"########### md=1.2, me=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH7[20] = {2.942825819628991E-01,3.312737399596715E-01,3.751133223021725E-01,4.275527410941541E-01,4.909344741522223E-01,5.684471503353977E-01,6.645184805024902E-01,7.854380297482468E-01,9.403770160575428E-01,1.143122465395806E+00,1.415159278445013E+00,1.791442802114547E+00,2.331917723283327E+00,3.146373617915388E+00,4.453677998995333E+00,6.742670411658336E+00,1.129636905548730E+01,2.243974179727523E+01,6.317465259496642E+01,5.762582762149498E+02};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH7);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH7);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH7);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH7);
      //std::cout<<"########### md=1.2, me=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      deedAmp.set_masses(md,me,mh,MW);
      ldouble dataCH8[20] = {2.772625027782976E+07,2.770917627770206E+07,2.769203841068243E+07,2.767483251617733E+07,2.765755339204564E+07,2.764019443568270E+07,2.762274712566826E+07,2.760520025474507E+07,2.758753876285757E+07,2.756974190380547E+07,2.755178025485960E+07,2.753361061728480E+07,2.751516684115469E+07,2.749634219005004E+07,2.747695249218554E+07,2.745665023757752E+07,2.743469171339845E+07,2.740914685957194E+07,2.737296358572724E+07,2.726347628536655E+07};
      i += deedAmp.test_2to2_amp2([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH8);
      i += deedAmp.test_2to2_amp2_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH8);
      i += deedAmp.test_2to2_amp2_boosts([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH8);
      i += deedAmp.test_2to2_amp2_boosts_and_rotations([&]() { return deedAmp.amp2(); }, md,me,me,md,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

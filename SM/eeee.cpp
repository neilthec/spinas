
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

//File:  SPINAS/SM/eeee.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeee.h"

namespace spinas {
  //Constructors
  eeee::eeee(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propA(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(me)), p2(particle(me)),
    p3(particle(me)), p4(particle(me)),
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
    preh = e*e*me*me/(4.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*me/MZ/MZ;//=preh!
  }
  void eeee::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(me);
    p4.set_mass(me);
    preh = e*e*me*me/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*me/MZ/MZ;
  }
  void eeee::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propPS[4], propPT[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenSA = propA.denominator(propPS);
    pDenSh = proph.denominator(propPS);
    pDenSZ = propZ.denominator(propPS);
    pDenTA = propA.denominator(propPT);
    pDenTh = proph.denominator(propPT);
    pDenTZ = propZ.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeee::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble one=1, two = 2, four=4, eight=8;
    cdouble amplitude(0,0);
    
    //Photon
    //S Channel
    // all ingoing: -(<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    // 34 outgoing:  (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude += two*e*e*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSA;

    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    // all ingoing:  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    // 34 outgoing:  (<12>[34] + <14>[23] + [12]<34> + [14]<23>)/pDenT
    amplitude += two*e*e*(
			  a12a.v(ds1,ds2)*s34s.v(ds3,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s12s.v(ds1,ds2)*a34a.v(ds3,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenTA;
    
    //Higgs
    //S Channel
    //EE^2 Me Me / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*me*me/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;

    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    //all ingoing:  - preh * ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 outgoing:  - preh * ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*me*me/MZ/MZ; // = preh
    //S Channel
    //all in:
    //+(EE^2 Me Me (gL-gR)^2 (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= + preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  + preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //34 out:
    //+ preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	      gL*gL*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;

    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    //all in:
    //= - preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //  - preZ (- gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) - gR^2 [14] <23>))/(t-MZ^2)
    //34 outgoing:
    //= - preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //  - preZ ( gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) + gR^2 [14] <23>))/(t-MZ^2)
    amplitude += 
      + preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      + two*preZ*(
	      gL*gL*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gL*gR*(s12s.v(ds1,ds2)*a34a.v(ds3,ds4)+s34s.v(ds3,ds4)*a12a.v(ds1,ds2))
	      + gR*gR*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenTZ;
    
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eeee::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    //M = amp(j1,j2,j3,j4);
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    return amp2/4.0;
  }

  



  //  Tests
  int test_eeee(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> e , E       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, pspatial=250\n";
      ldouble me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      eeee eeeeAmp = eeee(0.31333,me,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {5.903751877109695E+01,5.818183371014408E+00,1.814053801631326E+00,7.960388228797162E-01,4.141571289503319E-01,2.395442914138245E-01,1.495040625696446E-01,9.918050201236045E-02,6.938966335398465E-02,5.097889466341380E-02,3.921602240019052E-02,3.150010674998961E-02,2.633196285372663E-02,2.281449756028290E-02,2.039482675016671E-02,1.872428503093355E-02,1.757916040768066E-02,1.681418470833876E-02,1.633440151595338E-02,1.607769993956932E-02};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, pspatial=0.11\n";
      pspatial = 0.005;
      ldouble dataCH2[20] = {5.930325867429968E+01,5.994538060533914E+00,1.975917903811836E+00,9.300317435147173E-01,5.235520774076968E-01,3.293517267036250E-01,2.240075197987907E-01,1.617249691245188E-01,1.225482186803366E-01,9.674135051857594E-02,7.913870959452271E-02,6.681590141956169E-02,5.802970779981610E-02,5.169646802383646E-02,4.711872008180500E-02,4.383428873622961E-02,4.152924550954444E-02,3.998579947778291E-02,3.905008049892969E-02,3.861165497862376E-02};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      //std::cout<<"########### me=0.10, pspatial=0.05\n";
      me=0.10;
      pspatial = 0.05;
      eeeeAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {5.411094113774052E+02,5.707805183850071E+01,1.949015813018070E+01,9.423540985914370E+00,5.397421464654336E+00,3.417828275016637E+00,2.312677534834390E+00,1.640169343290845E+00,1.204621125352970E+00,9.089378666279678E-01,7.006882624484199E-01,5.496388745306172E-01,4.374173765893017E-01,3.523697516264196E-01,2.868344095315822E-01,2.356196743112483E-01,1.951150624854636E-01,1.627526081247944E-01,1.366695553572831E-01,1.154910524057082E-01};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      //std::cout<<"########### me=0.10, MW=0.11, pspatial=0.05\n";
      me=0.10;
      MW=0.11;
      pspatial = 0.05;
      eeeeAmp.set_masses(me,mh,MW);
      ldouble dataCH5[20] = {5.385422921713418E+02,5.623070362998340E+01,1.898707076386503E+01,9.068191059154865E+00,5.124282713819881E+00,3.197100323064331E+00,2.128300545981931E+00,1.482495780942978E+00,1.067398371101905E+00,7.878790165845272E-01,5.927235475012285E-01,4.524925803305621E-01,3.493540123528018E-01,2.720342384622367E-01,2.131473784382096E-01,1.677065374985379E-01,1.322565743957199E-01,1.043495466501740E-01,8.221694269277191E-02,6.455897962647159E-02};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      //std::cout<<"########### me=0.10, MW=0.006, pspatial=0.05\n";
      me=0.10;
      MW=0.006;
      pspatial = 0.05;
      eeeeAmp.set_masses(me,mh,MW);
      ldouble dataCH6[20] = {3.017284744230272E+03,2.721087597498562E+03,2.721372743208326E+03,2.727685959870128E+03,2.732798771461476E+03,2.736644017305740E+03,2.739570758621132E+03,2.741850694534890E+03,2.743667316973994E+03,2.745143699777569E+03,2.746364026274938E+03,2.747387329788103E+03,2.748256035736770E+03,2.749001331011343E+03,2.749646612288040E+03,2.750209753957868E+03,2.750704634960442E+03,2.751142189464394E+03,2.751531144637055E+03,2.751878548355030E+03};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      //std::cout<<"########### me=0.10, MW=0.11, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      MW=0.11;
      mh=0.125;
      pspatial = 0.05;
      eeeeAmp.set_masses(me,mh,MW);
      ldouble dataCH7[20] = {5.485819602404221E+02,5.942942570643509E+01,2.082314815836125E+01,1.032389713029461E+01,6.060137466395549E+00,3.931347085872752E+00,2.724475467483514E+00,1.978616141825107E+00,1.487980148006820E+00,1.149628521571057E+00,9.075180649635580E-01,7.290645209986131E-01,5.943009615694161E-01,4.904545529943030E-01,4.090563185131967E-01,3.443179213511226E-01,2.921788831055286E-01,2.497258325146152E-01,2.148267849955972E-01,1.858941803620416E-01};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      //std::cout<<"########### me=0.10, MW=0.006, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      MW=0.006;
      mh=0.125;
      pspatial = 0.05;
      eeeeAmp.set_masses(me,mh,MW);
      ldouble dataCH8[20] = {9.174753958992236E+03,6.479039441696164E+03,5.944760548136857E+03,5.689361449962438E+03,5.527721387093466E+03,5.410722453750483E+03,5.319483572060074E+03,5.245041642897833E+03,5.182499656965970E+03,5.128896001015356E+03,5.082296339150721E+03,5.041359409758318E+03,5.005110169995782E+03,4.972812535879424E+03,4.943893675194916E+03,4.917896681757509E+03,4.894449708104632E+03,4.873245078184316E+03,4.854024691654146E+03,4.836569532809621E+03};
      i += eeeeAmp.test_2to2_amp2([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      i += eeeeAmp.test_2to2_amp2_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      i += eeeeAmp.test_2to2_amp2_boosts([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      i += eeeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeeeAmp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

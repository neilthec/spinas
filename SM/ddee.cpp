
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

//File:  SPINAS/SM/ddee.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddee.h"

namespace spinas {
  //Constructors
  ddee::ddee(const ldouble& echarge, const ldouble& massd, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), md(massd), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propA(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(md)), p2(particle(md)),
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
    preh = e*e*md*me/(4.0*MW*MW*SW*SW);
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    gLe=-1.0+2.0*SW*SW;
    gRe=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLd-gRd)*(gLe-gRe)*md*me/MZ/MZ;//=-preh!
  }
  void ddee::set_masses(const ldouble& massd, const ldouble& masse, const ldouble& massh, const ldouble& massW){
    md=massd;
    me=masse;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p2.set_mass(md);
    p3.set_mass(me);
    p4.set_mass(me);
    preh = e*e*md*me/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLd-gRd)*(gLe-gRe)*md*me/MZ/MZ;//=-preh
  }
  void ddee::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenSA = propA.den(propP);
    pDenSh = proph.den(propP);
    pDenSZ = propZ.den(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddee::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*e*e*1.0/3.0*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSA;
    
    //Higgs
    //EE^2 Md Me / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*md*me/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //gLe=-1.0+2.0*SW*SW;
    //gRe=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLd-gRd)*(gLe-gRe)*md*me/MZ/MZ; // = -preh
    //all in:
    //
    //+(EE^2 Md Me (gLd-gRd)(gLe-gRe) (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 ( gLd gLe [23] <14> + gLe gRd [13] <24> + gLd gRe [24] <13> + gRd gRe [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= + preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  + preZ ( gLd gLe [23] <14> + gLe gRd [13] <24> + gLd gRe [24] <13> + gRd gRe [14] <23> )/(s-MZ^2)
    //34 out:
    //+ preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ ( gLd gLe [23] <14> + gLe gRd [13] <24> + gLd gRe [24] <13> + gRd gRe [14] <23> )/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	        gLd*gLe*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gLe*gRd*s13s.v(ds1,ds3)*a24a.v(ds2,ds4)
	      + gLd*gRe*s24s.v(ds2,ds4)*a13a.v(ds1,ds3)
	      + gRd*gRe*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ddee::amp2(){
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
    //Average over colors 1/3*1/3 = 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_ddee(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> e , E       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrame\n";
      //std::cout<<"########### md=0.0042, me=0.0005, pspatial=250\n";
      ldouble md=0.0042, me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrame.
      ddee ddeeAmp = ddee(0.31333,md,me,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.803660418913352E-03,3.425157392716273E-03,3.068252462070420E-03,2.732945626975794E-03,2.419236887432394E-03,2.127126243440222E-03,1.856613694999275E-03,1.607699242109556E-03,1.380382884771063E-03,1.174664622983796E-03,9.905444567477567E-04,8.280223860629437E-04,6.870984109293574E-04,5.677725313469977E-04,4.700447473158647E-04,3.939150588359583E-04,3.393834659072786E-04,3.064499685298255E-04,2.951145667035992E-04,3.053772604285994E-04};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH);
      //std::cout<<"########### md=0.0042, me=1.23, pspatial=1.25\n";
      me=1.23;
      pspatial = 1.25;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH2[20] = {7.125821685116313E-04,7.105625170442629E-04,7.087695548967375E-04,7.072032820690553E-04,7.058636985612162E-04,7.047508043732202E-04,7.038645995050673E-04,7.032050839567575E-04,7.027722577282909E-04,7.025661208196673E-04,7.025866732308869E-04,7.028339149619496E-04,7.033078460128553E-04,7.040084663836042E-04,7.049357760741961E-04,7.060897750846313E-04,7.074704634149094E-04,7.090778410650307E-04,7.109119080349951E-04,7.129726643248027E-04};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH2);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH2);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH2);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH2);
      //std::cout<<"########### md=1.23, me=0.0042, pspatial=0.005\n";
      md=1.23;
      me=0.0042;
      pspatial = 0.005;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH3[20] = {7.138806860352774E-04,7.138800783168108E-04,7.138795885608225E-04,7.138792167673125E-04,7.138789629362807E-04,7.138788270677272E-04,7.138788091616519E-04,7.138789092180548E-04,7.138791272369360E-04,7.138794632182954E-04,7.138799171621331E-04,7.138804890684490E-04,7.138811789372431E-04,7.138819867685156E-04,7.138829125622662E-04,7.138839563184951E-04,7.138851180372022E-04,7.138863977183876E-04,7.138877953620512E-04,7.138893109681931E-04};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH3);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH3);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH3);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH3);
      //std::cout<<"########### md=1.2, me=1.23, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      pspatial = 0.3;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH4[20] = {1.046019126663313E-03,1.045979782839618E-03,1.045945132351040E-03,1.045915175197580E-03,1.045889911379236E-03,1.045869340896010E-03,1.045853463747901E-03,1.045842279934909E-03,1.045835789457034E-03,1.045833992314277E-03,1.045836888506636E-03,1.045844478034113E-03,1.045856760896706E-03,1.045873737094417E-03,1.045895406627245E-03,1.045921769495190E-03,1.045952825698253E-03,1.045988575236432E-03,1.046029018109729E-03,1.046074154318142E-03};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH4);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH4);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH4);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH4);
      //std::cout<<"########### md=1.2, me=1.23, MW=2.11, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=2.11;
      pspatial = 0.3;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH5[20] = {9.234808560675727E-03,8.993484184364664E-03,8.754325580360638E-03,8.517332748663653E-03,8.282505689273708E-03,8.049844402190804E-03,7.819348887414940E-03,7.591019144946116E-03,7.364855174784332E-03,7.140856976929588E-03,6.919024551381884E-03,6.699357898141220E-03,6.481857017207596E-03,6.266521908581012E-03,6.053352572261469E-03,5.842349008248965E-03,5.633511216543502E-03,5.426839197145079E-03,5.222332950053696E-03,5.019992475269353E-03};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH5);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH5);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH5);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH5);
      //std::cout<<"########### md=1.2, me=1.23, MW=0.006, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH6[20] = {6.686841424751708E+06,6.686841424848984E+06,6.686841424946273E+06,6.686841425043576E+06,6.686841425140892E+06,6.686841425238223E+06,6.686841425335567E+06,6.686841425432923E+06,6.686841425530294E+06,6.686841425627679E+06,6.686841425725077E+06,6.686841425822489E+06,6.686841425919914E+06,6.686841426017352E+06,6.686841426114805E+06,6.686841426212272E+06,6.686841426309751E+06,6.686841426407244E+06,6.686841426504751E+06,6.686841426602271E+06};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH6);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH6);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH6);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH6);
      //std::cout<<"########### md=1.2, me=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH7[20] = {9.165901283411802E-03,8.931915572926183E-03,8.700095634747601E-03,8.470441468876062E-03,8.242953075311562E-03,8.017630454054101E-03,7.794473605103682E-03,7.573482528460302E-03,7.354657224123962E-03,7.137997692094663E-03,6.923503932372403E-03,6.711175944957184E-03,6.501013729849005E-03,6.293017287047866E-03,6.087186616553767E-03,5.883521718366708E-03,5.682022592486689E-03,5.482689238913711E-03,5.285521657647772E-03,5.090519848688873E-03};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH7);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH7);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH7);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH7);
      //std::cout<<"########### md=1.2, me=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      ddeeAmp.set_masses(md,me,mh,MW);
      ldouble dataCH8[20] = {6.699226258471670E+06,6.699226695210466E+06,6.699227131949277E+06,6.699227568688100E+06,6.699228005426937E+06,6.699228442165788E+06,6.699228878904652E+06,6.699229315643530E+06,6.699229752382422E+06,6.699230189121327E+06,6.699230625860246E+06,6.699231062599178E+06,6.699231499338125E+06,6.699231936077084E+06,6.699232372816057E+06,6.699232809555044E+06,6.699233246294045E+06,6.699233683033058E+06,6.699234119772086E+06,6.699234556511127E+06};
      i += ddeeAmp.test_2to2_amp2([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH8);
      i += ddeeAmp.test_2to2_amp2_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH8);
      i += ddeeAmp.test_2to2_amp2_boosts([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH8);
      i += ddeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddeeAmp.amp2(); }, md,md,me,me,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

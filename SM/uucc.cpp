
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

//File:  SPINAS/SM/uucc.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uucc.h"

namespace spinas {
  //Constructors
  uucc::uucc(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massc, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), mc(massc), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propAG(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(mu)), p2(particle(mu)),
    p3(particle(mc)), p4(particle(mc)),
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
    preh = e*e*mu*mc/(4.0*MW*MW*SW*SW);
    gL=1.0-4.0/3.0*SW*SW;
    gR=-4.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mc/MZ/MZ;//=preh!
  }
  void uucc::set_masses(const ldouble& massu, const ldouble& massc, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mc=massc;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(mc);
    p4.set_mass(mc);
    preh = e*e*mu*mc/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mc/MZ/MZ;//=preh
  }
  void uucc::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenSAG = propAG.denominator(propP);
    pDenSh = proph.denominator(propP);
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon so we can separate the color factor between the gluon^2 diagram
  // And the rest^2.
  // The cross term vanishes due to the trace of the adjoint rep.
  cdouble uucc::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*gs*gs*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSAG;

    return amplitude;
  }
  cdouble uucc::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*e*e*4.0/9.0*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSAG;
    
    //Higgs
    //EE^2 Me Mm / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gL=1.0-4.0/3.0*SW*SW;
    //gR=-4.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mc/MZ/MZ; // = preh
    //all in:
    //-(EE^2 Mu Mc (gL-gR)^2 (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //-(EE^2 (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= - preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  - preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //34 out:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //+ preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	      gL*gL*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uucc::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, nine = 9;
    cdouble M_rest, M_gluon;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest = amp_rest(j1,j2,j3,j4);
	    M_gluon = amp_gluon(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M_rest),2);// Color factor 3*3=9
	    //Cross term color factor 0 (Trace(Ta)*Trace(Ta)=0)
	    amp2 += two*std::pow(std::abs(M_gluon),2);//Color factor C^2*delta^ab*delta^ab = 1/4*8=2
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_uucc(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> c , C       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, mc=1.23, pspatial=250\n";
      ldouble mu=0.0042, mc=1.23, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      uucc uuccAmp = uucc(0.31333,1.238,mu,mc,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.008766180484270E+00,9.132452209670370E-01,8.282480707398938E-01,7.537747298028404E-01,6.898251981558768E-01,6.363994757990029E-01,5.934975627322187E-01,5.611194589555243E-01,5.392651644689196E-01,5.279346792724046E-01,5.271280033659793E-01,5.368451367496438E-01,5.570860794233981E-01,5.878508313872420E-01,6.291393926411758E-01,6.809517631851993E-01,7.432879430193124E-01,8.161479321435154E-01,8.995317305578081E-01,9.934393382621904E-01};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, mc=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {1.046188899933409E+00,1.043194445293522E+00,1.040532716969518E+00,1.038203714961397E+00,1.036207439269159E+00,1.034543889892804E+00,1.033213066832333E+00,1.032214970087744E+00,1.031549599659038E+00,1.031216955546216E+00,1.031217037749277E+00,1.031549846268220E+00,1.032215381103047E+00,1.033213642253757E+00,1.034544629720350E+00,1.036208343502826E+00,1.038204783601185E+00,1.040533950015426E+00,1.043195842745552E+00,1.046190461791560E+00};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH2);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH2);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH2);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, mc=0.0042, pspatial=0.005\n";
      mu=1.23;
      mc=0.0042;
      pspatial = 0.005;
      uuccAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH3[20] = {1.047811265418240E+00,1.047809708970209E+00,1.047808325662585E+00,1.047807115495368E+00,1.047806078468559E+00,1.047805214582157E+00,1.047804523836162E+00,1.047804006230575E+00,1.047803661765395E+00,1.047803490440622E+00,1.047803492256256E+00,1.047803667212298E+00,1.047804015308747E+00,1.047804536545604E+00,1.047805230922868E+00,1.047806098440539E+00,1.047807139098617E+00,1.047808352897103E+00,1.047809739835996E+00,1.047811299915296E+00};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH3);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH3);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH3);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, mc=1.23, pspatial=0.05\n";
      mu=1.2;
      mc=1.23;
      pspatial = 0.3;
      uuccAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH4[20] = {1.535346547445029E+00,1.535284560818502E+00,1.535229462882018E+00,1.535181253635578E+00,1.535139933079180E+00,1.535105501212826E+00,1.535077958036515E+00,1.535057303550247E+00,1.535043537754021E+00,1.535036660647839E+00,1.535036672231700E+00,1.535043572505604E+00,1.535057361469552E+00,1.535078039123542E+00,1.535105605467575E+00,1.535140060501651E+00,1.535181404225771E+00,1.535229636639934E+00,1.535284757744139E+00,1.535346767538388E+00};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH4);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH4);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH4);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=2.11, pspatial=0.05\n";
      mu=1.2;
      mc=1.23;
      MW=2.11;
      pspatial = 0.3;
      uuccAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH5[20] = {1.605106671998344E+00,1.603662983342511E+00,1.602232057454497E+00,1.600813894334302E+00,1.599408493981925E+00,1.598015856397366E+00,1.596635981580626E+00,1.595268869531704E+00,1.593914520250601E+00,1.592572933737316E+00,1.591244109991850E+00,1.589928049014202E+00,1.588624750804373E+00,1.587334215362362E+00,1.586056442688170E+00,1.584791432781796E+00,1.583539185643241E+00,1.582299701272504E+00,1.581072979669585E+00,1.579859020834485E+00};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH5);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH5);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH5);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=0.006, pspatial=0.05\n";
      mu=1.2;
      mc=1.23;
      MW=0.006;
      pspatial = 0.3;
      uuccAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH6[20] = {2.006052580641967E+07,2.006052580676554E+07,2.006052580711833E+07,2.006052580747803E+07,2.006052580784465E+07,2.006052580821818E+07,2.006052580859864E+07,2.006052580898601E+07,2.006052580938030E+07,2.006052580978150E+07,2.006052581018963E+07,2.006052581060467E+07,2.006052581102662E+07,2.006052581145550E+07,2.006052581189129E+07,2.006052581233399E+07,2.006052581278362E+07,2.006052581324016E+07,2.006052581370362E+07,2.006052581417400E+07};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH6);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH6);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH6);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=2.11, Mh=3.125, pspatial=0.05\n";
      mu=1.2;
      mc=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      uuccAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH7[20] = {1.604722490094385E+00,1.603319497443625E+00,1.601929267560684E+00,1.600551800445561E+00,1.599187096098257E+00,1.597835154518772E+00,1.596495975707104E+00,1.595169559663256E+00,1.593855906387225E+00,1.592555015879014E+00,1.591266888138620E+00,1.589991523166045E+00,1.588728920961289E+00,1.587479081524351E+00,1.586242004855231E+00,1.585017690953930E+00,1.583806139820448E+00,1.582607351454784E+00,1.581421325856938E+00,1.580248063026911E+00};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH7);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH7);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH7);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=0.006, Mh=3.125, pspatial=0.05\n";
      mu=1.2;
      mc=1.23;
      MW=0.006;
      pspatial = 0.3;
      uuccAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH8[20] = {2.009767533702833E+07,2.009767717051468E+07,2.009767900400795E+07,2.009768083750813E+07,2.009768267101523E+07,2.009768450452925E+07,2.009768633805018E+07,2.009768817157803E+07,2.009769000511280E+07,2.009769183865449E+07,2.009769367220309E+07,2.009769550575861E+07,2.009769733932105E+07,2.009769917289041E+07,2.009770100646668E+07,2.009770284004987E+07,2.009770467363997E+07,2.009770650723699E+07,2.009770834084094E+07,2.009771017445179E+07};
      i += uuccAmp.test_2to2_amp2([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH8);
      i += uuccAmp.test_2to2_amp2_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH8);
      i += uuccAmp.test_2to2_amp2_boosts([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH8);
      i += uuccAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuccAmp.amp2(); }, mu,mu,mc,mc,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

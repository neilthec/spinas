
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

//File:  SPINAS/SM/ddss.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddss.h"

namespace spinas {
  //Constructors
  ddss::ddss(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& masss, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), md(massd), ms(masss), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propAG(0,0), proph(mh,wh),
    p1(particle(md)), p2(particle(md)),
    p3(particle(ms)), p4(particle(ms)),
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
    preh = e*e*md*ms/(4.0*MW*MW*SW*SW);
    gL=-1.0+2.0/3.0*SW*SW;
    gR=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*md*ms/MZ/MZ;//=preh!
  }
  void ddss::set_masses(const ldouble& massd, const ldouble& masss, const ldouble& massh, const ldouble& massW){
    md=massd;
    ms=masss;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p2.set_mass(md);
    p3.set_mass(ms);
    p4.set_mass(ms);
    preh = e*e*md*ms/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*md*ms/MZ/MZ;//=preh
  }
  void ddss::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble ddss::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
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
  cdouble ddss::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*e*e*1.0/9.0*(
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
    //preZ0 = preZ*(gL-gR)*(gL-gR)*md*ms/MZ/MZ; // = preh
    //all in:
    //+(EE^2 Md Ms (gL-gR)^2 (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= - preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  - preZ 2(gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //34 out:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //+ preZ 2(gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
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
  ldouble ddss::amp2(){
    ldouble amp2 = 0, two=2, three = 3, nine = 9;
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
  int test_ddss(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> s , S       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, ms=1.23, pspatial=250\n";
      ldouble md=0.0042, ms=1.23, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      ddss ddssAmp = ddss(0.31333,1.238,md,ms,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.006147747025814E+00,9.108873218941875E-01,8.261358177880650E-01,7.518932347074471E-01,6.881595726523337E-01,6.349348316227249E-01,5.922190116186207E-01,5.600121126400209E-01,5.383141346869257E-01,5.271250777593350E-01,5.264449418572490E-01,5.362737269806674E-01,5.566114331295905E-01,5.874580603040180E-01,6.288136085039501E-01,6.806780777293868E-01,7.430514679803280E-01,8.159337792567737E-01,8.993250115587240E-01,9.932251648861787E-01};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH);
      //std::cout<<"########### md=0.0042, ms=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {1.042625366564081E+00,1.039641051900750E+00,1.036988330034000E+00,1.034667200963829E+00,1.032677664690239E+00,1.031019721213229E+00,1.029693370532799E+00,1.028698612648949E+00,1.028035447561680E+00,1.027703875270991E+00,1.027703895776881E+00,1.028035509079352E+00,1.028698715178403E+00,1.029693514074035E+00,1.031019905766246E+00,1.032677890255038E+00,1.034667467540410E+00,1.036988637622361E+00,1.039641400500894E+00,1.042625756176006E+00};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH2);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH2);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH2);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH2);
      //std::cout<<"########### md=1.23, ms=0.0042, pspatial=0.005\n";
      md=1.23;
      ms=0.0042;
      pspatial = 0.005;
      ddssAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH3[20] = {1.044241628903207E+00,1.044240076401132E+00,1.044238696449614E+00,1.044237489048655E+00,1.044236454198255E+00,1.044235591898413E+00,1.044234902149128E+00,1.044234384950403E+00,1.044234040302235E+00,1.044233868204626E+00,1.044233868657575E+00,1.044234041661082E+00,1.044234387215147E+00,1.044234905319771E+00,1.044235595974953E+00,1.044236459180694E+00,1.044237494936992E+00,1.044238703243849E+00,1.044240084101264E+00,1.044241637509238E+00};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH3);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH3);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH3);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH3);
      //std::cout<<"########### md=1.2, ms=1.23, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      pspatial = 0.3;
      ddssAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH4[20] = {1.530116063071206E+00,1.530054278963944E+00,1.529999360078574E+00,1.529951306415097E+00,1.529910117973513E+00,1.529875794753822E+00,1.529848336756024E+00,1.529827743980119E+00,1.529814016426107E+00,1.529807154093988E+00,1.529807156983762E+00,1.529814025095429E+00,1.529827758428989E+00,1.529848356984442E+00,1.529875820761788E+00,1.529910149761026E+00,1.529951343982158E+00,1.529999403425183E+00,1.530054328090100E+00,1.530116117976911E+00};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH4);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH4);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH4);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH4);
      //std::cout<<"########### md=1.2, ms=1.23, MW=2.11, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=2.11;
      pspatial = 0.3;
      ddssAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH5[20] = {1.833539646465952E+00,1.830050358337705E+00,1.826577530126277E+00,1.823121161831668E+00,1.819681253453876E+00,1.816257804992903E+00,1.812850816448748E+00,1.809460287821411E+00,1.806086219110892E+00,1.802728610317192E+00,1.799387461440310E+00,1.796062772480246E+00,1.792754543437001E+00,1.789462774310573E+00,1.786187465100964E+00,1.782928615808173E+00,1.779686226432201E+00,1.776460296973046E+00,1.773250827430710E+00,1.770057817805192E+00};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH5);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH5);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH5);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH5);
      //std::cout<<"########### md=1.2, ms=1.23, MW=0.006, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      ddssAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH6[20] = {2.006052580339286E+07,2.006052580355903E+07,2.006052580373210E+07,2.006052580391207E+07,2.006052580409895E+07,2.006052580429274E+07,2.006052580449343E+07,2.006052580470103E+07,2.006052580491554E+07,2.006052580513695E+07,2.006052580536527E+07,2.006052580560050E+07,2.006052580584263E+07,2.006052580609167E+07,2.006052580634762E+07,2.006052580661047E+07,2.006052580688023E+07,2.006052580715689E+07,2.006052580744046E+07,2.006052580773094E+07};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH6);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH6);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH6);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH6);
      //std::cout<<"########### md=1.2, ms=1.23, MW=2.11, Mh=3.125, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ddssAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH7[20] = {1.832711117479200E+00,1.829309298733163E+00,1.825923939903943E+00,1.822555040991543E+00,1.819202601995960E+00,1.815866622917196E+00,1.812547103755250E+00,1.809244044510122E+00,1.805957445181812E+00,1.802687305770321E+00,1.799433626275648E+00,1.796196406697793E+00,1.792975647036756E+00,1.789771347292538E+00,1.786583507465138E+00,1.783412127554556E+00,1.780257207560792E+00,1.777118747483847E+00,1.773996747323720E+00,1.770891207080411E+00};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH7);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH7);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH7);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH7);
      //std::cout<<"########### md=1.2, ms=1.23, MW=0.006, Mh=3.125, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      ddssAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH8[20] = {2.009768278978451E+07,2.009768383827189E+07,2.009768488676617E+07,2.009768593526737E+07,2.009768698377547E+07,2.009768803229048E+07,2.009768908081239E+07,2.009769012934121E+07,2.009769117787694E+07,2.009769222641957E+07,2.009769327496911E+07,2.009769432352556E+07,2.009769537208891E+07,2.009769642065917E+07,2.009769746923633E+07,2.009769851782040E+07,2.009769956641138E+07,2.009770061500927E+07,2.009770166361406E+07,2.009770271222575E+07};
      i += ddssAmp.test_2to2_amp2([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH8);
      i += ddssAmp.test_2to2_amp2_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH8);
      i += ddssAmp.test_2to2_amp2_boosts([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH8);
      i += ddssAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddssAmp.amp2(); }, md,md,ms,ms,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

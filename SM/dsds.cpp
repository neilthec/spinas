
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

//File:  SPINAS/SM/dsds.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dsds.h"

namespace spinas {
  //Constructors
  dsds::dsds(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& masss, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), md(massd), ms(masss), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propAG(0,0), proph(mh,wh), 
    p1(particle(md)), p2(particle(ms)),
    p3(particle(md)), p4(particle(ms)),
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
    preh = e*e*md*ms/(4.0*MW*MW*SW*SW);
    gL=-1.0+2.0/3.0*SW*SW;
    gR=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*md*ms/MZ/MZ;//=preh!
  }
  void dsds::set_masses(const ldouble& massd, const ldouble& masss, const ldouble& massh, const ldouble& massW){
    md=massd;
    ms=masss;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p2.set_mass(ms);
    p3.set_mass(md);
    p4.set_mass(ms);
    preh = e*e*md*ms/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*md*ms/MZ/MZ;//=preh
  }
  void dsds::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
      propP[j] = mom1[j]-mom3[j];
    pDenTAG = propAG.denominator(propP);
    pDenTh = proph.denominator(propP);
    pDenTZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon so we can separate the color factor between the gluon^2 diagram
  // And the rest^2.
  // The cross term vanishes due to the trace of the adjoint rep.
  cdouble dsds::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //dDSs
    //all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //dsDS: 4->2->3->4
    //- (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    //- (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = - two*gs*gs*(
			     a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			     + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			     )/pDenTAG;

    return amplitude;
  }
  cdouble dsds::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //dDSs
    //all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //dsDS: 4->2->3->4
    //- (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    //- (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = - two*e*e*1.0/9.0*(
				   a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				   + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				   )/pDenTAG;
    
    //Higgs
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    //dDSs all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //dsDS: 4->2->3->4
    //- preh ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //- preh ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //Defined above:
    //gL=1.0-4.0/3.0*SW*SW;
    //gR=-4.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*md*ms/MZ/MZ; // = preh
    //dDSs all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //dsDS: 4->2->3->4
    //+ preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //- preZ (gL^2 [34] <12> - gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    //34 out:
    //+ preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //- preZ (gL^2 [34] <12> + gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    amplitude += 
      + preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      - two*preZ*(
	      gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      + gL*gR*(s14s.v(ds1,ds4)*a23a.v(ds2,ds3)+s23s.v(ds2,ds3)*a14a.v(ds1,ds4))
	      + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenTZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dsds::amp2(){
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
  int test_dsds(){
    int n=0;//Number of fails
    std::cout<<"\t* d , s  -> d , s       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, ms=1.23, pspatial=250\n";
      ldouble md=0.0042, ms=1.23, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      dsds dsdsAmp = dsds(0.31333,1.238,md,ms,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.263394469316299E+03,3.456174769736424E+02,1.185194972165383E+02,5.760339690760652E+01,3.320903051621151E+01,2.120115705035533E+01,1.449061879394666E+01,1.040320305362270E+01,7.753620961824529E+00,5.953367215872361E+00,4.684452337730268E+00,3.763404049563808E+00,3.078663632369774E+00,2.559425161525436E+00,2.159112070952319E+00,1.846144351930207E+00,1.598547005225988E+00,1.400680800889099E+00,1.241194341989236E+00,1.111704193961353E+00};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH);
      //std::cout<<"########### md=0.0042, ms=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {4.704145982906134E+03,4.968128649491355E+02,1.698735795899364E+02,8.225953260855248E+01,4.719689639566374E+01,2.994656214654735E+01,2.031049308973026E+01,1.444340686506870E+01,1.064156974552139E+01,8.059320406940742E+00,6.239840306770859E+00,4.919618867471847E+00,3.938465732834692E+00,3.194729329176217E+00,2.621555498254969E+00,2.173623893347821E+00,1.819403895951154E+00,1.536459605616443E+00,1.308509191652674E+00,1.123529941791817E+00};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH2);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH2);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH2);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH2);
      //std::cout<<"########### md=1.23, ms=0.0042, pspatial=0.005\n";
      md=1.23;
      ms=0.0042;
      pspatial = 0.005;
      dsdsAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH3[20] = {8.548734675619254E+07,9.214742378606446E+06,3.215120976504688E+06,1.588232158605500E+06,9.292436479582082E+05,6.009430045574707E+05,4.151448519064252E+05,3.004660655413832E+05,2.250873974924911E+05,1.731181699660438E+05,1.359207918232367E+05,1.084810725600110E+05,8.773106228534071E+04,7.171104091402360E+04,5.912341507515087E+04,4.908245960871349E+04,4.096758828616977E+04,3.433401792424262E+04,2.885657313707010E+04,2.429331600657057E+04};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH3);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH3);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH3);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH3);
      //std::cout<<"########### md=1.2, ms=1.23, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      pspatial = 0.3;
      dsdsAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH4[20] = {2.813855648434557E+05,3.094292465074352E+04,1.102381731064237E+04,5.565569338194255E+03,3.331342176210225E+03,2.206387511326846E+03,1.562813108662947E+03,1.161183536682032E+03,8.942059862933042E+02,7.080141778945174E+02,5.731724723285356E+02,4.725016186571167E+02,3.954336098209148E+02,3.351810835421326E+02,2.872251695148670E+02,2.484641822241379E+02,2.167127131022536E+02,1.903955623794577E+02,1.683545194717316E+02,1.497228562092099E+02};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH4);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH4);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH4);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH4);
      //std::cout<<"########### md=1.2, ms=1.23, MW=2.11, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=2.11;
      pspatial = 0.3;
      dsdsAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH5[20] = {2.813858891322124E+05,3.094303580662996E+04,1.102388582374693E+04,5.565619564402606E+03,3.331382234272388E+03,2.206421091681502E+03,1.562842198550763E+03,1.161209328516793E+03,8.942292516971066E+02,7.080354448723477E+02,5.731921180740187E+02,4.725199220079515E+02,3.954507827321481E+02,3.351972909026748E+02,2.872405421307369E+02,2.484788256135728E+02,2.167267136256432E+02,1.904089916178753E+02,1.683674374220859E+02,1.497353136601338E+02};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH5);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH5);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH5);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH5);
      //std::cout<<"########### md=1.2, ms=1.23, MW=0.006, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      dsdsAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH6[20] = {2.034353488432129E+07,2.009163620925182E+07,2.007160539920136E+07,2.006611775667979E+07,2.006387200541231E+07,2.006274151009241E+07,2.006209491979598E+07,2.006169150705300E+07,2.006142341212500E+07,2.006123649016412E+07,2.006110115662624E+07,2.006100014730172E+07,2.006092284312193E+07,2.006086242470526E+07,2.006081435238220E+07,2.006077551054863E+07,2.006074370424236E+07,2.006071735156099E+07,2.006069528950123E+07,2.006067664785032E+07};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH6);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH6);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH6);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH6);
      //std::cout<<"########### md=1.2, ms=1.23, MW=2.11, Mh=3.125, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      dsdsAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH7[20] = {2.813855381779339E+05,3.094291910259524E+04,1.102381597477952E+04,5.565569799717885E+03,3.331343630240893E+03,2.206389592080125E+03,1.562815619236102E+03,1.161186358982933E+03,8.942090439523431E+02,7.080174186965831E+02,5.731758590107266E+02,4.725051237074583E+02,3.954372123575985E+02,3.351847673486016E+02,2.872289217524657E+02,2.484679925600429E+02,2.167165731016743E+02,1.903994650714997E+02,1.683584590326001E+02,1.497268277266475E+02};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH7);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH7);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH7);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH7);
      //std::cout<<"########### md=1.2, ms=1.23, MW=0.006, Mh=3.125, pspatial=0.05\n";
      md=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      dsdsAmp.set_masses(md,ms,mh,MW);
      ldouble dataCH8[20] = {2.757289869068480E+07,2.741319669700015E+07,2.742571775094280E+07,2.744420089847230E+07,2.746303461224242E+07,2.748164528499917E+07,2.749999919340371E+07,2.751813463659767E+07,2.753609155793227E+07,2.755390144559402E+07,2.757158777829328E+07,2.758916794406679E+07,2.760665492850181E+07,2.762405856045032E+07,2.764138639050727E+07,2.765864430547621E+07,2.767583696052118E+07,2.769296808641854E+07,2.771004071084001E+07,2.772705731986801E+07};
      i += dsdsAmp.test_2to2_amp2([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH8);
      i += dsdsAmp.test_2to2_amp2_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH8);
      i += dsdsAmp.test_2to2_amp2_boosts([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH8);
      i += dsdsAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dsdsAmp.amp2(); }, md,ms,md,ms,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

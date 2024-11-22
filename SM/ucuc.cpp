
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

//File:  SPINAS/SM/ucuc.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ucuc.h"

namespace spinas {
  //Constructors
  ucuc::ucuc(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massc, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), mc(massc), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propAG(0,0), proph(mh,wh),
    p1(particle(mu)), p2(particle(mc)),
    p3(particle(mu)), p4(particle(mc)),
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
    //std::cout<<"preZ0="<<preZ0<<std::endl;
  }
  void ucuc::set_masses(const ldouble& massu, const ldouble& massc, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mc=massc;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mc);
    p3.set_mass(mu);
    p4.set_mass(mc);
    preh = e*e*mu*mc/(4.0*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mc/MZ/MZ;//=preh
    //std::cout<<"preZ0="<<preZ0<<std::endl;
  }
  void ucuc::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble ucuc::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);

    //Gluon
    //uUCc
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //ucUC: 4->2->3->4
    //all in:
    // - (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //34 out:
    // - (<12>[34] + <14>[23] + [12]<34> + [14]<23>)
    amplitude = - two*gs*gs*(
			     + a12a.v(ds1,ds2)*s34s.v(ds3,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) 
			     + s12s.v(ds1,ds2)*a34a.v(ds3,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			     )/pDenTAG;
    
    return amplitude;
  }
  cdouble ucuc::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //uUCc:
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //ucUC: 4->2->3->4
    //all in:
    // - (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //34 out:
    // - (<12>[34] + <14>[23] + [12]<34> + [14]<23>)  
    amplitude += - two*e*e*4.0/9.0*(
				   + a12a.v(ds1,ds2)*s34s.v(ds3,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) 
				   + s12s.v(ds1,ds2)*a34a.v(ds3,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
				   )/pDenTAG;
    
    //Higgs
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    //uUCc:
    //all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //ucUC: 4->2->3->4
    //- preh ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //- preh ([13]-<13>) ([24]-<24>)/(t-Mh^2)    
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //Defined above:
    //gL=1.0-4.0/3.0*SW*SW;
    //gR=-4.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mc/MZ/MZ; // = preh
    //uUCc:
    //all in:
    //= - preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  - 2 preZ (gL^2 [23] <14> + gLgR( [13] <24> + [24] <13> ) + gR^2 [14] <23>)/(s-MZ^2)
    //ucUC: 4->2->3->4
    //  + preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //  - 2 preZ (gL^2 [34] <12> - gLgR( [14] <23> + [23] <14> ) + gR^2 [12] <34>)/(t-MZ^2)
    //34 out:
    //  + preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //  - 2 preZ (gL^2 <12> [34] + gLgR( <23> [14] + <14> [23] ) + gR^2 <34> [12])/(t-MZ^2)
    amplitude += 
      + preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      - two*preZ*(
		  + gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
		  + gL*gR*(s14s.v(ds1,ds4)*a23a.v(ds2,ds3)+s23s.v(ds2,ds3)*a14a.v(ds1,ds4))
		  + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
		  )/pDenTZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ucuc::amp2(){
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
  int test_ucuc(){
    int n=0;//Number of fails
    std::cout<<"\t* u , c  -> u , c       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, mc=1.23, pspatial=250\n";
      ldouble mu=0.0042, mc=1.23, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      ucuc ucucAmp = ucuc(0.31333,1.238,mu,mc,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.273713217753536E+03,3.465470109943232E+02,1.188148731739309E+02,5.774233037087929E+01,3.328836449041349E+01,2.125200551357041E+01,1.452580167084529E+01,1.042890692179141E+01,7.773177922598976E+00,5.968723239199983E+00,4.696816775890225E+00,3.773566425512172E+00,3.087160422402236E+00,2.566633177835786E+00,2.165303464925449E+00,1.851520417225207E+00,1.603259716091448E+00,1.404846969332991E+00,1.244905178665090E+00,1.115031920365472E+00};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, mc=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {4.720227486078775E+03,4.985112711692337E+02,1.704543165369626E+02,8.254075306387703E+01,4.735825180234669E+01,3.004894542441789E+01,2.037993429925762E+01,1.449279063112023E+01,1.067795632522521E+01,8.086879085574921E+00,6.261178691755747E+00,4.936443729416020E+00,3.951936215423402E+00,3.205657069626309E+00,2.630523585699770E+00,2.181060484336015E+00,1.825629356049793E+00,1.541717600660258E+00,1.312987726092752E+00,1.127375918959854E+00};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH2);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH2);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH2);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, mc=0.0042, pspatial=0.005\n";
      mu=1.23;
      mc=0.0042;
      pspatial = 0.005;
      ucucAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH3[20] = {8.577959165492864E+07,9.246243665704541E+06,3.226112108408187E+06,1.593661649214682E+06,9.324203369788296E+05,6.029973732355579E+05,4.165640556814824E+05,3.014932312946217E+05,2.258568756223545E+05,1.737099873985653E+05,1.363854472320981E+05,1.088519232330419E+05,8.803097749384740E+04,7.195619047939795E+04,5.932553289864313E+04,4.925025167376444E+04,4.110763905570954E+04,3.445139133657040E+04,2.895522149406945E+04,2.437636452722841E+04};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH3);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH3);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH3);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, mc=1.23, pspatial=0.3\n";
      mu=1.2;
      mc=1.23;
      pspatial = 0.3;
      ucucAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH4[20] = {2.823475025459295E+05,3.104870538229006E+04,1.106150307264030E+04,5.584595666772516E+03,3.342730629726911E+03,2.213930221684462E+03,1.568155710765713E+03,1.165153137308310E+03,8.972629036278347E+02,7.104345841588971E+02,5.751319117751810E+02,4.741169069615420E+02,3.967854351207948E+02,3.363269308588620E+02,2.882070757660276E+02,2.493135809777164E+02,2.174535669393127E+02,1.910464489395429E+02,1.689300569623107E+02,1.502346998606402E+02};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH4);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH4);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH4);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      mc=1.23;
      MW=2.11;
      pspatial = 0.3;
      ucucAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH5[20] = {2.823479294491429E+05,3.104885342213059E+04,1.106159530238747E+04,5.584663951500585E+03,3.342785591372488E+03,2.213976688598125E+03,1.568196282993216E+03,1.165189375017679E+03,8.972958164617492E+02,7.104648630332028E+02,5.751600502971939E+02,4.741432700578049E+02,3.968103002137543E+02,3.363505137929421E+02,2.882295477936968E+02,2.493350802661939E+02,2.174742065904084E+02,1.910663227181287E+02,1.689492434765688E+02,1.502532656704646E+02};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH5);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH5);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH5);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      mc=1.23;
      MW=0.006;
      pspatial = 0.3;
      ucucAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH6[20] = {2.034427658753384E+07,2.009170632143851E+07,2.007162627916218E+07,2.006612619998469E+07,2.006387578408141E+07,2.006274315489204E+07,2.006209546746721E+07,2.006169144782643E+07,2.006142300164315E+07,2.006123587145220E+07,2.006110041398631E+07,2.006099933241037E+07,2.006092198860132E+07,2.006086155159399E+07,2.006081347444411E+07,2.006077463685560E+07,2.006074284077362E+07,2.006071650222107E+07,2.006069445678326E+07,2.006067583327700E+07};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH6);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH6);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH6);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      mc=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ucucAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH7[20] = {2.823465273837498E+05,3.104838832979886E+04,1.106131760262496E+04,5.584466561380332E+03,3.342632816585941E+03,2.213852304460002E+03,1.568091553134922E+03,1.165099057618452E+03,8.972165197555933E+02,7.103942663384938E+02,5.750964959845643E+02,4.740855329232961E+02,3.967574491027076E+02,3.363018244742997E+02,2.881844458776470E+02,2.492931025354782E+02,2.174349740494148E+02,1.910295213422426E+02,1.689146101644029E+02,1.502205777980797E+02};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH7);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH7);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH7);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, mc=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      mc=1.23;
      MW=0.006;
      pspatial = 0.3;
      ucucAmp.set_masses(mu,mc,mh,MW);
      ldouble dataCH8[20] = {2.748952704412622E+07,2.738540909632766E+07,2.740910649022921E+07,2.743238537089864E+07,2.745388479897373E+07,2.747419229429995E+07,2.749372108661034E+07,2.751271815871458E+07,2.753133397188685E+07,2.754966401461401E+07,2.756777139253312E+07,2.758569934129716E+07,2.760347842506187E+07,2.762113084405845E+07,2.763867310824420E+07,2.765611775218356E+07,2.767347446754476E+07,2.769075087062362E+07,2.770795303464790E+07,2.772508586663137E+07};
      i += ucucAmp.test_2to2_amp2([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH8);
      i += ucucAmp.test_2to2_amp2_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH8);
      i += ucucAmp.test_2to2_amp2_boosts([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH8);
      i += ucucAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ucucAmp.amp2(); }, mu,mc,mu,mc,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

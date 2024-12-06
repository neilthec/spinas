
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

//File:  SPINAS/SM/AdZd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AdZd.h"

namespace spinas {

  AdZd::AdZd(const ldouble& echarge, const ldouble& massd, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qd(-1.0/3.0), md(massd), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(massd,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(md);
    p3=particle(MZ);
    p4=particle(md);
    //Spinor Products
    s12s= sproduct(SQUARE,&p1,&p2,2);
    s13s= sproduct(SQUARE,&p1,&p3,2);
    s14s= sproduct(SQUARE,&p1,&p4,2);
    a23a= sproduct(ANGLE,&p2,&p3,2);
    a34a= sproduct(ANGLE,&p3,&p4,2);
    s123a= sproduct(SQUARE,&p1,&p2,&p3,2);
    s132a= sproduct(SQUARE,&p1,&p3,&p2,2);
    s134a= sproduct(SQUARE,&p1,&p3,&p4,2);
    s143a= sproduct(SQUARE,&p1,&p4,&p3,2);
    //Spinor Products
    a12a= sproduct(ANGLE,&p1,&p2,2);
    a13a= sproduct(ANGLE,&p1,&p3,2);
    a14a= sproduct(ANGLE,&p1,&p4,2);
    s23s= sproduct(SQUARE,&p2,&p3,2);
    s34s= sproduct(SQUARE,&p3,&p4,2);
    a123s= sproduct(ANGLE,&p1,&p2,&p3,2);
    a132s= sproduct(ANGLE,&p1,&p3,&p2,2);
    a134s= sproduct(ANGLE,&p1,&p3,&p4,2);
    a143s= sproduct(ANGLE,&p1,&p4,&p3,2);
    //Couplings
    preTU = 2.0*e*e*Qd/(2.0*MW*SW);
    gLd=-2.0*Qd*SW*SW-1.0;
    gRd=-2.0*Qd*SW*SW;
  }
  void AdZd::set_masses(const ldouble& massd, const ldouble& massW){
    md=massd;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p2.set_mass(md);
    p3.set_mass(MZ);
    p4.set_mass(md);
    prope.set_mass(md);
    //Couplings
    preTU = 2.0*e*e*Qd/(2.0*MW*SW);
  }
  void AdZd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Spinor Products
    s12s.update();
    s13s.update();
    s14s.update();
    a23a.update();
    a34a.update();
    s123a.update();
    s132a.update();
    s134a.update();
    s143a.update();
    //Spinor Products
    a12a.update();
    a13a.update();
    a14a.update();
    s23s.update();
    s34s.update();
    a123s.update();
    a132s.update();
    a134s.update();
    a143s.update();
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=prope.denominator(propSP);
    pDenU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AdZd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      if(ds1>0){
	//preTU = e*e/(MW*SW);
	//AdZD: all in:
	// preTU (gLd<23>(MZ[14][123>-md[13][134>)+gRd<34>(-md[13][132>+MZ[12][143>))/((t-md^2)(u-md^2))
	amplitude += normFactor*preTU*(
				       +gLd*a23a.v(ds2,ds3a)*(MZ*s14s.v(ds4)*s123a.v(ds3b)+md*s13s.v(ds3b)*s134a.v(ds4))
				       +gRd*a34a.v(ds3a,ds4)*(md*s13s.v(ds3b)*s132a.v(ds2)+MZ*s12s.v(ds2)*s143a.v(ds3b))
				       )/pDenU/pDenS;

	
      }
      else if(ds1<0){
	//preTU = e*e/(MW*SW);
	
	amplitude += normFactor*preTU*(
				       +gRd*s23s.v(ds2,ds3a)*(MZ*a14a.v(ds4)*a123s.v(ds3b)+md*a13a.v(ds3b)*a134s.v(ds4))
				       +gLd*s34s.v(ds3a,ds4)*(md*a13a.v(ds3b)*a132s.v(ds2)+MZ*a12a.v(ds2)*a143s.v(ds3b))
				       )/pDenU/pDenS;
      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AdZd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3
    return amp2/12.0;
  }
  //set_momenta(...) must be called before amp2_Aplus().
  ldouble AdZd::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3
    return amp2/6.0;
  }  
  //set_momenta(...) must be called before amp2_Aminus().
  ldouble AdZd::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3
    return amp2/6.0;
  }



  //  Tests
  int test_AdZd(){
    int n=0;//Number of fails
    std::cout<<"\t* A , d  -> Z , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, MW=80.385, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      AdZd AdZdAmp = AdZd(EE,md,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {2.313423400145495E-03,2.209862244489111E-03,2.106739668551023E-03,2.004135414098921E-03,1.902149801421387E-03,1.800910825369659E-03,1.700584405195763E-03,1.601389554313011E-03,1.503621464714256E-03,1.407687782486475E-03,1.314167794177035E-03,1.223913391772175E-03,1.138230809659905E-03,1.059230118793286E-03,9.905559872240809E-04,9.390927880471578E-04,9.195927496607935E-04,9.704116673194628E-04,1.232187467204270E-03,2.970661823129956E-03};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      ldouble dataCHm[20] = {2.272729200722160E-03,2.388484147604428E-03,2.517849177291446E-03,2.663298850291751E-03,2.827946323373748E-03,3.015763555172437E-03,3.231899380954589E-03,3.483150263271069E-03,3.778676651645045E-03,4.131128690652945E-03,4.558482901529678E-03,5.087175344905365E-03,5.757741313728409E-03,6.635660896571393E-03,7.834036064523150E-03,9.566503873873207E-03,1.229085801600982E-02,1.719736208576740E-02,2.865031591083022E-02,8.592841727897131E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      ldouble dataCH[20] = {2.293076300433827E-03,2.299173196046769E-03,2.312294422921234E-03,2.333717132195336E-03,2.365048062397567E-03,2.408337190271048E-03,2.466241893075176E-03,2.542269908792040E-03,2.641149058179651E-03,2.769408236569710E-03,2.936325347853356E-03,3.155544368338770E-03,3.447986061694157E-03,3.847445507682339E-03,4.412296025873615E-03,5.252798330960182E-03,6.605225382835308E-03,9.083886876543433E-03,1.494125168901725E-02,4.444953955105063E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {2.577643933090830E-03,2.487166286815016E-03,2.397355300329849E-03,2.308332184510282E-03,2.220249430458145E-03,2.133301595787670E-03,2.047740884815883E-03,1.963900207244711E-03,1.882228267460399E-03,1.803344704862181E-03,1.728130059663656E-03,1.657879243969718E-03,1.594577789723782E-03,1.541433094781124E-03,1.503985210209492E-03,1.492698677499880E-03,1.529997515960365E-03,1.674184183466987E-03,2.139034414329302E-03,4.848531646547875E-03};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {2.053706527770364E-03,2.157585208492770E-03,2.273677355685376E-03,2.404203599609189E-03,2.551957636397289E-03,2.720503836975601E-03,2.914462682170360E-03,3.139933206666658E-03,3.405135848509673E-03,3.721422641176393E-03,4.104925421689497E-03,4.579367479746383E-03,5.181124526137363E-03,5.968957324285480E-03,7.044361737921088E-03,8.599052175373505E-03,1.104384522298739E-02,1.544686399294268E-02,2.572455980761538E-02,7.712499163302601E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      ldouble dataCH2[20] = {2.315675230430597E-03,2.322375747653893E-03,2.335516328007612E-03,2.356267892059736E-03,2.386103533427717E-03,2.426902716381635E-03,2.481101783493121E-03,2.551916706955684E-03,2.643682057985036E-03,2.762383673019287E-03,2.916527740676577E-03,3.118623361858051E-03,3.387851157930573E-03,3.755195209533302E-03,4.274173474065291E-03,5.045875426436692E-03,6.286921369473879E-03,8.560524088204836E-03,1.393179711097234E-02,4.098676163978694E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, MW=80.385, pspatial=95\n";
      md = 125;
      pspatial = 250;
      AdZdAmp.set_masses(md,MW);
      ldouble dataCH4p[20] = {2.328373000142898E-03,2.269918317584830E-03,2.231394838289576E-03,2.216305776879420E-03,2.229024933860817E-03,2.275084530434844E-03,2.361585096542566E-03,2.497792256037230E-03,2.696027767546154E-03,2.973038774769789E-03,3.352173136083131E-03,3.866972590135187E-03,4.567388383650272E-03,5.531149184282125E-03,6.886027813742648E-03,8.857404312253611E-03,1.188206071182846E-02,1.692628118844484E-02,2.660472331134024E-02,5.096369826938851E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {2.288443200258409E-03,2.431687800144111E-03,2.602980376441149E-03,2.806892407212238E-03,3.049044360491749E-03,3.336424312904258E-03,3.677829934921720E-03,4.084493245898360E-03,4.570982517571952E-03,5.156535565219359E-03,5.867084627910595E-03,6.738427860359455E-03,7.821376325497230E-03,9.190458331189704E-03,1.095936328777115E-02,1.330991052404758E-02,1.654987934898841E-02,2.123542462893869E-02,2.842642160585934E-02,3.977591391085477E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      ldouble dataCH4[20] = {2.308408100200653E-03,2.350803058864471E-03,2.417187607365362E-03,2.511599092045829E-03,2.639034647176283E-03,2.805754421669551E-03,3.019707515732143E-03,3.291142750967795E-03,3.633505142559053E-03,4.064787169994575E-03,4.609628881996863E-03,5.302700225247320E-03,6.194382354573751E-03,7.360803757735915E-03,8.922695550756896E-03,1.108365741815059E-02,1.421597003040843E-02,1.908085290869176E-02,2.751557245859979E-02,4.536980609012164E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      //std::cout<<"\n# md=125, MW=80.385, pspatial=125.1\n";
      md = 125;
      pspatial = 125.1;
      AdZdAmp.set_masses(md,MW);
      ldouble dataCH3p[20] = {2.610423875070498E-03,2.626667254145967E-03,2.662541413450876E-03,2.720934162542806E-03,2.805325146179231E-03,2.919944685754139E-03,3.069986369597693E-03,3.261895754969735E-03,3.503768789419749E-03,3.805911556116338E-03,4.181642498593244E-03,4.648468209931755E-03,5.229850986724589E-03,5.957943971465550E-03,6.877966993209850E-03,8.055484414776876E-03,9.589076978796207E-03,1.163365078612175E-02,1.444628724821790E-02,1.848432944734575E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {2.156199651079690E-03,2.300768385439339E-03,2.462333679697385E-03,2.642770416141756E-03,2.844222853903117E-03,3.069149110863212E-03,3.320372152095128E-03,3.601136745394292E-03,3.915170074412586E-03,4.266740036484538E-03,4.660697801150014E-03,5.102476057597153E-03,5.597983233975173E-03,6.153268584653158E-03,6.773691638316147E-03,7.462011664862122E-03,8.214062824709366E-03,9.008800553373013E-03,9.784422522576362E-03,1.037715850063497E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2_Aminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      ldouble dataCH3[20] = {2.383311763075094E-03,2.463717819792653E-03,2.562437546574131E-03,2.681852289342281E-03,2.824774000041174E-03,2.994546898308676E-03,3.195179260846411E-03,3.431516250182013E-03,3.709469431916168E-03,4.036325796300438E-03,4.421170149871629E-03,4.875472133764454E-03,5.413917110349881E-03,6.055606278059354E-03,6.825829315762998E-03,7.758748039819499E-03,8.901569901752787E-03,1.032122566974738E-02,1.211535488539713E-02,1.443074397399036E-02};
      i += AdZdAmp.test_2to2_amp2([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      i += AdZdAmp.test_2to2_amp2_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      i += AdZdAmp.test_2to2_amp2_boosts([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      i += AdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

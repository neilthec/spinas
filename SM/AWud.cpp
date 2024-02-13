
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2024 Neil Christensen

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

//File:  SPINAS/SM/AWud.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AWud.h"

namespace spinas {

  AWud::AWud(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), mu(massu), md(massd), MW(massW), SW(sinW), WW(widthW) {
    constexpr ldouble sqrt2 = std::sqrt(2);
    propW = propagator(MW,WW);
    propu = propagator(mu,0);
    propd = propagator(md,0);
    p1=particle(0);
    p2=particle(MW);
    p3=particle(mu);
    p4=particle(md);
    //Spinor Products
    s12s= sproduct(SQUARE,&p1,&p2);
    s13s= sproduct(SQUARE,&p1,&p3);
    a24a= sproduct(ANGLE,&p2,&p4);
    s142a= sproduct(SQUARE,&p1,&p4,&p2);
    s143a= sproduct(SQUARE,&p1,&p4,&p3);
    //Spinor Products
    a12a= sproduct(ANGLE,&p1,&p2);
    a14a= sproduct(ANGLE,&p1,&p4);
    s23s= sproduct(SQUARE,&p2,&p3);
    s231a= sproduct(SQUARE,&p2,&p3,&p1);
    s431a= sproduct(SQUARE,&p4,&p3,&p1);
    //prefactor
    pre = sqrt2*e*e/(MW*SW);
  }
  void AWud::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    p2.set_mass(MW);
    p3.set_mass(mu);
    p4.set_mass(md);
    propW.set_mass(MW);
    propu.set_mass(mu);
    propd.set_mass(md);
    pre = sqrt2*e*e/(MW*SW);
  }
  void AWud::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1, two=2, three=3;
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Spinor Products
    s12s.update();
    s13s.update();
    a24a.update();
    s142a.update();
    s143a.update();
    //Spinor Products
    a12a.update();
    a14a.update();
    s23s.update();
    s231a.update();
    s431a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);
    pDenT=propu.denominator(propTP);
    pDenU=propd.denominator(propUP);
    prop = (
	     +one/pDenS*(one/pDenU-one/pDenT)
	     + two/three/pDenT*(one/pDenU-one/pDenS)
	     - one/three/pDenU*(one/pDenT-one/pDenS)
	     )/three;
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AWud::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);

      if(ds1>0){
      
	//STU Diagram
	//pre = sqrt(2)*e*e/(MW*SW);
	//AW+Ud all in:
	//S W: +    pre <24>(-MW[13][142>+[12](Md^2 [13]-Mu[143>))/(s-MW^2)*(1/(u-Md^2) or -1/(t-Mu^2))
	//T u: +2/3 pre <24>(-MW[13][142>+[12](Md^2 [13]-Mu[143>))/(t-Mu^2)*(1/(u-Md^2) or -1/(s-MW^2))
	//U d: -1/3 pre <24>(-MW[13][142>+[12](Md^2 [13]-Mu[143>))/(u-Md^2)*(1/(t-Mu^2) or -1/(s-MW^2))
	amplitude += normFactor*pre*a24a.v(ds2a,ds4)*(-MW*s13s.v(ds3)*s142a.v(ds2b)+s12s.v(ds2b)*(-md*md*s13s.v(ds3)+mu*s143a.v(ds3)))*prop;
	
	
      }
      else if(ds1<0){

	//STU Diagram
	amplitude += normFactor*pre*s23s.v(ds2a,ds3)*(MW*a14a.v(ds4)*s231a.v(ds2b)+a12a.v(ds2b)*(mu*mu*a14a.v(ds4)-md*s431a.v(ds4)))*prop;
	
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AWud::amp2(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;


    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    //Color factor 3
	    amp2 += three*std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AWud::amp2_Aplus(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(1,j2,j3,j4);
	  //Color factor 3
	  amp2 += three*std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/three;
  }

  



  



  //  Tests
  int test_AWud(){
    int n=0;//Number of fails
    std::cout<<"\t* A , W+ -> u , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      AWud AWudAmp = AWud(EE,mu,md,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.267787906037086E-01,1.976271350514920E-01,9.585465005370086E-02,5.480162319857441E-02,3.373470997959619E-02,2.155210451461579E-02,1.399233231550287E-02,9.076958988570118E-03,5.771559118319317E-03,3.494400200855540E-03,1.909039051862807E-03,8.303907899301059E-04,1.852440618534192E-04,8.323814334116979E-06,4.780684090621213E-04,2.029606965693192E-03,5.683247893134909E-03,1.417968191899982E-02,3.767500022159385E-02,1.678124737591666E-01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      ldouble dataCHp[20] = {3.653932618889176E-03,4.347770446387526E-03,5.193176500618757E-03,5.758485702910375E-03,6.008974727631038E-03,5.949740944672197E-03,5.599636373907372E-03,4.987953633624118E-03,4.156146946766907E-03,3.162486138863000E-03,2.090368261547287E-03,1.062810367641339E-03,2.686931490118320E-04,1.331649460608096E-05,8.241597625723272E-04,3.697691403144457E-03,1.076930718501295E-02,2.759114233033009E-02,7.452115510750416E-02,3.347812553696825E-01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {1.046160624354637E+00,2.930356438121826E-01,1.465749786447161E-01,8.645976640834536E-02,5.488619836437452E-02,3.609994166743666E-02,2.405362993489026E-02,1.593970965539732E-02,1.029025006501710E-02,6.279727131859177E-03,3.430701599095401E-03,1.480523508131968E-03,3.253002074826283E-04,1.430911656746394E-05,8.007682555374149E-04,3.302159948369257E-03,8.966381788424830E-03,2.168268867426649E-02,5.586336956384493E-02,2.415573967764732E-01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      ldouble dataCH2p[20] = {7.628632224670011E-02,3.380420944191961E-02,2.497739363942367E-02,2.068180293044667E-02,1.769970302708733E-02,1.517363563317472E-02,1.279712234408074E-02,1.045542662916672E-02,8.123616174383310E-03,5.834555133430665E-03,3.673905166722517E-03,1.792250843786030E-03,4.372243507884834E-04,2.100543125247846E-05,1.264955189685594E-03,5.539439210055444E-03,1.578793993263200E-02,3.967049637748794E-02,1.052824124907572E-01,4.655003476273998E-01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      AWudAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {4.955386008267213E-01,2.083549637769480E-01,1.158776111717062E-01,7.179672718189153E-02,4.687013073801372E-02,3.136192147405411E-02,2.111166702275707E-02,1.405363765218121E-02,9.058574423541231E-03,5.473283988859336E-03,2.915882379323829E-03,1.180956727331740E-03,2.019486455264448E-04,5.559603061840632E-05,1.019614746518944E-03,3.740559531068072E-03,9.714608396930191E-03,2.294375667666784E-02,5.823944699356514E-02,2.492139439749798E-01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      ldouble dataCH3p[20] = {4.918065466002233E-01,1.820305598897551E-01,1.018873687066943E-01,6.626620387794423E-02,4.627140522975643E-02,3.341592505468113E-02,2.436124956610571E-02,1.755717105417992E-02,1.221477904477014E-02,7.925459395051310E-03,4.504691021186834E-03,1.932056509092962E-03,3.471322489271420E-04,9.961368457233844E-05,1.889722158404394E-03,7.119736283362036E-03,1.886587754845118E-02,4.519762639125378E-02,1.157897514436019E-01,4.979084003216981E-01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {1.781703806459979E-01,2.003979619050844E-01,2.253513609478253E-01,2.535263214572240E-01,2.855465110493468E-01,3.222076324405987E-01,3.645411392513530E-01,4.139085348461038E-01,4.721448696098857E-01,5.417842349206073E-01,6.264276710144694E-01,7.313707676913767E-01,8.647333306240219E-01,1.039631780720772E+00,1.278721347428400E+00,1.624794219081551E+00,2.169544791248932E+00,3.151365341656844E+00,5.444266144160003E+00,1.691121618992118E+01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      ldouble dataCH4p[20] = {1.750085989089947E-02,5.527513905327299E-02,9.809964127869453E-02,1.468961715042248E-01,2.028243884255907E-01,2.673638216193711E-01,3.424323399181908E-01,4.305614832578167E-01,5.351632696544314E-01,6.609494590150331E-01,8.146156068986806E-01,1.006007962955979E+00,1.250223851314107E+00,1.571650775428910E+00,2.012411582262966E+00,2.652069105529632E+00,3.661107343835269E+00,5.482752286667250E+00,9.741932698102621E+00,3.105532239473116E+01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=250\n";
      MW=1;
      AWudAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {1.072854089824838E+03,4.791769900176317E+02,2.821428609482993E+02,1.842515898865616E+02,1.260819710387350E+02,8.785958062843518E+01,6.113698058909221E+01,4.172090725352362E+01,2.731956562731075E+01,1.660766007882367E+01,8.810690411609084E+00,3.512481713708699E+00,5.789410681909526E-01,1.672955923801706E-01,2.843633190613400E+00,9.924685706522293E+00,2.445658732002864E+01,5.458633487655766E+01,1.304650666464367E+02,5.240803334856409E+02};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      ldouble dataCH5p[20] = {1.045568513883333E+03,7.187021778813055E+02,4.743315644777250E+02,3.266887241602029E+02,2.305531877647288E+02,1.639286480987088E+02,1.156998836261261E+02,7.979217928551968E+01,5.267627559596194E+01,3.223010136180344E+01,1.718911772644572E+01,6.882708760812416E+00,1.138633058879575E+00,3.300693190808542E-01,5.625741082972278E+00,1.968148730873630E+01,4.860134340260900E+01,1.086787585928923E+02,2.601816600589304E+02,1.046714965365820E+03};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      //std::cout<<"\n# mu=40, md=10, MW=1, pspatial=50\n";
      mu=40;
      md=10;
      AWudAmp.set_masses(mu,md,MW);
      pspatial = 50;
      ldouble dataCH6[20] = {5.562414050481776E+01,3.821112311073423E+01,2.712617225443829E+01,1.953847898904659E+01,1.408463757973561E+01,1.003393990428822E+01,6.964212120200670E+00,4.618908765060986E+00,2.838182007122827E+00,1.523775377295358E+00,6.212591816753744E-01,1.127470609748048E-01,1.810094371447662E-02,4.066107985862395E-01,1.427583218796996E+00,3.384889597123570E+00,6.936072551979389E+00,1.373532867381786E+01,2.932832309950583E+01,9.011852913170651E+01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      //std::cout<<"\n# A+\n";
      ldouble dataCH6p[20] = {1.639503162854304E+01,2.170717181046166E+01,2.091059906423877E+01,1.811253226260387E+01,1.479871445954067E+01,1.154699160138275E+01,8.587891208382350E+00,6.012594500328762E+00,3.858026339932550E+00,2.145292394661177E+00,9.000388041458606E-01,1.671773531817268E-01,2.733990025566899E-02,6.227787825765286E-01,2.206622360053652E+00,5.249016077007283E+00,1.069700486686444E+01,2.072622245232255E+01,4.149656485495009E+01,9.526838915611259E+01};
      i += AWudAmp.test_2to2_amp2([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      i += AWudAmp.test_2to2_amp2_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      i += AWudAmp.test_2to2_amp2_boosts([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      i += AWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AWudAmp.amp2_Aplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }





}

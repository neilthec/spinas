
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

//File:  SPINAS/SM/gWud.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gWud.h"

namespace spinas {

  gWud::gWud(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW):
    e(echarge), gs(gscharge), mu(massu), md(massd), MW(massW), SW(sinW) {
    constexpr ldouble sqrt2 = std::sqrt(2);
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
    pre = sqrt2*e*gs/(MW*SW);
  }
  void gWud::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    p2.set_mass(MW);
    p3.set_mass(mu);
    p4.set_mass(md);
    propu.set_mass(mu);
    propd.set_mass(md);
    pre = sqrt2*e*gs/(MW*SW);
  }
  void gWud::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1;
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
    pDenT=propu.denominator(propTP);
    pDenU=propd.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gWud::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    constexpr ldouble one=1, two=2, three=3;
    cdouble amplitude(0,0);
    int ds2a, ds2b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);

      if(ds1>0){
      
	//TU Diagram
	//pre = sqrt(2)*e*gs/(MW*SW);
	//all in: - pre <24>(-MW[13][142>+[12](Md^2 [13]-Mu[143>))/(t-Mu^2)/(u-Md^2)
	amplitude += normFactor*pre*a24a.v(ds2a,ds4)*(MW*s13s.v(ds3)*s142a.v(ds2b)+s12s.v(ds2b)*(md*md*s13s.v(ds3)-mu*s143a.v(ds3)))/pDenT/pDenU;
	
      }
      else if(ds1<0){

	//TU Diagram
	amplitude += normFactor*pre*s23s.v(ds2a,ds3)*(-MW*a14a.v(ds4)*s231a.v(ds2b)+a12a.v(ds2b)*(-mu*mu*a14a.v(ds4)+md*s431a.v(ds4)))/pDenT/pDenU;
	
      }
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble gWud::amp2(){
    constexpr ldouble three=3, four=4;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    //Color factor Tr(Ta,Ta) = 4
	    amp2 += four*std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/8
    return amp2/48.0;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble gWud::amp2_gplus(){
    constexpr ldouble one=1, three=3, four=4;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(1,j2,j3,j4);
	  //Color factor Tr(Ta,Ta) = 4
	  amp2 += four*std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    //Average over initial colors 1/8
    return amp2/24.0;
  }


  



  //  Tests
  int test_gWud(){
    int n=0;//Number of fails
    std::cout<<"\t* g , W+ -> u , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,gs=1.238,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      gWud gWudAmp = gWud(EE,gs,mu,md,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.592713647798128E+00,1.468853394156237E+00,8.500315520431206E-01,5.898450955609778E-01,4.499599038683798E-01,3.655464261093249E-01,3.118681120401530E-01,2.776213476589802E-01,2.571258416047644E-01,2.474944947676303E-01,2.474944947530804E-01,2.571258415593099E-01,2.776213475766203E-01,3.118681119085416E-01,3.655464259044592E-01,4.499599035410448E-01,5.898450949966552E-01,8.500315509085712E-01,1.468853390960206E+00,4.592713618862534E+00};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH);
      ldouble dataCHp[20] = {2.309019804081973E-02,3.231457752764851E-02,4.605268370790316E-02,6.198018145214457E-02,8.014883461063489E-02,1.009138823137345E-01,1.248075006128169E-01,1.525579669987445E-01,1.851584224665634E-01,2.239863393311113E-01,2.710026472400843E-01,3.290932576705731E-01,4.026847250428834E-01,4.989287198558811E-01,6.301789657543218E-01,8.197709680776405E-01,1.117710002983559E+00,1.654010410132242E+00,2.905392189359308E+00,9.162336961213230E+00};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCHp);
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=81\n";
      pspatial = 81;
      ldouble dataCH2[20] = {6.610974679907356E+00,2.177972171207703E+00,1.299815466491087E+00,9.305904866441964E-01,7.320824325219054E-01,6.122930893491844E-01,5.361193530908622E-01,4.875205092954631E-01,4.584357791854787E-01,4.447681427668814E-01,4.447681425862296E-01,4.584357786211188E-01,4.875205082728864E-01,5.361193514567865E-01,6.122930868055801E-01,7.320824284577281E-01,9.305904796375886E-01,1.299815452404571E+00,2.177972131525904E+00,6.610974320644110E+00};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH2);
      ldouble dataCH2p[20] = {4.820741032069530E-01,2.512480272924283E-01,2.214975766351876E-01,2.226039908877896E-01,2.360819665622161E-01,2.573608656790696E-01,2.852286732233075E-01,3.197821682672080E-01,3.619111573748696E-01,4.132383774131674E-01,4.762979028756059E-01,5.549603951186770E-01,6.552588434394646E-01,7.870100244871485E-01,9.672253019614648E-01,1.228082882900338E+00,1.638576957916766E+00,2.378133310432857E+00,4.104696192423911E+00,1.273987418925262E+01};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH2p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      gWudAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {3.376760254490281E+00,1.673617852814595E+00,1.113490355658148E+00,8.400275661180000E-01,6.822317041793488E-01,5.833461303944317E-01,5.193233362879244E-01,4.784766643549087E-01,4.547428224583208E-01,4.450539384154302E-01,4.482486695175341E-01,4.647042446255404E-01,4.964532625235959E-01,5.478577830199475E-01,6.273020903818139E-01,7.513372391873429E-01,9.560778339840513E-01,1.336077099690137E+00,2.238834883465498E+00,6.793660322936635E+00};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH3);
      ldouble dataCH3p[20] = {3.351328830260934E+00,1.462166243927212E+00,9.790554126127078E-01,7.753199922114065E-01,6.735167823858258E-01,6.215515395102718E-01,5.992594231036031E-01,5.977595871919139E-01,6.131851259165435E-01,6.444498265938218E-01,6.924908120960714E-01,7.602605919867023E-01,8.533602048074783E-01,9.816192231928026E-01,1.162622121987808E+00,1.430086316887397E+00,1.856713786675262E+00,2.631980212861764E+00,4.451177819534895E+00,1.357315922925233E+01};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH3p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {5.141007970738404E+00,5.717655419169828E+00,6.358052750390573E+00,7.073808685802049E+00,7.879527524588099E+00,8.793842087723768E+00,9.840905740584677E+00,1.105260057357862E+01,1.247189765460525E+01,1.415813738488681E+01,1.619564473790717E+01,1.870842568359595E+01,2.188662040157539E+01,2.603737400572849E+01,3.169120035028305E+01,3.985015720286491E+01,5.266143741387666E+01,7.570743442338085E+01,1.294547055152235E+02,3.980291241818096E+02};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH4);
      ldouble dataCH4p[20] = {5.049776504247177E-01,1.577082897195358E+00,2.767778687561348E+00,4.098649039377870E+00,5.596847761803000E+00,7.297019035464809E+00,9.244071565100487E+00,1.149728429394575E+01,1.413655417488461E+01,1.727221399598816E+01,2.106105074496449E+01,2.573363066106254E+01,3.164348346670208E+01,3.936168536312434E+01,4.987461792906119E+01,6.504538822642903E+01,8.886618798124803E+01,1.317159593385225E+02,2.316453669214553E+02,7.309304449271654E+02};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH4p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=250\n";
      MW=1;
      gWudAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {7.339545633102079E+03,3.864608148479344E+03,2.722538049324045E+03,2.165168056578138E+03,1.843616905971035E+03,1.642134946926608E+03,1.511693366580199E+03,1.428469303659195E+03,1.380103965640947E+03,1.360343654025241E+03,1.366819351740668E+03,1.400295967330374E+03,1.464909248711849E+03,1.569537549723511E+03,1.731247451089652E+03,1.983730446721180E+03,2.400502670403283E+03,3.174040738797079E+03,5.011731231570623E+03,1.428374087362452E+04};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH5);
      ldouble dataCH5p[20] = {7.152881172717868E+03,5.796401644552543E+03,4.577063293203948E+03,3.838968176238936E+03,3.371233422088411E+03,3.063899917571422E+03,2.860833899653200E+03,2.731979870159897E+03,2.661050246439831E+03,2.639987430380258E+03,2.666580897780219E+03,2.743880283977605E+03,2.881112069070497E+03,3.096651758358793E+03,3.425037358029335E+03,3.933904484798876E+03,4.770397974847871E+03,6.319361942797659E+03,9.994710347505890E+03,2.852807933927018E+04};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH5p);
      //std::cout<<"\n# mu=80, md=0.0075, MW=1, pspatial=50\n";
      pspatial = 50;
      ldouble dataCH6[20] = {1.540125756724649E+03,1.524742965112971E+03,1.516051872253968E+03,1.514344814282833E+03,1.520112309337081E+03,1.534092304913569E+03,1.557344704427398E+03,1.591364653825984E+03,1.638257595871353E+03,1.701016756401846E+03,1.783978080732569E+03,1.893598340893143E+03,2.039857657234458E+03,2.238958519660149E+03,2.518971048497860E+03,2.933007181124116E+03,3.594981259267404E+03,4.801198498073660E+03,7.638300809661729E+03,2.188656970179591E+04};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH6);
      ldouble dataCH6p[20] = {1.578635479278066E+01,4.701097378750048E+01,7.844966653254411E+01,1.105471652784770E+02,1.438060828595358E+02,1.788160257273157E+02,2.162939816005052E+02,2.571426327143960E+02,3.025377971929054E+02,3.540646722398313E+02,4.139390647536542E+02,4.853838032696723E+02,5.733053549159384E+02,6.855941216133941E+02,8.358423596959116E+02,1.049685085715449E+03,1.382010119313345E+03,1.975471037023961E+03,3.353050024479847E+03,1.022135857019656E+04};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH6p);
      //std::cout<<"\n# mu=0.0042, md=0.0042, MW=80.385, pspatial=250\n";
      mu = 0.0042;
      md = 0.0042;
      MW = 80.385;
      gWudAmp.set_masses(mu,md,MW);
      pspatial = 250;
      ldouble dataCH7[20] = {4.592713634170686E+00,1.468853389382679E+00,8.500315490255421E-01,5.898450932836159E-01,4.499599019909377E-01,3.655464244753756E-01,3.118681105635959E-01,2.776213462860743E-01,2.571258402983417E-01,2.474944934994213E-01,2.474944934994213E-01,2.571258402983416E-01,2.776213462860742E-01,3.118681105635959E-01,3.655464244753757E-01,4.499599019909375E-01,5.898450932836157E-01,8.500315490255418E-01,1.468853389382678E+00,4.592713634170678E+00};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH7);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH7);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH7);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2(); }, 0,MW,mu,md,pspatial,dataCH7);
      ldouble dataCH7p[20] = {2.309019812933420E-02,3.231457758901838E-02,4.605268377002356E-02,6.198018152035332E-02,8.014883468815232E-02,1.009138824036671E-01,1.248075007188116E-01,1.525579671254735E-01,1.851584226203141E-01,2.239863395206691E-01,2.710026474781734E-01,3.290932579763690E-01,4.026847254466749E-01,4.989287204083802E-01,6.301789665470843E-01,8.197709692937227E-01,1.117710005046878E+00,1.654010414281060E+00,2.905392201176337E+00,9.162337070212020E+00};
      i += gWudAmp.test_2to2_amp2([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH7p);
      i += gWudAmp.test_2to2_amp2_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH7p);
      i += gWudAmp.test_2to2_amp2_boosts([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH7p);
      i += gWudAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gWudAmp.amp2_gplus(); }, 0,MW,mu,md,pspatial,dataCH7p);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  



}

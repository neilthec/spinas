
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

//File:  SPINAS/SM/gdZd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gdZd.h"

namespace spinas {

  gdZd::gdZd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qd(-1.0/3.0), gs(gscharge), md(massd), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), propd(massd,0) {
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
    preTU = 2.0*e*gs/(2.0*MW*SW);
    gL=-2.0*Qd*SW*SW-1.0;
    gR=-2.0*Qd*SW*SW;
  }
  void gdZd::set_masses(const ldouble& massd, const ldouble& massW){
    md=massd;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p2.set_mass(md);
    p3.set_mass(MZ);
    p4.set_mass(md);
    propd.set_mass(md);
    //Couplings
    preTU = 2.0*e*gs/(2.0*MW*SW);
  }
  void gdZd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenS=propd.denominator(propSP);
    pDenU=propd.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gdZd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //preTU = 2.0*e*gs/(2.0*MW*SW);
      
      if(ds1>0){
	//Same as photon
	amplitude += normFactor*preTU*(
				       +gL*a23a.v(ds2,ds3a)*(MZ*s14s.v(ds4)*s123a.v(ds3b)+md*s13s.v(ds3b)*s134a.v(ds4))
				       +gR*a34a.v(ds3a,ds4)*(md*s13s.v(ds3b)*s132a.v(ds2)+MZ*s12s.v(ds2)*s143a.v(ds3b))
				       )/pDenU/pDenS;
	
      }
      else if(ds1<0){
	amplitude += normFactor*preTU*(
				       +gR*s23s.v(ds2,ds3a)*(MZ*a14a.v(ds4)*a123s.v(ds3b)+md*a13a.v(ds3b)*a134s.v(ds4))
				       +gL*s34s.v(ds3a,ds4)*(md*a13a.v(ds3b)*a132s.v(ds2)+MZ*a12a.v(ds2)*a143s.v(ds3b))
				       )/pDenU/pDenS;

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble gdZd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 // It's the same for both diagrams
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/3=1/24
    return amp2/96.0;
  }
  //set_momenta(...) must be called before amp2_gplus().
  ldouble gdZd::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3
    return amp2/48.0;
  }  
  //set_momenta(...) must be called before amp2_gminus().
  ldouble gdZd::amp2_gminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3
    return amp2/48.0;
  }



  //  Tests
  int test_gdZd(){
    int n=0;//Number of fails
    std::cout<<"\t* g , d  -> Z , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, MW=80.385, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333, gs=1.238, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      gdZd gdZdAmp = gdZd(EE,gs,md,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {5.417310334288764E-02,5.174802664169696E-02,4.933322009874962E-02,4.693055101556145E-02,4.454236857791501E-02,4.217167002285005E-02,3.982234065765014E-02,3.749950908794453E-02,3.521008778221341E-02,3.296362253030039E-02,3.077367840204685E-02,2.866020402968529E-02,2.665378731618550E-02,2.480383948985689E-02,2.319570721876132E-02,2.199060087843827E-02,2.153397128152646E-02,2.272399057411408E-02,2.885395686516778E-02,6.956356105461241E-02};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCHp);
      ldouble dataCHm[20] = {5.322017312238509E-02,5.593079008057320E-02,5.896011239215258E-02,6.236608648497933E-02,6.622161270376632E-02,7.061970183313433E-02,7.568092340868358E-02,8.156442921737662E-02,8.848472818942248E-02,9.673804694265443E-02,1.067453390917645E-01,1.191256540699019E-01,1.348282010074792E-01,1.553863177296956E-01,1.834484968417933E-01,2.240174465931635E-01,2.878132560741729E-01,4.027081568547908E-01,6.709003308889204E-01,2.012173400274124E+00};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCHm);
      ldouble dataCH[20] = {5.369663823263637E-02,5.383940836113508E-02,5.414666624545111E-02,5.464831875027039E-02,5.538199064084066E-02,5.639568592799219E-02,5.775163203316686E-02,5.953196915266058E-02,6.184740798581794E-02,6.485083473647742E-02,6.875950874690566E-02,7.389292904979358E-02,8.074099416183235E-02,9.009507860977624E-02,1.033221020302773E-01,1.230040237358009E-01,1.546736136778497E-01,2.127160737144524E-01,3.498771438770441E-01,1.040868480664368E+00};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {6.036031759673340E-02,5.824161555472386E-02,5.613852458923881E-02,5.405388307792016E-02,5.199126188298783E-02,4.995521693213804E-02,4.795165405762312E-02,4.598836895809695E-02,4.407586887969739E-02,4.222866382919965E-02,4.046737329025981E-02,3.882232003357083E-02,3.733999897534840E-02,3.609551729029507E-02,3.521860555820837E-02,3.495431044352151E-02,3.582773198424011E-02,3.920413045890393E-02,5.008946152017252E-02,1.135373688764366E-01};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {4.809135066147852E-02,5.052386280151468E-02,5.324237593091731E-02,5.629889022942357E-02,5.975882527795259E-02,6.370564743830841E-02,6.824755384607709E-02,7.352736472010250E-02,7.973757687688492E-02,8.714401925316898E-02,9.612444875851353E-02,1.072343902589373E-01,1.213256485471518E-01,1.397742159768363E-01,1.649568066049139E-01,2.013627691835807E-01,2.586121366845885E-01,3.617169945460853E-01,6.023883206250086E-01,1.806024846896809E+00};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH2m);
      ldouble dataCH2[20] = {5.422583412910596E-02,5.438273917811928E-02,5.469045026007806E-02,5.517638665367187E-02,5.587504358047021E-02,5.683043218522323E-02,5.809960395185011E-02,5.975786683909973E-02,6.190672287829116E-02,6.468634154118431E-02,6.829591102438667E-02,7.302835514625408E-02,7.933282376125009E-02,8.793486663356570E-02,1.000877060815611E-01,1.181585398135511E-01,1.472199343344143E-01,2.004605625024946E-01,3.262388910725905E-01,9.597811078866229E-01};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, MW=80.385, pspatial=95\n";
      md = 125;
      pspatial = 250;
      gdZdAmp.set_masses(md,MW);
      ldouble dataCH4p[20] = {5.452317597790258E-02,5.315435107585709E-02,5.225225229668068E-02,5.189891391380267E-02,5.219675658520657E-02,5.327532753982714E-02,5.530089886701143E-02,5.849044234914340E-02,6.313249483747205E-02,6.961922178961562E-02,7.849735665026224E-02,9.055234149370689E-02,1.069538774869440E-01,1.295221256715046E-01,1.612491238536008E-01,2.074125640500032E-01,2.782405083426862E-01,3.963602944340967E-01,6.229989829200843E-01,1.193409599345149E+00};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {5.358814559156379E-02,5.694248380411524E-02,6.095361745004677E-02,6.572859617427232E-02,7.139903366914163E-02,7.812856872084807E-02,8.612321511468658E-02,9.564599143393397E-02,1.070380408045430E-01,1.207498523825494E-01,1.373886776840240E-01,1.577927969540241E-01,1.831520455521149E-01,2.152116421545362E-01,2.566338353469306E-01,3.116762622261320E-01,3.875461466451344E-01,4.972668871946428E-01,6.656574301210328E-01,9.314268606062426E-01};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH4m);
      ldouble dataCH4[20] = {5.405566078473319E-02,5.504841743998617E-02,5.660293487336372E-02,5.881375504403749E-02,6.179789512717410E-02,6.570194813033760E-02,7.071205699084901E-02,7.706821689153869E-02,8.508526782100752E-02,9.518453708608249E-02,1.079430171671431E-01,1.241725692238655E-01,1.450529615195295E-01,1.723668839130204E-01,2.089414796002657E-01,2.595444131380676E-01,3.328933274939103E-01,4.468135908143698E-01,6.443282065205586E-01,1.062418229975696E+00};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH4);
      //std::cout<<"\n# md=125, MW=80.385, pspatial=125.1\n";
      md = 125;
      pspatial = 125.1;
      gdZdAmp.set_masses(md,MW);
      ldouble dataCH3p[20] = {6.112792078788584E-02,6.150828966166733E-02,6.234834969531349E-02,6.371572431028172E-02,6.569189584786385E-02,6.837592513628341E-02,7.188942968719550E-02,7.638334418878472E-02,8.204725028150185E-02,8.912248403401377E-02,9.792092152480852E-02,1.088525121285293E-01,1.224666690731648E-01,1.395163178763726E-01,1.610603647774523E-01,1.886341210395309E-01,2.245460377475915E-01,2.724235288067042E-01,3.382866326888235E-01,4.328448866351339E-01};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gplus(); }, 0,md,MZ,md,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {5.049141740266485E-02,5.387676268192640E-02,5.766011396207456E-02,6.188537509224187E-02,6.660275788038024E-02,7.186982371971676E-02,7.775267757773308E-02,8.432730171565045E-02,9.168097505748193E-02,9.991363834185137E-02,1.091388907088040E-01,1.194839486604295E-01,1.310871690097095E-01,1.440902062769812E-01,1.586185637753454E-01,1.747368549315207E-01,1.923475288786162E-01,2.109577880740060E-01,2.291204162768782E-01,2.430004295042004E-01};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2_gminus(); }, 0,md,MZ,md,pspatial,dataCH3m);
      ldouble dataCH3[20] = {5.580966909527534E-02,5.769252617179686E-02,6.000423182869403E-02,6.280054970126180E-02,6.614732686412204E-02,7.012287442800008E-02,7.482105363246429E-02,8.035532295221758E-02,8.686411266949189E-02,9.451806118793257E-02,1.035299061168063E-01,1.141682303944794E-01,1.267769190414372E-01,1.418032620766769E-01,1.598394642763988E-01,1.816854879855259E-01,2.084467833131038E-01,2.416906584403551E-01,2.837035244828508E-01,3.379226580696671E-01};
      i += gdZdAmp.test_2to2_amp2([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      i += gdZdAmp.test_2to2_amp2_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      i += gdZdAmp.test_2to2_amp2_boosts([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      i += gdZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gdZdAmp.amp2(); }, 0,md,MZ,md,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

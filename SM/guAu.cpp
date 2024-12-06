
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

//File:  SPINAS/SM/guAu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/guAu.h"

namespace spinas {

  guAu::guAu(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu):
    e(echarge), Qu(2.0/3.0), gs(gscharge), mu(massu), prop(massu,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(mu);
    p3=particle(0);
    p4=particle(mu);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    s143a = sproduct(SQUARE,&p1,&p4,&p3,2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1,2);
  }
  void guAu::set_masses(const ldouble& massu){
    mu=massu;
    p2.set_mass(mu);
    p4.set_mass(mu);
    prop.set_mass(mu);
  }
  void guAu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s143a.update();
    s341a.update();
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS = prop.denominator(propSP);
    pDenU = prop.denominator(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble guAu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds3>0){
      //gAUu all in:      + mu[12]^2<34>
      //guAU: 4->2->3->4: - mu[13]^2<24>
      //34 out:           + mu[13]^2<24>
      return + 2.0*e*Qu*gs*mu*s13s.v()*s13s.v()*a24a.v(ds2,ds4)/pDenU/pDenS;
    }
    else if(ds1<0&&ds3<0){
      //gAUu all in:      + mu<12>^2[34]
      //guAU: 4->2->3->4: - mu<13>^2[24]
      //34 out:           - mu<13>^2[24]
      return - 2.0*e*Qu*gs*mu*a13a.v()*a13a.v()*s24s.v(ds2,ds4)/pDenU/pDenS;
    }
    else if(ds1>0&&ds3<0){
      //gAUu all in:      + ([31]<42>+[41]<32>)[132>
      //guAU: 4->2->3->4: - ([14]<23>-[12]<34>)[143>
      //34 out:           + ([14]<23>+[12]<34>)[143>
      return + 2.0*e*Qu*gs*(s14s.v(ds4)*a23a.v(ds2)+s12s.v(ds2)*a34a.v(ds4))*s143a.v()/pDenU/pDenS;
    }
    else if(ds1<0&&ds3>0){
      //gAUu all in:      + (<31>[42]+<41>[32])[231>
      //guAU: 4->2->3->4: - (<14>[23]-<12>[34])[341>
      //34 out:           - (<14>[23]+<12>[34])[341>
      return - 2.0*e*Qu*gs*(a14a.v(ds4)*s23s.v(ds2)+a12a.v(ds2)*s34s.v(ds4))*s341a.v()/pDenU/pDenS;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble guAu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 //The color factor is the same for both diagrams
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/3=1/24
    return amp2/96.0;
  }

  //g+, u -> A, u
  ldouble guAu::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3
    return amp2/48.0;
  }

  //g-, u -> A, u
  ldouble guAu::amp2_gminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3
    return amp2/48.0;
  }



  //  Tests
  int test_guAu(){
    int n=0;//Number of fails
    std::cout<<"\t* g , u  -> A , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  mu=0.0042, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333, gs=1.238;
      guAu guAuAmp = guAu(EE,gs,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.459751335790951E-02,4.471878096229907E-02,4.498128836106490E-02,4.541071551681868E-02,4.603936947929705E-02,4.690846958781835E-02,4.807144831927541E-02,4.959883654316487E-02,5.158569760532697E-02,5.416330945953626E-02,5.751822498746414E-02,6.192478666214221E-02,6.780365294219140E-02,7.583434902994315E-02,8.719060027475228E-02,1.040894434472951E-01,1.312816716265705E-01,1.811193468136855E-01,2.988933629817089E-01,8.922217648656300E-01};
      i += guAuAmp.test_2to2_amp2([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guAuAmp.test_2to2_amp2_rotations([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guAuAmp.test_2to2_amp2_boosts([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      i += guAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH);
      //std::cout<<"\n#  g+ , u -> A , u\n";
      ldouble dataCHpp[20] = {4.459751335790951E-02,4.471878096229907E-02,4.498128836106490E-02,4.541071551681868E-02,4.603936947929705E-02,4.690846958781835E-02,4.807144831927541E-02,4.959883654316487E-02,5.158569760532697E-02,5.416330945953626E-02,5.751822498746414E-02,6.192478666214221E-02,6.780365294219140E-02,7.583434902994315E-02,8.719060027475228E-02,1.040894434472951E-01,1.312816716265705E-01,1.811193468136855E-01,2.988933629817089E-01,8.922217648656300E-01};
      i += guAuAmp.test_2to2_amp2([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCHpp);
      i += guAuAmp.test_2to2_amp2_rotations([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCHpp);
      i += guAuAmp.test_2to2_amp2_boosts([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCHpp);
      i += guAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCHpp);
      //std::cout<<"\n#  g- , u -> A , u\n";
      ldouble dataCHpm[20] = {4.459751335790951E-02,4.471878096229907E-02,4.498128836106490E-02,4.541071551681868E-02,4.603936947929705E-02,4.690846958781835E-02,4.807144831927541E-02,4.959883654316487E-02,5.158569760532697E-02,5.416330945953626E-02,5.751822498746414E-02,6.192478666214221E-02,6.780365294219140E-02,7.583434902994315E-02,8.719060027475228E-02,1.040894434472951E-01,1.312816716265705E-01,1.811193468136855E-01,2.988933629817089E-01,8.922217648656300E-01};
      i += guAuAmp.test_2to2_amp2([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCHpm);
      i += guAuAmp.test_2to2_amp2_rotations([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCHpm);
      i += guAuAmp.test_2to2_amp2_boosts([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCHpm);
      i += guAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  mu=0.0042, pspatial=0.006\n";
      pspatial = 0.006;
      ldouble dataCH2[20] = {4.436875822548277E-02,4.398535762283213E-02,4.367105642901919E-02,4.343797229257763E-02,4.330120757091645E-02,4.327982103894069E-02,4.339820266612876E-02,4.368805951468625E-02,4.419135238955135E-02,4.496475539712769E-02,4.608663746294708E-02,4.766838359199867E-02,4.987352540521237E-02,5.295168662176965E-02,5.730247428066550E-02,6.360479021357489E-02,7.310371910986031E-02,8.832833407161066E-02,1.152105507350654E-01,1.710754391092357E-01};
      i += guAuAmp.test_2to2_amp2([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guAuAmp.test_2to2_amp2_rotations([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guAuAmp.test_2to2_amp2_boosts([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      i += guAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guAuAmp.amp2(); }, 0,mu,0,mu,pspatial,dataCH2);
      //std::cout<<"\n#  g+ , u -> A , u\n";
      ldouble dataCH2pp[20] = {4.436875822548277E-02,4.398535762283213E-02,4.367105642901919E-02,4.343797229257763E-02,4.330120757091645E-02,4.327982103894069E-02,4.339820266612876E-02,4.368805951468625E-02,4.419135238955135E-02,4.496475539712769E-02,4.608663746294708E-02,4.766838359199867E-02,4.987352540521237E-02,5.295168662176965E-02,5.730247428066550E-02,6.360479021357489E-02,7.310371910986031E-02,8.832833407161066E-02,1.152105507350654E-01,1.710754391092357E-01};
      i += guAuAmp.test_2to2_amp2([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2pp);
      i += guAuAmp.test_2to2_amp2_rotations([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2pp);
      i += guAuAmp.test_2to2_amp2_boosts([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2pp);
      i += guAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guAuAmp.amp2_gplus(); }, 0,mu,0,mu,pspatial,dataCH2pp);
      //std::cout<<"\n#  g- , u -> A , u\n";
      ldouble dataCH2pm[20] = {4.436875822548277E-02,4.398535762283213E-02,4.367105642901919E-02,4.343797229257763E-02,4.330120757091645E-02,4.327982103894069E-02,4.339820266612876E-02,4.368805951468625E-02,4.419135238955135E-02,4.496475539712769E-02,4.608663746294708E-02,4.766838359199867E-02,4.987352540521237E-02,5.295168662176965E-02,5.730247428066550E-02,6.360479021357489E-02,7.310371910986031E-02,8.832833407161066E-02,1.152105507350654E-01,1.710754391092357E-01};
      i += guAuAmp.test_2to2_amp2([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2pm);
      i += guAuAmp.test_2to2_amp2_rotations([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2pm);
      i += guAuAmp.test_2to2_amp2_boosts([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2pm);
      i += guAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guAuAmp.amp2_gminus(); }, 0,mu,0,mu,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

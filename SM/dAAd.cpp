
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

//File:  SPINAS/SM/dAAd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dAAd.h"

namespace spinas {

  dAAd::dAAd(const ldouble& echarge, const ldouble& massd):
    e(echarge), Qd(-1.0/3.0), md(massd), prop(massd,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(md);
    p2=particle(0);
    p3=particle(0);
    p4=particle(md);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s342a = sproduct(SQUARE,&p3,&p4,&p2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
  }
  void dAAd::set_masses(const ldouble& massd){
    md=massd;
    p1.set_mass(md);
    p4.set_mass(md);
    prop.set_mass(md);
  }
  void dAAd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s342a.update();
    s243a.update();
    //Propagator Momentum
    ldouble propTP[4], propSP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propSP[j] = mom1[j]+mom2[j];
    }
    pDenT = prop.denominator(propTP);
    pDenS = prop.denominator(propSP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dAAd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds3>0&&ds2>0){
      //AADd all in:      + md[12]^2<34>/tu
      //dAAD: 4->1->3->4: - md[23]^2<14>/st
      //34 out:           + md[23]^2<14>/st
      return + 2.0*e*e*Qd*Qd*md*s23s.v()*s23s.v()*a14a.v(ds1,ds4)/pDenS/pDenT;
    }
    else if(ds3<0&&ds2<0){
      //AADd all in:      + md<12>^2[34]/tu
      //dAAD: 4->1->3->4: - md<23>^2[14]/st
      //34 out:           - md<23>^2[14]/st
      return - 2.0*e*e*Qd*Qd*md*a23a.v()*a23a.v()*s14s.v(ds1,ds4)/pDenS/pDenT;
    }
    else if(ds3>0&&ds2<0){
      //AADd all in:      + ([31]<42>+[41]<32>)[132>/tu
      //dAAD: 4->1->3->4: - ([34]<12>+[13]<24>)[342>/st
      //34 out:           - ([34]<12>-[13]<24>)[342>/st
      return - 2.0*e*e*Qd*Qd*(s34s.v(ds4)*a12a.v(ds1)-s13s.v(ds1)*a24a.v(ds4))*s342a.v()/pDenS/pDenT;
    }
    else if(ds3<0&&ds2>0){
      //AADd all in:      + (<31>[42]+<41>[32])[231>/tu
      //dAAD: 4->1->3->4: + (-<34>[12]-<13>[24])[243>/st
      //34 out:           - (<34>[12]-<13>[24])[243>/st
      return - 2.0*e*e*Qd*Qd*(a34a.v(ds4)*s12s.v(ds1)-a13a.v(ds1)*s24s.v(ds4))*s243a.v()/pDenS/pDenT;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble dAAd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2=1/4
    //average over initial colors 1/3
    return amp2/12.0;
  }

  //d, A+ -> A, d
  ldouble dAAd::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2
    //average over initial colors 1/3
    return amp2/6.0;
  }

  //d, A- -> A, d
  ldouble dAAd::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,-2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2
    //average over initial colors 1/3
    return amp2/6.0;
  }



  //  Tests
  int test_dAAd(){
    int n=0;//Number of fails
    std::cout<<"\t* d , A  -> A , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  md=0.0075, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333;
      dAAd dAAdAmp = dAAd(EE,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {9.525422063396554E-03,3.191006494487596E-03,1.933642844049899E-03,1.401572330417499E-03,1.111266196203725E-03,9.308529618335877E-04,8.096128274143672E-04,7.238765529798987E-04,6.611133644245504E-04,6.140686030530894E-04,5.782512896720440E-04,5.507325248067996E-04,5.295206568142811E-04,5.132141530565722E-04,5.007981105888547E-04,4.915195369036069E-04,4.848079831314445E-04,4.802233887159518E-04,4.774208413291041E-04,4.761261798896345E-04};
      i += dAAdAmp.test_2to2_amp2([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dAAdAmp.test_2to2_amp2_rotations([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dAAdAmp.test_2to2_amp2_boosts([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dAAdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      //std::cout<<"\n#  d , A+  -> A , d\n";
      ldouble dataCHpp[20] = {9.525422063396554E-03,3.191006494487596E-03,1.933642844049899E-03,1.401572330417499E-03,1.111266196203725E-03,9.308529618335877E-04,8.096128274143672E-04,7.238765529798987E-04,6.611133644245504E-04,6.140686030530894E-04,5.782512896720440E-04,5.507325248067996E-04,5.295206568142811E-04,5.132141530565722E-04,5.007981105888547E-04,4.915195369036069E-04,4.848079831314445E-04,4.802233887159518E-04,4.774208413291041E-04,4.761261798896345E-04};
      i += dAAdAmp.test_2to2_amp2([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dAAdAmp.test_2to2_amp2_rotations([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dAAdAmp.test_2to2_amp2_boosts([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dAAdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCHpp);
      //std::cout<<"\n#  d , A-  -> A , d\n";
      ldouble dataCHpm[20] = {9.525422063396554E-03,3.191006494487596E-03,1.933642844049899E-03,1.401572330417499E-03,1.111266196203725E-03,9.308529618335877E-04,8.096128274143672E-04,7.238765529798987E-04,6.611133644245504E-04,6.140686030530894E-04,5.782512896720440E-04,5.507325248067996E-04,5.295206568142811E-04,5.132141530565722E-04,5.007981105888547E-04,4.915195369036069E-04,4.848079831314445E-04,4.802233887159518E-04,4.774208413291041E-04,4.761261798896345E-04};
      i += dAAdAmp.test_2to2_amp2([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dAAdAmp.test_2to2_amp2_rotations([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dAAdAmp.test_2to2_amp2_boosts([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dAAdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  md=0.0075, pspatial=0.008\n";
      pspatial = 0.008;
      ldouble dataCH2[20] = {1.268028815924107E-03,9.216971055679197E-04,7.337431583592343E-04,6.226291384950725E-04,5.532389390034049E-04,5.084341934708469E-04,4.790452458856723E-04,4.598070451541048E-04,4.475222511930997E-04,4.401585997620001E-04,4.363745162918414E-04,4.352556203606773E-04,4.361614224875212E-04,4.386325413475884E-04,4.423326065373548E-04,4.470107860130930E-04,4.524769782962673E-04,4.585850053039990E-04,4.652209883251504E-04,4.722951584060094E-04};
      i += dAAdAmp.test_2to2_amp2([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dAAdAmp.test_2to2_amp2_rotations([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dAAdAmp.test_2to2_amp2_boosts([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dAAdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAAdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      //std::cout<<"\n#  d , A+  -> A , d\n";
      ldouble dataCH2pp[20] = {1.268028815924107E-03,9.216971055679197E-04,7.337431583592343E-04,6.226291384950725E-04,5.532389390034049E-04,5.084341934708469E-04,4.790452458856723E-04,4.598070451541048E-04,4.475222511930997E-04,4.401585997620001E-04,4.363745162918414E-04,4.352556203606773E-04,4.361614224875212E-04,4.386325413475884E-04,4.423326065373548E-04,4.470107860130930E-04,4.524769782962673E-04,4.585850053039990E-04,4.652209883251504E-04,4.722951584060094E-04};
      i += dAAdAmp.test_2to2_amp2([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dAAdAmp.test_2to2_amp2_rotations([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dAAdAmp.test_2to2_amp2_boosts([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dAAdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAAdAmp.amp2_Aplus(); }, md,0,0,md,pspatial,dataCH2pp);
      //std::cout<<"\n#  d , A-  -> A , d\n";
      ldouble dataCH2pm[20] = {1.268028815924107E-03,9.216971055679197E-04,7.337431583592343E-04,6.226291384950725E-04,5.532389390034049E-04,5.084341934708469E-04,4.790452458856723E-04,4.598070451541048E-04,4.475222511930997E-04,4.401585997620001E-04,4.363745162918414E-04,4.352556203606773E-04,4.361614224875212E-04,4.386325413475884E-04,4.423326065373548E-04,4.470107860130930E-04,4.524769782962673E-04,4.585850053039990E-04,4.652209883251504E-04,4.722951584060094E-04};
      i += dAAdAmp.test_2to2_amp2([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dAAdAmp.test_2to2_amp2_rotations([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dAAdAmp.test_2to2_amp2_boosts([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dAAdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dAAdAmp.amp2_Aminus(); }, md,0,0,md,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}

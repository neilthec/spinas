
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

//File:  SPINAS/SM/dggd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dggd.h"

namespace spinas {

  dggd::dggd(const ldouble& gcharge, const ldouble& massd):
    gs(gcharge), md(massd), propd(massd,0), propg(0,0){
    constexpr ldouble two=2;
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
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
  }
  void dggd::set_masses(const ldouble& massd){
    md=massd;
    p1.set_mass(md);
    p4.set_mass(md);
    propd.set_mass(md);
  }
  void dggd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS = std::real(propd.denominator(propSP));
    pDenT = std::real(propd.denominator(propTP));
    pDenU = std::real(propg.denominator(propUP));
    
  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dggd::amp_Num(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    ////////////////////////////////////////////////////
    //Missing color factor includig averaging over color
    //Missing symmetry factor when identical
    ////////////////////////////////////////////////////
    if(ds3>0&&ds2>0){
      //ggDd all in:   + md[12]^2<34>
      //Dggd 1<->3:    + md[23]^2<14>
      //34 out:        - md[23]^2<14>
      return - 2.0*gs*gs*md*s23s.v()*s23s.v()*a14a.v(ds1,ds4);
    }
    else if(ds3<0&&ds2<0){
      //ggDd all in:   + md<12>^2[34]
      //Dggd 1<->3:    + md<23>^2[14]
      //34 out:        + md<23>^2[14]
      return + 2.0*gs*gs*md*a23a.v()*a23a.v()*s14s.v(ds1,ds4);
    }
    else if(ds3>0&&ds2<0){
      //ggDd all in:   + ([31]<42>+[41]<32>)[132>
      //Dggd 1<->3:    - ([13]<24>+[34]<12>)[312>
      //34 out:        + ([13]<24>-[34]<12>)[312>
      return + 2.0*gs*gs*(s13s.v(ds1)*a24a.v(ds4)-s34s.v(ds4)*a12a.v(ds1))*s312a.v();
    }
    else if(ds3<0&&ds2>0){
      //ggDd all in:   + (<31>[42]+<41>[32])[231>
      //Dggd 1<->3:    - (<13>[24]+<34>[12])[213>
      //34 out:        - (<13>[24]-<34>[12])[213>
      return - 2.0*gs*gs*(a13a.v(ds1)*s24s.v(ds4)-a34a.v(ds4)*s12s.v(ds1))*s213a.v();
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble dggd::amp2(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_Num = amp_Num(j1,j2,j3,j4);
	    //1<->3:  S<->U
	    amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenU*pDenT,2)+6.0/std::pow(pDenU*pDenS,2)-2.0/3.0/std::pow(pDenT*pDenS,2));
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/3=1/24
    return amp2/96.0;
  }

  //D, g+ -> g, D
  ldouble dggd::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M_Num = amp_Num(j1,2,j3,j4);
	  amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenU*pDenT,2)+6.0/std::pow(pDenU*pDenS,2)-2.0/3.0/std::pow(pDenT*pDenS,2));
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3*1/8=1/24
    return amp2/48.0;
  }

  //D, g- -> g, D
  ldouble dggd::amp2_gminus(){
    ldouble amp2 = 0;
    cdouble M_Num;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M_Num = amp_Num(j1,-2,j3,j4);
	  amp2 += std::pow(std::abs(M_Num),2)*(6.0/std::pow(pDenU*pDenT,2)+6.0/std::pow(pDenU*pDenS,2)-2.0/3.0/std::pow(pDenT*pDenS,2));
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3*1/8=1/24
    return amp2/48.0;
  }
  



  //  Tests
  int test_dggd(){
    int n=0;//Number of fails
    std::cout<<"\t* D , g  -> g , D       :";
    {//amp^2
      int i=0;
      // md=0.0075, pspatial=250
      ldouble md=0.0075;
      ldouble gg=1.238;
      dggd dggdAmp = dggd(gg,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.425860525151625E+01,1.675908728203342E+01,1.159850868919754E+01,9.705338369637635E+00,8.983807684842580E+00,8.890385698846716E+00,9.251709468536607E+00,1.003456987673026E+01,1.128817621489624E+01,1.313910801515463E+01,1.581728982539627E+01,1.972050322087015E+01,2.555187423354871E+01,3.462305364810209E+01,4.958452768963280E+01,7.642511971706138E+01,1.310339897785782E+02,2.675433731152843E+02,7.770025652758417E+02,7.333310275906781E+03};
      i += dggdAmp.test_2to2_amp2([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dggdAmp.test_2to2_amp2_rotations([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dggdAmp.test_2to2_amp2_boosts([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      i += dggdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH);
      //Helicities D, g+ -> g, D
      ldouble dataCHpp[20] = {4.425860525151625E+01,1.675908728203342E+01,1.159850868919754E+01,9.705338369637635E+00,8.983807684842580E+00,8.890385698846716E+00,9.251709468536607E+00,1.003456987673026E+01,1.128817621489624E+01,1.313910801515463E+01,1.581728982539627E+01,1.972050322087015E+01,2.555187423354871E+01,3.462305364810209E+01,4.958452768963280E+01,7.642511971706138E+01,1.310339897785782E+02,2.675433731152843E+02,7.770025652758417E+02,7.333310275906781E+03};
      i += dggdAmp.test_2to2_amp2([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dggdAmp.test_2to2_amp2_rotations([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dggdAmp.test_2to2_amp2_boosts([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCHpp);
      i += dggdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCHpp);
      //Helicities D, g- -> g, D
      ldouble dataCHpm[20] = {4.425860525151625E+01,1.675908728203342E+01,1.159850868919754E+01,9.705338369637635E+00,8.983807684842580E+00,8.890385698846716E+00,9.251709468536607E+00,1.003456987673026E+01,1.128817621489624E+01,1.313910801515463E+01,1.581728982539627E+01,1.972050322087015E+01,2.555187423354871E+01,3.462305364810209E+01,4.958452768963280E+01,7.642511971706138E+01,1.310339897785782E+02,2.675433731152843E+02,7.770025652758417E+02,7.333310275906781E+03};
      i += dggdAmp.test_2to2_amp2([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dggdAmp.test_2to2_amp2_rotations([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dggdAmp.test_2to2_amp2_boosts([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCHpm);
      i += dggdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCHpm);
      //Close to threshold
      pspatial = 0.008;
      ldouble dataCH2[20] = {8.845631130707270E+00,7.324785496487646E+00,6.698819763549957E+00,6.588326469869709E+00,6.849394210449703E+00,7.440472681671003E+00,8.379639758381193E+00,9.734452029698776E+00,1.162929103175043E+01,1.426924511295844E+01,1.798874831170677E+01,2.334674011419745E+01,3.132017149287364E+01,4.372642740050046E+01,6.423951915774859E+01,1.011747020946623E+02,1.766400229443928E+02,3.661260517494143E+02,1.076488046517356E+03,1.026120413069680E+04};
      i += dggdAmp.test_2to2_amp2([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dggdAmp.test_2to2_amp2_rotations([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dggdAmp.test_2to2_amp2_boosts([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      i += dggdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dggdAmp.amp2(); }, md,0,0,md,pspatial,dataCH2);
      //Helicities D, g+ -> g, D
      ldouble dataCH2pp[20] = {8.845631130707270E+00,7.324785496487646E+00,6.698819763549957E+00,6.588326469869709E+00,6.849394210449703E+00,7.440472681671003E+00,8.379639758381193E+00,9.734452029698776E+00,1.162929103175043E+01,1.426924511295844E+01,1.798874831170677E+01,2.334674011419745E+01,3.132017149287364E+01,4.372642740050046E+01,6.423951915774859E+01,1.011747020946623E+02,1.766400229443928E+02,3.661260517494143E+02,1.076488046517356E+03,1.026120413069680E+04};
      i += dggdAmp.test_2to2_amp2([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dggdAmp.test_2to2_amp2_rotations([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dggdAmp.test_2to2_amp2_boosts([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCH2pp);
      i += dggdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dggdAmp.amp2_gplus(); }, md,0,0,md,pspatial,dataCH2pp);
      //Helicities D, g- -> g, D
      ldouble dataCH2pm[20] = {8.845631130707270E+00,7.324785496487646E+00,6.698819763549957E+00,6.588326469869709E+00,6.849394210449703E+00,7.440472681671003E+00,8.379639758381193E+00,9.734452029698776E+00,1.162929103175043E+01,1.426924511295844E+01,1.798874831170677E+01,2.334674011419745E+01,3.132017149287364E+01,4.372642740050046E+01,6.423951915774859E+01,1.011747020946623E+02,1.766400229443928E+02,3.661260517494143E+02,1.076488046517356E+03,1.026120413069680E+04};
      i += dggdAmp.test_2to2_amp2([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dggdAmp.test_2to2_amp2_rotations([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dggdAmp.test_2to2_amp2_boosts([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCH2pm);
      i += dggdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dggdAmp.amp2_gminus(); }, md,0,0,md,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



  
  

}


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

//File:  SPINAS/SM/gggg.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gggg.h"

namespace spinas {

  gggg::gggg(const ldouble& gcharge):
    gs(gcharge), prop(){
    p1=particle(0);
    p2=particle(0);
    p3=particle(0);
    p4=particle(0);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s41s = sproduct(SQUARE,&p4,&p1,2);
    a41a = sproduct(ANGLE,&p4,&p1,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s1432s = sproduct(SQUARE,&p1,&p4,&p3,&p2,2);
    a1432a = sproduct(ANGLE,&p1,&p4,&p3,&p2,2);
    s3214s = sproduct(SQUARE,&p3,&p2,&p1,&p4,2);
    a3214a = sproduct(ANGLE,&p3,&p2,&p1,&p4,2);
    a3124a = sproduct(ANGLE,&p3,&p1,&p2,&p4,2);
  }
  void gggg::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s41s.update();
    a41a.update();
    s1432s.update();
    a1432a.update();
    s3214s.update();
    a3214a.update();
    a3124a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS = std::real(prop.denominator(propSP));
    pDenT = std::real(prop.denominator(propTP));
    pDenU = std::real(prop.denominator(propUP));
  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  ////////////////////////////////////////////////////
  //This amplitude is not yet in agreement with CH
  //It is an attempt to test the BCFW amps.
  //They seem to be wrong.
  //My guess is much closer but still wrong.
  //I also do not have the color factor right.
  ////////////////////////////////////////////////////
  cdouble gggg::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr cdouble zero = cdouble(0,0), one = cdouble(1,0);
    constexpr ldouble two=2, oneHalf = 1.0/2.0;
    //constexpr ldouble sqrt2 = std::sqrt(two);

    //++++
    if(ds1>0&&ds2>0&&ds3>0&&ds4>0)
      return zero;
    //+++- & ++-+ & +-++ & -+++
    else if((ds1>0&&ds2>0&&ds3>0&&ds4<0) ||
	    (ds1>0&&ds2>0&&ds3<0&&ds4>0) ||
	    (ds1>0&&ds2<0&&ds3>0&&ds4>0) ||
	    (ds1<0&&ds2>0&&ds3>0&&ds4>0))
      return zero;
    //++--
    else if(ds1>0&&ds2>0&&ds3<0&&ds4<0){
      //Square: s^4/den^2
      //all in:   [12]^4/[13][23][24][41]
      //34 out:   [12]^4/[13][23][24][41]
      //return sqrt2*gs*gs*s12s.v()*s12s.v()*a34a.v()*a34a.v()/pDenT/pDenU;
      //The denominator has to be done separately because it has a different color factor for each denominator combination.
      return two*gs*gs*s12s.v()*s12s.v()*a34a.v()*a34a.v();
    }
    //+-+-
    else if(ds1>0&&ds2<0&&ds3>0&&ds4<0){
      //Square: t^4/den^2
      //all in:   [13]^4/[12][23][34][41]
      //34 out:   [13]^4/[12][23][34][41]
      return two*gs*gs*s13s.v()*s13s.v()*a24a.v()*a24a.v();
    }
    //+--+
    else if(ds1>0&&ds2<0&&ds3<0&&ds4>0){
      //Square: u^4/den^2
      //all in:   [14]^4/[12][23][34][41]
      //34 out:   [14]^4/[12][23][34][41]
      return two*gs*gs*s41s.v()*s41s.v()*a23a.v()*a23a.v();
    }
    //-++-
    else if(ds1<0&&ds2>0&&ds3>0&&ds4<0){
      //Square: u^4/den^2
      //all in:   [23]^4/[12][23][34][41]
      //34 out:   [23]^4/[12][23][34][41]
      return two*gs*gs*s23s.v()*s23s.v()*a41a.v()*a41a.v();
    }
    //-+-+
    else if(ds1<0&&ds2>0&&ds3<0&&ds4>0){
      //Square: t^4/den^2
      //all in:   [24]^4/[12][23][34][41]
      //34 out:   [24]^4/[12][23][34][41]
      return two*gs*gs*s24s.v()*s24s.v()*a13a.v()*a13a.v();
    }
    //--++
    else if(ds1<0&&ds2<0&&ds3>0&&ds4>0){
      //Square: s^4/den^2
      //all in:   [34]^4/[12][23][34][41]
      //34 out:   [34]^4/[12][23][34][41]
      return two*gs*gs*s34s.v()*s34s.v()*a12a.v()*a12a.v();
    }
    //+--- & -+-- & --+- & ---+
    else if((ds1>0&&ds2<0&&ds3<0&&ds4<0) ||
	    (ds1<0&&ds2>0&&ds3<0&&ds4<0) ||
	    (ds1<0&&ds2<0&&ds3>0&&ds4<0) ||
	    (ds1<0&&ds2<0&&ds3<0&&ds4>0))
      return zero;
    //----
    else if(ds1<0&&ds2<0&&ds3<0&&ds4<0)
      return zero;
    throw std::runtime_error("Gluon helicity combination not caught.  Double helicities =  " + std::to_string(ds1) + " , " + std::to_string(ds2) + " , " + std::to_string(ds3) + " , " + std::to_string(ds4));
  }

  
  //set_momenta(...) must be called before amp2_Feynman().
  //This expression comes from Ellis, Stirling and Weber: QCD and Collider Physics, Table 7.1
  //This agrees with CalcHEP.
  ldouble gggg::amp2_Feynman(){
    ldouble pre = 9.0/4.0*gs*gs*gs*gs;
    constexpr ldouble three = 3.0;
    cdouble M2 =  pre*(three - pDenT*pDenU/pDenS/pDenS - pDenS*pDenU/pDenT/pDenT - pDenS*pDenT/pDenU/pDenU);
    return std::abs(M2);
  }
 
  //set_momenta(...) must be called before amp2().
  ////////////////////////////////////////////////////
  //Not confident of color factor
  //Missing symmetry factor
  ////////////////////////////////////////////////////
  ldouble gggg::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2)*
	  (36.0/std::pow(pDenS*pDenT,2)+36.0/std::pow(pDenS*pDenU,2)+36.0/std::pow(pDenT*pDenU,2));;
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/8=1/64
    //Symmetry factor 1/2
    return amp2/512.0;
  }
  
  ldouble gggg::amp2_gplus_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=4)
      for(int j4=-2;j4<=2;j4+=4){
	M = amp(2,2,j3,j4);
	amp2 += std::pow(std::abs(M),2)*
	  (36.0/std::pow(pDenS*pDenT,2)+36.0/std::pow(pDenS*pDenU,2)+36.0/std::pow(pDenT*pDenU,2));
      }
    //Average over initial colors 1/8*1/8=1/64
    //Symmetry factor 1/2
    return amp2/128.0;
  }
  
  ldouble gggg::amp2_gplus_gminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=4)
      for(int j4=-2;j4<=2;j4+=4){
	M = amp(2,-2,j3,j4);
	amp2 += std::pow(std::abs(M),2)*
	  (36.0/std::pow(pDenS*pDenT,2)+36.0/std::pow(pDenS*pDenU,2)+36.0/std::pow(pDenT*pDenU,2));
      }
    //Average over initial colors 1/8*1/8=1/64
    //Symmetry factor 1/2
    return amp2/128.0;
  }
  
  ldouble gggg::amp2_gminus_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=4)
      for(int j4=-2;j4<=2;j4+=4){
	M = amp(-2,2,j3,j4);
	amp2 += std::pow(std::abs(M),2)*
	  (36.0/std::pow(pDenS*pDenT,2)+36.0/std::pow(pDenS*pDenU,2)+36.0/std::pow(pDenT*pDenU,2));;
      }
    //Average over initial colors 1/8*1/8=1/64
    //Symmetry factor 1/2
    return amp2/128.0;
  }
  
  ldouble gggg::amp2_gminus_gminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-2;j3<=2;j3+=4)
      for(int j4=-2;j4<=2;j4+=4){
	M = amp(-2,-2,j3,j4);
	amp2 += std::pow(std::abs(M),2)*
	  (36.0/std::pow(pDenS*pDenT,2)+36.0/std::pow(pDenS*pDenU,2)+36.0/std::pow(pDenT*pDenU,2));;
      }
    //Average over initial colors 1/8*1/8=1/64
    //Symmetry factor 1/2
    return amp2/128.0;
  }


  //  Tests
  int test_gggg(){
    int n=0;//Number of fails
    std::cout<<"\t* g , g  -> g , g       :";
    {//amp^2
      int i=0;
      // pspatial=250
      ldouble gs=1.238;
      gggg ggggAmp = gggg(gs);
      ldouble pspatial=250;
      //std::cout<<"\nHelicity Sum\n";
      ldouble dataCH[20] = {8.260847583994650E+03,8.850814728298411E+02,3.121142589333567E+02,1.588296207164410E+02,9.782394221087002E+01,6.823561926740386E+01,5.224180771889487E+01,4.318081597788433E+01,3.818303995970413E+01,3.594416051657328E+01,3.594416051657328E+01,3.818303995970411E+01,4.318081597788431E+01,5.224180771889485E+01,6.823561926740385E+01,9.782394221087000E+01,1.588296207164409E+02,3.121142589333566E+02,8.850814728298401E+02,8.260847583994611E+03};
      i += ggggAmp.test_2to2_amp2([&]() { return ggggAmp.amp2_Feynman(); }, 0,0,0,0,pspatial,dataCH);
      i += ggggAmp.test_2to2_amp2([&]() { return ggggAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += ggggAmp.test_2to2_amp2_rotations([&]() { return ggggAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += ggggAmp.test_2to2_amp2_boosts([&]() { return ggggAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += ggggAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggggAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      //Helicity + +
      //std::cout<<"\n + + Helicities\n";
      ldouble dataCH2[20] = {8.678781778884946E+03,1.021959816556874E+03,3.934810725118587E+02,2.169524544764772E+02,1.435091164389990E+02,1.064517592584796E+02,8.573009104505115E+01,7.366456569988097E+01,6.687408695662225E+01,6.379436156965200E+01,6.379436156965199E+01,6.687408695662222E+01,7.366456569988094E+01,8.573009104505115E+01,1.064517592584795E+02,1.435091164389990E+02,2.169524544764770E+02,3.934810725118587E+02,1.021959816556873E+03,8.678781778884906E+03};
      i += ggggAmp.test_2to2_amp2([&]() { return ggggAmp.amp2_gplus_gplus(); }, 0,0,0,0,pspatial,dataCH2);
      i += ggggAmp.test_2to2_amp2_rotations([&]() { return ggggAmp.amp2_gplus_gplus(); }, 0,0,0,0,pspatial,dataCH2);
      i += ggggAmp.test_2to2_amp2_boosts([&]() { return ggggAmp.amp2_gplus_gplus(); }, 0,0,0,0,pspatial,dataCH2);
      i += ggggAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggggAmp.amp2_gplus_gplus(); }, 0,0,0,0,pspatial,dataCH2);
      //Helicity + -
      //std::cout<<"\n + - Helicities\n";
      ldouble dataCH3[20] = {7.842913389104354E+03,7.482031291028081E+02,2.307474453548546E+02,1.007067869564048E+02,5.213876798274104E+01,3.001947927632816E+01,1.875352439273858E+01,1.269706625588769E+01,9.491992962786007E+00,8.093959463494574E+00,8.093959463494572E+00,9.491992962786000E+00,1.269706625588768E+01,1.875352439273857E+01,3.001947927632814E+01,5.213876798274102E+01,1.007067869564047E+02,2.307474453548546E+02,7.482031291028072E+02,7.842913389104318E+03};
      i += ggggAmp.test_2to2_amp2([&]() { return ggggAmp.amp2_gplus_gminus(); }, 0,0,0,0,pspatial,dataCH3);
      i += ggggAmp.test_2to2_amp2_rotations([&]() { return ggggAmp.amp2_gplus_gminus(); }, 0,0,0,0,pspatial,dataCH3);
      i += ggggAmp.test_2to2_amp2_boosts([&]() { return ggggAmp.amp2_gplus_gminus(); }, 0,0,0,0,pspatial,dataCH3);
      i += ggggAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggggAmp.amp2_gplus_gminus(); }, 0,0,0,0,pspatial,dataCH3);
      //Helicity - +
      //std::cout<<"\n - + Helicities\n";
      ldouble dataCH4[20] = {7.842913389104354E+03,7.482031291028081E+02,2.307474453548546E+02,1.007067869564048E+02,5.213876798274104E+01,3.001947927632816E+01,1.875352439273858E+01,1.269706625588769E+01,9.491992962786007E+00,8.093959463494574E+00,8.093959463494572E+00,9.491992962786000E+00,1.269706625588768E+01,1.875352439273857E+01,3.001947927632814E+01,5.213876798274102E+01,1.007067869564047E+02,2.307474453548546E+02,7.482031291028072E+02,7.842913389104318E+03};
      i += ggggAmp.test_2to2_amp2([&]() { return ggggAmp.amp2_gminus_gplus(); }, 0,0,0,0,pspatial,dataCH4);
      i += ggggAmp.test_2to2_amp2_rotations([&]() { return ggggAmp.amp2_gminus_gplus(); }, 0,0,0,0,pspatial,dataCH4);
      i += ggggAmp.test_2to2_amp2_boosts([&]() { return ggggAmp.amp2_gminus_gplus(); }, 0,0,0,0,pspatial,dataCH4);
      i += ggggAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggggAmp.amp2_gminus_gplus(); }, 0,0,0,0,pspatial,dataCH4);
      //Helicity - -
      //std::cout<<"\n - - Helicities\n";
      ldouble dataCH5[20] = {8.678781778884946E+03,1.021959816556874E+03,3.934810725118587E+02,2.169524544764772E+02,1.435091164389990E+02,1.064517592584796E+02,8.573009104505115E+01,7.366456569988097E+01,6.687408695662225E+01,6.379436156965200E+01,6.379436156965199E+01,6.687408695662222E+01,7.366456569988094E+01,8.573009104505115E+01,1.064517592584795E+02,1.435091164389990E+02,2.169524544764770E+02,3.934810725118587E+02,1.021959816556873E+03,8.678781778884906E+03};
      i += ggggAmp.test_2to2_amp2([&]() { return ggggAmp.amp2_gminus_gminus(); }, 0,0,0,0,pspatial,dataCH5);
      i += ggggAmp.test_2to2_amp2_rotations([&]() { return ggggAmp.amp2_gminus_gminus(); }, 0,0,0,0,pspatial,dataCH5);
      i += ggggAmp.test_2to2_amp2_boosts([&]() { return ggggAmp.amp2_gminus_gminus(); }, 0,0,0,0,pspatial,dataCH5);
      i += ggggAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ggggAmp.amp2_gminus_gminus(); }, 0,0,0,0,pspatial,dataCH5);
      
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  


}


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

//File:  SPINAS/SM/ulnd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ulnd.h"

namespace spinas {
  //Constructors
  ulnd::ulnd(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& masstau, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mu(massu), md(massd), ml(masstau), MW(massW), WW(widthW), SW(sinW), prop(massW,widthW),
    p1(particle(mu)), p2(particle(ml)),
    p3(particle(0)), p4(particle(md)),
    //[14], <23>, <34>, [12], <12>
    a14a(sproduct(ANGLE,&p1,&p4)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    a34a(sproduct(ANGLE,&p3,&p4))
  {
    preW = e*e/(4.0*MW*MW*SW*SW);
  }
  void ulnd::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& masstau, const ldouble& massW){
    mu=massu;
    md=massd;
    ml=masstau;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(ml);
    p4.set_mass(md);
    prop.set_mass(massW);
    preW = e*e/(4.0*MW*MW*SW*SW);
  }
  void ulnd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a14a.update();
    s14s.update();
    a23a.update();
    s12s.update();
    a12a.update();
    a34a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++){
      propP[j] = mom1[j]-mom4[j];
    }
    pDenU = prop.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ulnd::amp(const int& ds1, const int& ds2, const int& ds4){
    //preW = e*e/(4.0*MW*MW*SW*SW);
    //uDnlL all ingoing:
    //+ preW ( 2MW^2 [14]<23> + mlmu<34><12> - mlmd<34>[12] )/(s−MW^2)
    //uLNlD: 2<->4
    //+ preW (- 2MW^2 [12]<34> - mlmu<23><14> + mlmd<23>[14] )/(u−MW^2)
    //34 out:
    //- preW ( 2MW^2 [12]<34> + mlmu<23><14> + mlmd<23>[14] )/(u−MW^2)
    return - preW*2.0*(
		       + 2.0*MW*MW*s12s.v(ds1,ds2)*a34a.v(ds4)
		       + mu*ml*a23a.v(ds2)*a14a.v(ds1,ds4)
		       + ml*md*a23a.v(ds2)*s14s.v(ds1,ds4)
		       )/pDenU;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ulnd::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble three=3;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,j2,j4);
	  amp2 += three*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over initial colors 1/3
    return amp2/12.0;
  }
  



  //  Tests
  int test_ulnd(){
    int n=0;//Number of fails
    std::cout<<"\t* u , l  -> nl, d       :";
    {//amp^2
      int i=0;
      // mu=0.0042, md=0.0075, mtau=1.777, pspatial=250
      ldouble mu=0.0042, md=0.0075, mtau=1.777, MW=80.385, WW=0;
      ldouble sinW=0.474;
      ulnd ulndAmp = ulnd(0.31333,mu,md,mtau,MW,WW,sinW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.765461314731104E-02,5.279819380806557E-02,5.882179737072814E-02,6.593825337503889E-02,7.442884816403622E-02,8.467155358388845E-02,9.718383175801965E-02,1.126892978230035E-01,1.322247426233121E-01,1.573180632789198E-01,1.902963539283075E-01,2.348455399025207E-01,2.970868916602365E-01,3.877977932801755E-01,5.274115484947456E-01,7.586184039743973E-01,1.183343382669650E+00,2.097821867031888E+00,4.693717514737354E+00,1.846354424791671E+01};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH);
      // mu=0.0042, md=0.0075, mtau=1.777, pspatial=1.8
      pspatial = 1.8;
      ldouble dataCH2[20] = {3.324085538994853E-07,3.324885570333439E-07,3.325685890530723E-07,3.326486499725781E-07,3.327287398057775E-07,3.328088585665949E-07,3.328890062689630E-07,3.329691829268232E-07,3.330493885541250E-07,3.331296231648264E-07,3.332098867728937E-07,3.332901793923020E-07,3.333705010370343E-07,3.334508517210825E-07,3.335312314584464E-07,3.336116402631348E-07,3.336920781491646E-07,3.337725451305611E-07,3.338530412213585E-07,3.339335664355989E-07};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH2);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH2);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH2);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH2);
      // mu=0.0042, md=0.0075, mtau=0.005, pspatial=250
      mu=0.0042;
      md=0.0075;
      mtau=0.005;
      pspatial = 250;
      ulndAmp.set_masses(mu,md,mtau,MW);
      ldouble dataCH3[20] = {4.765398015784574E-02,5.279749068595890E-02,5.882101178950035E-02,6.593736994319352E-02,7.442784741048661E-02,8.467041051237466E-02,9.718251374095767E-02,1.126877614542187E-01,1.322229288766073E-01,1.573158898063454E-01,1.902937023173211E-01,2.348422335855611E-01,2.970826554006200E-01,3.877921735251572E-01,5.274037424117062E-01,7.586068476922729E-01,1.183324589650575E+00,2.097786290486642E+00,4.693627842322704E+00,1.846307395956124E+01};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH3);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH3);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH3);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH3);
      // mu=0.0042, md=0.0075, mtau=0.005, pspatial=0.001
      pspatial = 0.001;
      ldouble dataCH4[20] = {1.705954733390640E-18,1.705954733486770E-18,1.705954733582900E-18,1.705954733679030E-18,1.705954733775160E-18,1.705954733871290E-18,1.705954733967420E-18,1.705954734063550E-18,1.705954734159680E-18,1.705954734255809E-18,1.705954734351940E-18,1.705954734448069E-18,1.705954734544199E-18,1.705954734640329E-18,1.705954734736459E-18,1.705954734832589E-18,1.705954734928719E-18,1.705954735024849E-18,1.705954735120979E-18,1.705954735217109E-18};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH4);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH4);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH4);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH4);
      // mu=0.0042, md=0.0075, mtau=0.005, MW=0.006, pspatial=250
      MW=0.006;
      pspatial=250;
      ulndAmp.set_masses(mu,md,mtau,MW);
      ldouble dataCH5[20] = {6.722370322670571E-02,7.279896758897715E-02,7.935705943456038E-02,8.714332979051814E-02,9.648475716221772E-02,1.078248445664051E-01,1.217772484008065E-01,1.392105776709527E-01,1.613869284701800E-01,1.901969080117866E-01,2.285763569572835E-01,2.812850188429024E-01,3.564566808215167E-01,4.689361720629358E-01,6.482128109572346E-01,9.599180780768048E-01,1.575695013655833E+00,3.072032965935234E+00,8.503185524140635E+00,7.639259197357448E+01};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH5);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH5);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH5);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH5);
      // mu=0.0042, md=0.0075, mtau=0.005, MW=0.006, pspatial=0.001
      pspatial=0.001;
      ldouble dataCH6[20] = {7.329474958194877E-02,7.404331839995912E-02,7.481141541286615E-02,7.559979651276250E-02,7.640925639353401E-02,7.724063103647769E-02,7.809480038718851E-02,7.897269124093924E-02,7.987528035555676E-02,8.080359781278854E-02,8.175873065137652E-02,8.274182679754785E-02,8.375409932142491E-02,8.479683105099579E-02,8.587137957881455E-02,8.697918270057475E-02,8.812176432918375E-02,8.930074093302848E-02,9.051782855285544E-02,9.177485045818114E-02};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH6);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH6);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH6);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH6);
      // mu=0.0042, md=0.0075, mtau=0.005, MW=0.0006, pspatial=0.001
      MW=0.0006;
      ulndAmp.set_masses(mu,md,mtau,MW);
      ldouble dataCH7[20] = {3.882478896695871E+03,3.138991353553574E+03,2.570790374958549E+03,2.128332720677232E+03,1.778237300210656E+03,1.497378656771067E+03,1.269346729849616E+03,1.082250159774909E+03,9.273115715532997E+02,7.979452260786248E+02,6.891370465992479E+02,5.970190970589757E+02,5.185719974025496E+02,4.514132597626629E+02,3.936444095636587E+02,3.437390102051986E+02,3.004595909340115E+02,2.627952885535734E+02,2.299145279854658E+02,2.011287531268409E+02};
      i += ulndAmp.test_2to2_amp2([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH7);
      i += ulndAmp.test_2to2_amp2_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH7);
      i += ulndAmp.test_2to2_amp2_boosts([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH7);
      i += ulndAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ulndAmp.amp2(); }, mu,mtau,0,md,pspatial,dataCH7);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}

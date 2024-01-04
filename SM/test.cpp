
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

//File:  SPINAS/SM/test.cpp

#include <cstdlib>
#include<iostream>
#include <sstream>
#include <complex>
using namespace std;

#include "spinas.h"
#include "include/eemmQED.h"
#include "include/eeeeQED.h"
#include "include/uucc.h"
#include "include/ucuc.h"
#include "include/uuuu.h"
#include "include/uuuu2.h"
#include "include/uuss.h"
#include "include/usus.h"
#include "include/uudd.h"
#include "include/udud.h"
#include "include/uuee.h"
#include "include/ueeu.h"
#include "include/uunn.h"
#include "include/unnu.h"
#include "include/ddss.h"
#include "include/dsds.h"
#include "include/dddd.h"
#include "include/dddd2.h"
#include "include/ddee.h"
#include "include/deed.h"
#include "include/ddnn.h"
#include "include/dnnd.h"
#include "include/eemm.h"
#include "include/emem.h"
#include "include/eeee.h"
#include "include/eeee2.h"
#include "include/eenmnm.h"
#include "include/enmenm.h"
#include "include/eenene.h"
#include "include/enenee.h"
#include "include/nenenmnm.h"
#include "include/nenmnenm.h"
#include "include/nenenene.h"
#include "include/nenenene2.h"
#include "include/udtb.h"
#include "include/ubdt.h"
#include "include/udnl.h"
#include "include/ulnd.h"
#include "include/menn.h"
#include "include/mnen.h"
#include "include/eehh.h"
#include "include/eheh.h"
#include "include/eeAh.h"
#include "include/eAeh.h"
#include "include/eeZh.h"
#include "include/eZeh.h"
#include "include/eeAA.h"
#include "include/eAAe.h"
#include "include/AAee.h"
#include "include/AZee.h"
#include "include/AeZe.h"
#include "include/eeZZ.h"
#include "include/eZZe.h"
#include "include/eeWW.h"
#include "include/eWWe.h"
#include "include/neneZh.h"
#include "include/neZneh.h"
#include "include/neneZZ.h"
#include "include/neZZne.h"
#include "include/neneWW.h"
#include "include/neWWne.h"
#include "include/uuhh.h"
#include "include/uhuh.h"
#include "include/uuAh.h"
#include "include/uAuh.h"
#include "include/uuZh.h"
#include "include/uZuh.h"
#include "include/uuAA.h"
#include "include/AAuu.h"
#include "include/uAAu.h"
#include "include/AZuu.h"
#include "include/AuZu.h"
#include "include/uuZZ.h"
#include "include/uZZu.h"
#include "include/uuWW.h"
#include "include/uWWu.h"
#include "include/uguh.h"
#include "include/hguu.h"
#include "include/gZuu.h"
#include "include/guZu.h"
#include "include/gAuu.h"
#include "include/guAu.h"
#include "include/gguu.h"
#include "include/gugu.h"
#include "include/ddhh.h"
#include "include/dhdh.h"
#include "include/ddAh.h"
#include "include/dAdh.h"
#include "include/ddZh.h"
#include "include/dZdh.h"
#include "include/ddAA.h"
#include "include/AAdd.h"
#include "include/dAAd.h"
#include "include/AZdd.h"
#include "include/AdZd.h"
#include "include/ddZZ.h"
#include "include/dZZd.h"
#include "include/ddWW.h"
#include "include/WdWd.h"
#include "include/dgdh.h"
#include "include/hgdd.h"
#include "include/gZdd.h"
#include "include/gdZd.h"
#include "include/gAdd.h"
#include "include/dAgd.h"
#include "include/ggdd.h"
#include "include/dggd.h"
#include "include/enWh.h"
#include "include/hnWe.h"
#include "include/enAW.h"
#include "include/AnWe.h"
#include "include/enZW.h"
#include "include/ZnWe.h"
#include "include/udWh.h"
#include "include/uhWd.h"
#include "include/udAW.h"
#include "include/AWud.h"
#include "include/gWud.h"
#include "include/hhhh.h"
#include "include/hhZZ.h"
#include "include/hZZh.h"
#include "include/hhWW.h"
#include "include/hWWh.h"
#include "include/AhWW.h"
#include "include/AWWh.h"
#include "include/ZhWW.h"
#include "include/ZWWh.h"
#include "include/AAWW.h"
#include "include/AWAW.h"
#include "include/AZWW.h"
#include "include/AWZW.h"
#include "include/ZZZZ.h"
#include "include/ZZWW.h"
#include "include/ZWZW.h"
#include "include/WWWW.h"
#include "include/WWWW2.h"
#include "include/gggg.h"



int main(){
 
  
  stringstream message;
  cout<<"sizeof(float)="<<sizeof(float)<<endl;
  cout<<"sizeof(double)="<<sizeof(double)<<endl;
  cout<<"sizeof(long double)="<<sizeof(long double)<<endl;
  cout<<"sizeof(ldouble)="<<sizeof(ldouble)<<endl;

  int i=0;
  int j=0;
  int n=0;

  cout<<"==============================================================================\n";
  cout<<"Testing SM Processes\n";
  cout<<"=============================================================================="<<endl;
  

  j = spinas::test_eemmQED(); if(j>0){message<<"    eemmQED: "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeeeQED(); if(j>0){message<<"    eeeeQED: "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uucc(); if(j>0){message<<"    uucc   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ucuc(); if(j>0){message<<"    ucuc   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuuu(); if(j>0){message<<"    uuuu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuuu2(); if(j>0){message<<"    uuuu2  : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuss(); if(j>0){message<<"    uuss   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_usus(); if(j>0){message<<"    usus   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uudd(); if(j>0){message<<"    uudd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_udud(); if(j>0){message<<"    udud   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuee(); if(j>0){message<<"    uuee   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ueeu(); if(j>0){message<<"    ueeu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uunn(); if(j>0){message<<"    uunn   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_unnu(); if(j>0){message<<"    unnu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddss(); if(j>0){message<<"    ddss   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dsds(); if(j>0){message<<"    dsds   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dddd(); if(j>0){message<<"    dddd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dddd2(); if(j>0){message<<"    dddd2   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddee(); if(j>0){message<<"    ddee   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_deed(); if(j>0){message<<"    deed   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddnn(); if(j>0){message<<"    ddnn   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dnnd(); if(j>0){message<<"    dnnd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eemm(); if(j>0){message<<"    eemm   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_emem(); if(j>0){message<<"    emem   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeee(); if(j>0){message<<"    eeee   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeee2(); if(j>0){message<<"    eeee2  : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eenmnm(); if(j>0){message<<"    eenmnm   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_enmenm(); if(j>0){message<<"    enmenm   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eenene(); if(j>0){message<<"    eenene   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_enenee(); if(j>0){message<<"    enenee   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_nenenmnm(); if(j>0){message<<"    nenenmnm   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_nenmnenm(); if(j>0){message<<"    nenmnenm   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_nenenene(); if(j>0){message<<"    nenenene   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_nenenene2(); if(j>0){message<<"    nenenene2   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_udtb(); if(j>0){message<<"    udtb   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ubdt(); if(j>0){message<<"    ubdt   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_udnl(); if(j>0){message<<"    udnl   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ulnd(); if(j>0){message<<"    ulnd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_menn(); if(j>0){message<<"    menn   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_mnen(); if(j>0){message<<"    mnen   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eehh(); if(j>0){message<<"    eehh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eheh(); if(j>0){message<<"    eheh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeAh(); if(j>0){message<<"    eeAh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eAeh(); if(j>0){message<<"    eAeh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeZh(); if(j>0){message<<"    eeZh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eZeh(); if(j>0){message<<"    eZeh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeAA(); if(j>0){message<<"    eeAA   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eAAe(); if(j>0){message<<"    eAAe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AAee(); if(j>0){message<<"    AAee   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AZee(); if(j>0){message<<"    AZee   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AeZe(); if(j>0){message<<"    AeZe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeZZ(); if(j>0){message<<"    eeZZ   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eZZe(); if(j>0){message<<"    eZZe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eeWW(); if(j>0){message<<"    eeWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_eWWe(); if(j>0){message<<"    eWWe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_neneZh(); if(j>0){message<<"    neneZh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_neZneh(); if(j>0){message<<"    neZneh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_neneZZ(); if(j>0){message<<"    neneZZ   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_neZZne(); if(j>0){message<<"    neZZne   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_neneWW(); if(j>0){message<<"    neneWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_neWWne(); if(j>0){message<<"    neWWne   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuhh(); if(j>0){message<<"    uuhh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uhuh(); if(j>0){message<<"    uhuh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuAh(); if(j>0){message<<"    uuAh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uAuh(); if(j>0){message<<"    uAuh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuZh(); if(j>0){message<<"    uuZh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uZuh(); if(j>0){message<<"    uZuh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuAA(); if(j>0){message<<"    uuAA   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AAuu(); if(j>0){message<<"    AAuu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uAAu(); if(j>0){message<<"    uAAu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AZuu(); if(j>0){message<<"    AZuu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AuZu(); if(j>0){message<<"    AuZu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuZZ(); if(j>0){message<<"    uuZZ   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uZZu(); if(j>0){message<<"    uZZu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uuWW(); if(j>0){message<<"    uuWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uWWu(); if(j>0){message<<"    uWWu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uguh(); if(j>0){message<<"    uguh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hguu(); if(j>0){message<<"    hguu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gZuu(); if(j>0){message<<"    gZuu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_guZu(); if(j>0){message<<"    guZu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gAuu(); if(j>0){message<<"    gAuu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_guAu(); if(j>0){message<<"    guAu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gguu(); if(j>0){message<<"    gguu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gugu(); if(j>0){message<<"    gugu   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddhh(); if(j>0){message<<"    ddhh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dhdh(); if(j>0){message<<"    dhdh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddAh(); if(j>0){message<<"    ddAh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dAdh(); if(j>0){message<<"    dAdh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddZh(); if(j>0){message<<"    ddZh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dZdh(); if(j>0){message<<"    dZdh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddAA(); if(j>0){message<<"    ddAA   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AAdd(); if(j>0){message<<"    AAdd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dAAd(); if(j>0){message<<"    dAAd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AZdd(); if(j>0){message<<"    AZdd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AdZd(); if(j>0){message<<"    AdZd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddZZ(); if(j>0){message<<"    ddZZ   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dZZd(); if(j>0){message<<"    dZZd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ddWW(); if(j>0){message<<"    ddWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_WdWd(); if(j>0){message<<"    WdWd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dgdh(); if(j>0){message<<"    dgdh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hgdd(); if(j>0){message<<"    hgdd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gZdd(); if(j>0){message<<"    gZdd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gdZd(); if(j>0){message<<"    gdZd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gAdd(); if(j>0){message<<"    gAdd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dAgd(); if(j>0){message<<"    dAgd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ggdd(); if(j>0){message<<"    ggdd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_dggd(); if(j>0){message<<"    dggd   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_enWh(); if(j>0){message<<"    enWh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hnWe(); if(j>0){message<<"    hnWe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_enAW(); if(j>0){message<<"    enAW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AnWe(); if(j>0){message<<"    AnWe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_enZW(); if(j>0){message<<"    enZW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ZnWe(); if(j>0){message<<"    ZnWe   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_udWh(); if(j>0){message<<"    udWh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_uhWd(); if(j>0){message<<"    uhWd   : "<<j<<" failed tests."<<endl;n++;}
  //j = spinas::test_udAW(); if(j>0){message<<"    udAW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AWud(); if(j>0){message<<"    AWud   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gWud(); if(j>0){message<<"    gWud   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hhhh(); if(j>0){message<<"    hhhh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hhZZ(); if(j>0){message<<"    hhZZ   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hZZh(); if(j>0){message<<"    hZZh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hhWW(); if(j>0){message<<"    hhWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_hWWh(); if(j>0){message<<"    hWWh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AhWW(); if(j>0){message<<"    AhWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AWWh(); if(j>0){message<<"    AWWh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ZhWW(); if(j>0){message<<"    ZhWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ZWWh(); if(j>0){message<<"    ZWWh   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AAWW(); if(j>0){message<<"    AAWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AWAW(); if(j>0){message<<"    AWAW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AZWW(); if(j>0){message<<"    AZWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_AWZW(); if(j>0){message<<"    AWZW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ZZZZ(); if(j>0){message<<"    ZZZZ   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ZZWW(); if(j>0){message<<"    ZZWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_ZWZW(); if(j>0){message<<"    ZWZW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_WWWW(); if(j>0){message<<"    WWWW   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_WWWW2(); if(j>0){message<<"    WWWW2   : "<<j<<" failed tests."<<endl;n++;}
  j = spinas::test_gggg(); if(j>0){message<<"    gggg   : "<<j<<" failed tests."<<endl;n++;}
  
  
  cout<<"==============================================================================\n";
  cout<<"All Tests:     ";
  if(n>0){
    cout<<n<<" modules had failed tests:"<<endl;
    cout<<message.str();
    cout<<"                                                                        Fail!\n";
  }
  else    cout<<"                                                           Pass\n";
  cout<<"=============================================================================="<<endl;

  
}

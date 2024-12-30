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

//File:  SPINAS/source/starproduct.cpp

#include <iostream>
#include <sstream>
#include <cmath>

#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "sproduct.h"
#include "schain.h"
#include "supleproduct.h"
#include "starproduct.h"

namespace spinas {
  //Constructors
  starproduct::starproduct(){}

  //4 chains for spin-1
  starproduct::starproduct(schain *s0, schain *s1, particle *pp, schain *s2, schain *s3):
    N(2) {
    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
    p = pp;

    //Test that the chains all have the same Lorentz index exposed.
    check();
    //Update
    update();
  }

  //Check that the star product makes sense
  void starproduct::check(){
    bool isDotted = s[0]->is_Lindex_dotted();
    //Check that the chains on one side all have the same Lorentz index exposed.
    for(int i=1;i<N;i++)
      if(s[i]->is_Lindex_dotted()!=isDotted)
        throw std::runtime_error("Incorrect usage of starproduct.  All chains on one side must have the same Lorentz index exposed.");
    //Check that the chains on the other side have the opposite Lorentz index exposed.
    for(int i=N;i<2*N;i++)
      if(s[i]->is_Lindex_dotted()==isDotted)
        throw std::runtime_error("Incorrect usage of starproduct.  The chains on the other side must have the opposite Lorentz index exposed.");

    if(isDotted)
      isMomUpper = false;
    else
      isMomUpper = true;

    //Determine which chains are massive
    for(int i=0;i<2*N;i++)
      isMassive[i] = s[i]->is_right_massive();
  }

  //Update
  //Must be run after masses or momenta of particles is updated.
  void starproduct::update(){
    //Update isCalculated
    for(int i=0;i<N+1;i++)
      for(int j=0;j<N+1;j++)
        for(int k=0;k<N+1;k++)
          for(int l=0;l<N+1;l++)
            isCalculated[i][j][k][l] = false;
    if(isMomUpper)
      pMat = p->umat(N+1);
    else
      pMat = p->lmat(N+1);
  } 


  //products
  cdouble starproduct::v(){
    //Check that there are no spin indices.
    int nmassive=0;
    for(int i=0;i<N;i++)
      if(isMassive[i]) nmassive++;
    if(nmassive>0){
      throw std::runtime_error("Incorrect usage of starproduct.v().  All the particles at the ends must be massless.");
      return cdouble(0,0);
    }

    //Check whether it is already calculted.
    if(isCalculated[0][0][0][0]) return product[0][0][0][0];

    //Calcualte it.
    cvector vec[4];
    for(int i=0;i<2*N;i++)
      vec[i] = s[i]->v();
    product[0][0][0][0] = v(vec);
    isCalculated[0][0][0][0] = true;

    return product[0][0][0][0];
  }

   cdouble starproduct::v(const int& ds){
    //Check that there is 1 spin index.
    int nmassive=0, imassive=0;
    for(int i=0;i<2*N;i++)
      if(isMassive[i]){ 
        imassive=i;
        nmassive++;
      }
    if(nmassive>1||nmassive==0){
      throw std::runtime_error("Incorrect usage of starproduct.v(int j).  1 particle at the end of a chain must be massive.");
      return cdouble(0,0);
    }
    
    //Get index from spin
    int jmassive=(ds+2)/2;
    //Check whether it is already calculted.
    if(imassive==0&&isCalculated[jmassive][0][0][0]) return product[jmassive][0][0][0];
    if(imassive==1&&isCalculated[0][jmassive][0][0]) return product[0][jmassive][0][0];
    if(imassive==2&&isCalculated[0][0][jmassive][0]) return product[0][0][jmassive][0];
    if(imassive==3&&isCalculated[0][0][0][jmassive]) return product[0][0][0][jmassive];

    //Calcualte it.
    cvector vec[4];
    for(int i=0;i<2*N;i++){
      if(i==imassive)
        vec[i] = s[i]->v(ds);
      else
        vec[i] = s[i]->v();
    }
    if(imassive==0){
      product[jmassive][0][0][0] = v(vec);
      isCalculated[jmassive][0][0][0] = true;
      return product[jmassive][0][0][0];
    }
    if(imassive==1){
      product[0][jmassive][0][0] = v(vec);
      isCalculated[0][jmassive][0][0] = true;
      return product[0][jmassive][0][0];
    }
    if(imassive==2){
      product[0][0][jmassive][0] = v(vec);
      isCalculated[0][0][jmassive][0] = true;
      return product[0][0][jmassive][0];
    }
    if(imassive==3){
      product[0][0][0][jmassive] = v(vec);
      isCalculated[0][0][0][jmassive] = true;
      return product[0][0][0][jmassive];
    }
    return cdouble(0,0);
  }

  cdouble starproduct::v(const int& ds0, const int& ds1){
    //Check that there are 2 spin indices.
    int nmassive=0, imassive[2];
    int imi=0;
    for(int i=0;i<2*N;i++)
      if(isMassive[i]){ 
        imassive[imi++]=i;
        nmassive++;
      }
    if(nmassive>2||nmassive<2){
      throw std::runtime_error("Incorrect usage of starproduct.v(int j1, int j2).  2 particles at the end of chains must be massive.");
      return cdouble(0,0);
    }
    
    //Get index from spin
    int jmassive0=(ds0+2)/2;
    int jmassive1=(ds1+2)/2;
    //Check whether it is already calculted.
    if(imassive[0]==0&&imassive[1]==1&&isCalculated[jmassive0][jmassive1][0][0]) return product[jmassive0][jmassive1][0][0];
    if(imassive[0]==0&&imassive[1]==2&&isCalculated[jmassive0][0][jmassive1][0]) return product[jmassive0][0][jmassive1][0];
    if(imassive[0]==0&&imassive[1]==3&&isCalculated[jmassive0][0][0][jmassive1]) return product[jmassive0][0][0][jmassive1];
    if(imassive[0]==1&&imassive[1]==2&&isCalculated[0][jmassive0][jmassive1][0]) return product[0][jmassive0][jmassive1][0];
    if(imassive[0]==1&&imassive[1]==3&&isCalculated[0][jmassive0][0][jmassive1]) return product[0][jmassive0][0][jmassive1];
    if(imassive[0]==2&&imassive[1]==3&&isCalculated[0][0][jmassive0][jmassive1]) return product[0][0][jmassive0][jmassive1];


    //Calcualte it.
    cvector vec[4];
    for(int i=0;i<2*N;i++){
      if(i==imassive[0]) 
        vec[i] = s[i]->v(ds0);
      else if(i==imassive[1])
        vec[i] = s[i]->v(ds1);
      else
        vec[i] = s[i]->v();
    }
    if(imassive[0]==0&&imassive[1]==1){
      product[jmassive0][jmassive1][0][0] = v(vec);
      isCalculated[jmassive0][jmassive1][0][0] = true;
      return product[jmassive0][jmassive1][0][0];
    }
    if(imassive[0]==0&&imassive[1]==2){
      product[jmassive0][0][jmassive1][0] = v(vec);
      isCalculated[jmassive0][0][jmassive1][0] = true;
      return product[jmassive0][0][jmassive1][0];
    }
    if(imassive[0]==0&&imassive[1]==3){
      product[jmassive0][0][0][jmassive1] = v(vec);
      isCalculated[jmassive0][0][0][jmassive1] = true;
      return product[jmassive0][0][0][jmassive1];
    }
    if(imassive[0]==1&&imassive[1]==2){
      product[0][jmassive0][jmassive1][0] = v(vec);
      isCalculated[0][jmassive0][jmassive1][0] = true;
      return product[0][jmassive0][jmassive1][0];
    }
    if(imassive[0]==1&&imassive[1]==3){
      product[0][jmassive0][0][jmassive1] = v(vec);
      isCalculated[0][jmassive0][0][jmassive1] = true;
      return product[0][jmassive0][0][jmassive1];
    }
    if(imassive[0]==2&&imassive[1]==3){
      product[0][0][jmassive0][jmassive1] = v(vec);
      isCalculated[0][0][jmassive0][jmassive1] = true;
      return product[0][0][jmassive0][jmassive1];
    }
    return cdouble(0,0);
  }

  cdouble starproduct::v(const int& ds0, const int& ds1, const int& ds2){
    //Check that there are 2 spin indices.
    int nmassive=0, imassive[3];
    int imi=0;
    for(int i=0;i<2*N;i++)
      if(isMassive[i]){ 
        imassive[imi++]=i;
        nmassive++;
      }
    if(nmassive>3||nmassive<3){
      throw std::runtime_error("Incorrect usage of starproduct.v(int j1, int j2, int j3).  3 particles at the end of chains must be massive.");
      return cdouble(0,0);
    }
    
    //Get index from spin
    int jmassive0=(ds0+2)/2;
    int jmassive1=(ds1+2)/2;
    int jmassive2=(ds2+2)/2;
    //Check whether it is already calculted.
    if(imassive[0]==0&&imassive[1]==1&&imassive[2]==2&&isCalculated[jmassive0][jmassive1][jmassive2][0]) return product[jmassive0][jmassive1][jmassive2][0];
    if(imassive[0]==0&&imassive[1]==1&&imassive[2]==3&&isCalculated[jmassive0][jmassive1][0][jmassive2]) return product[jmassive0][jmassive1][0][jmassive2];
    if(imassive[0]==0&&imassive[1]==2&&imassive[2]==3&&isCalculated[jmassive0][0][jmassive1][jmassive2]) return product[jmassive0][0][jmassive1][jmassive2];
    if(imassive[0]==1&&imassive[1]==2&&imassive[2]==3&&isCalculated[0][jmassive0][jmassive1][jmassive2]) return product[0][jmassive0][jmassive1][jmassive2];

    //Calcualte it.
    cvector vec[4];
    for(int i=0;i<2*N;i++){
      if(i==imassive[0]) 
        vec[i] = s[i]->v(ds0);
      else if(i==imassive[1])
        vec[i] = s[i]->v(ds1);
      else if(i==imassive[2])
        vec[i] = s[i]->v(ds2);
      else
        vec[i] = s[i]->v();
    }
    if(imassive[0]==0&&imassive[1]==1&&imassive[2]==2){
      product[jmassive0][jmassive1][jmassive2][0] = v(vec);
      isCalculated[jmassive0][jmassive1][jmassive2][0] = true;
      return product[jmassive0][jmassive1][jmassive2][0];
    }
    if(imassive[0]==0&&imassive[1]==1&&imassive[2]==3){
      product[jmassive0][jmassive1][0][jmassive2] = v(vec);
      isCalculated[jmassive0][jmassive1][0][jmassive2] = true;
      return product[jmassive0][jmassive1][0][jmassive2];
    }
    if(imassive[0]==0&&imassive[1]==2&&imassive[2]==3){
      product[jmassive0][0][jmassive1][jmassive2] = v(vec);
      isCalculated[jmassive0][0][jmassive1][jmassive2] = true;
      return product[jmassive0][0][jmassive1][jmassive2];
    }
    if(imassive[0]==1&&imassive[1]==2&&imassive[2]==3){
      product[0][jmassive0][jmassive1][jmassive2] = v(vec);
      isCalculated[0][jmassive0][jmassive1][jmassive2] = true;
      return product[0][jmassive0][jmassive1][jmassive2];
    }
    return cdouble(0,0);
  }


  cdouble starproduct::v(const int& ds0, const int& ds1, const int& ds2, const int& ds3){
    //Check that there is 1 spin index.
    int nmassive=0;
    for(int i=0;i<2*N;i++)
      if(isMassive[i]){ 
        nmassive++;
      }
    if(nmassive>4||nmassive<4){
      throw std::runtime_error("Incorrect usage of starproduct.v(int j1, int j2, int j3, int j4).  4 particles at the end must be massive.");
      return cdouble(0,0);
    }
    
    //Get index from spin
    int jmassive0=(ds0+2)/2;
    int jmassive1=(ds1+2)/2;
    int jmassive2=(ds2+2)/2;
    int jmassive3=(ds3+2)/2;
    //Check whether it is already calculted.
    if(isCalculated[jmassive0][jmassive1][jmassive2][jmassive3]) return product[jmassive0][jmassive1][jmassive2][jmassive3];

    //Calcualte it.
    cvector vec[4];
    vec[0] = s[0]->v(ds0);
    vec[1] = s[1]->v(ds1);
    vec[2] = s[2]->v(ds2);
    vec[3] = s[3]->v(ds3);
    product[jmassive0][jmassive1][jmassive2][jmassive3] = v(vec);
    isCalculated[jmassive0][jmassive1][jmassive2][jmassive3] = true;
    return product[jmassive0][jmassive1][jmassive2][jmassive3];
  }


  //Calculate the triple product
  cdouble starproduct::v(cvector vec[4]){
    if(N!=2)
      throw std::runtime_error("Only spin-1 star products are implemented.");
    //The left chains are facing the wrong way.  (All schains are right facing.)  
    //Use the metric to switch their direction.
    //This means 0 and 2 are switched and 1 gets a minus sign.
    //Also, the matrix needs to be transposed.  See paper.
    return (vec[0].get(2)) * (-vec[1].get(1)) 
            * ( pMat.get(0,2) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,2) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,2) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,2) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,2) * vec[2].get(1) * vec[3].get(0) )
         + (-vec[0].get(1)) * (vec[1].get(0)) 
            * ( pMat.get(0,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,0) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,0) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,0) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,0) * vec[2].get(1) * vec[3].get(0) )
         + (vec[0].get(0)) * (vec[1].get(2)) 
            * ( pMat.get(0,1) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,1) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,1) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,1) * vec[2].get(1) * vec[3].get(0) )
         - (vec[0].get(0)) * (-vec[1].get(1)) 
            * ( pMat.get(0,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,0) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,0) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,0) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,0) * vec[2].get(1) * vec[3].get(0) )
         - (vec[0].get(2)) * (vec[1].get(0)) 
            * ( pMat.get(0,1) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,1) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,1) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,1) * vec[2].get(1) * vec[3].get(0) )
         - (-vec[0].get(1)) * (vec[1].get(2)) 
            * ( pMat.get(0,2) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,2) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,2) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,2) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,2) * vec[2].get(1) * vec[3].get(0) )
        ;
        /*vec[0].get(0) * vec[1].get(1) 
            * ( pMat.get(2,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(2,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(2,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(2,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,2) * vec[2].get(1) * vec[3].get(0) )
         + vec[0].get(1) * vec[1].get(2) 
            * ( pMat.get(0,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(0,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(0,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(0,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(0,2) * vec[2].get(1) * vec[3].get(0) )
         + vec[0].get(2) * vec[1].get(0) 
            * ( pMat.get(1,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(1,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(1,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(1,2) * vec[2].get(1) * vec[3].get(0) )
         - vec[0].get(2) * vec[1].get(1) 
            * ( pMat.get(0,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(0,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(0,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(0,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(0,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(0,2) * vec[2].get(1) * vec[3].get(0) )
         - vec[0].get(0) * vec[1].get(2) 
            * ( pMat.get(1,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(1,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(1,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(1,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(1,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(1,2) * vec[2].get(1) * vec[3].get(0) )
         - vec[0].get(1) * vec[1].get(0) 
            * ( pMat.get(2,0) * vec[2].get(1) * vec[3].get(2)
              + pMat.get(2,1) * vec[2].get(2) * vec[3].get(0)
              + pMat.get(2,2) * vec[2].get(0) * vec[3].get(1)
              - pMat.get(2,0) * vec[2].get(2) * vec[3].get(1)
              - pMat.get(2,1) * vec[2].get(0) * vec[3].get(2)
              - pMat.get(2,2) * vec[2].get(1) * vec[3].get(0) )
        ;*/
  }

}
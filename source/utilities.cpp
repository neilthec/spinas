
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

//File:  SPINAS/source/utilities.cpp

#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <random>
#include <fstream>
#include "types.h"
//#include "aliases.h"
#include "utilities.h"


namespace spinas {


  // double to Mathematica String
  std::string double_to_Mathematica_string(const ldouble& x) {
    std::stringstream tmpStr;
    std::string num;
    size_t epos;
    
    std::stringstream xStr;
    xStr << std::setprecision(std::numeric_limits<ldouble>::digits10 + 1) << x;
    num=xStr.str();
    epos=num.find("e");
    if(epos != std::string::npos)
      num.replace(epos,1," 10^");
    tmpStr<<num;
    
    return tmpStr.str();
  }
  
  

  //Rotate a Momentum by theta around axis given by normalized u[3]={ux,uy,uz}.
  void rotate_momentum(ldouble mom[4], const ldouble u[3], const ldouble& theta){
    ldouble cost = std::cos(theta), ux=u[0], uy=u[1], uz=u[2];
    ldouble sint = std::sqrt(1.0-cost*cost);
    ldouble R[4][4] = {{1,0,0,0},
		       {0, cost + ux*ux*(1 - cost), ux*uy*(1 - cost) - uz*sint, ux*uz*(1 - cost) + uy*sint},
		       {0, ux*uy*(1 - cost) + uz*sint, cost + uy*uy*(1 - cost), uy*uz*(1 - cost) - ux*sint},
		       {0, uz*ux*(1 - cost) - uy*sint, uz*uy*(1 - cost) + ux*sint,    cost + uz*uz*(1 - cost)}};
    ldouble nMom[4] = {0,0,0,0};
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++)
	nMom[i] += R[i][j]*mom[j];

    for(int i=0;i<4;i++) mom[i]=nMom[i];
  }
  //Boost a momentum by v[3]={vx,vy,vz};
  void boost_momentum(ldouble mom[4], const ldouble v[3]){
    ldouble vx=v[0], vy=v[1], vz=v[2];
    ldouble v2=vx*vx+vy*vy+vz*vz;
    ldouble gamma = 1./std::sqrt(1-v2);
    ldouble B[4][4] = {{gamma, -gamma*vx, -gamma*vy, -gamma*vz},
		       {-gamma*vx, 1 + (gamma - 1)*vx*vx/v2, (gamma - 1)*vx*vy/v2, (gamma - 1)*vx*vz/v2},
		       {-gamma*vy, (gamma - 1)*vx*vy/v2, 1 + (gamma - 1)*vy*vy/v2, (gamma - 1)*vy*vz/v2},
		       {-gamma*vz, (gamma - 1)*vx*vz/v2, (gamma - 1)*vz*vy/v2, 1 + (gamma - 1)*vz*vz/v2}};
    ldouble nMom[4] = {0,0,0,0};
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++)
	nMom[i] += B[i][j]*mom[j];

    for(int i=0;i<4;i++) mom[i]=nMom[i];
  }


  
  //Random Generation Used in Tests
  std::random_device rd;
  std::mt19937 gen(rd());  
  
  int choose_random_int(int begin, int end){
    std::uniform_int_distribution<> disInt(begin, end);
    return disInt(gen);
  }


  ldouble choose_random_ldouble(ldouble begin, ldouble end){
    std::uniform_real_distribution<ldouble> dis(begin, end);
    ldouble ld = dis(gen);
    return ld;
  }
  
  ldouble choose_random_momentum(ldouble p[4], ldouble begin, ldouble end){
    ldouble p2 = 0;
    for(int j=1;j<4;j++){
      p[j] = choose_random_ldouble(begin,end);
      p2 += p[j]*p[j];
    }
    // Ensure the range for p[0] is positive
    ldouble mass = choose_random_ldouble(0,end);
    p[0] = std::sqrt(mass*mass + p2);
    return mass;
  }
  
  void choose_random_massless_momentum(ldouble p[4], ldouble begin, ldouble end){
    ldouble p2 = 0;
    for(int j=1;j<4;j++){
      p[j] = choose_random_ldouble(begin,end);
      p2 += p[j]*p[j];
    }
    p[0] = std::sqrt(p2);
  }
  
  cdouble choose_random_cdouble(ldouble begin, ldouble end){
    ldouble ld1 = choose_random_ldouble(begin,end);
    ldouble ld2 = choose_random_ldouble(begin,end);
    return cdouble(ld1, ld2);
  }

  //Tools for File Importing

  //Break up string into multiple based on delimiter
  std::vector<std::string> split(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);

        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }

        return tokens;
    }

    std::string remove_whitespace(const std::string& original_string) {
        std::string new_string = original_string;
        new_string.erase(std::remove_if(new_string.begin(), new_string.end(), ::isspace), new_string.end());
        return new_string;
    }

    std::vector<std::string> remove_whitespace(const std::vector<std::string>& original_string_vector) {
        std::vector<std::string> new_string_vector;
        for (auto it = original_string_vector.begin(); it != original_string_vector.end(); ++it) 
            new_string_vector.push_back(remove_whitespace(*it));
        return new_string_vector;
    }

    //Removes spaces from start and end of string
    std::string trim(const std::string& str) {
        const std::string whitespace = " \t\n\r\f\v";

        size_t start = str.find_first_not_of(whitespace);
        //If no non-whitespace character is found, return an empty string
        if (start == std::string::npos) {
            return "";
        }

        size_t end = str.find_last_not_of(whitespace);

        return str.substr(start, end - start + 1);
    }

    std::vector<std::string> event_string_processor (const std::string& directory, const int& file_number) {
      // Construct filename (adjust naming convention as needed)
      std::string filename = directory + "/events" + std::to_string(file_number) + ".lhe";
      std::ifstream file(filename);
      
      std::vector<std::string> event_strings;
      if (!file.is_open()) {
          std::cerr << "Error: Could not open file " << filename << std::endl;
          return event_strings;
      }

      std::string line;
      std::string current_event;
      bool recording = false;

      while (std::getline(file, line)) {
          // Use your trim utility to handle potential leading spaces around tags
          std::string trimmed = trim(line);

          if (trimmed == "<event>") {
              recording = true;
              current_event.clear();
              continue; 
          }

          if (trimmed == "</event>") {
              recording = false;
              if (!current_event.empty()) {
                  event_strings.push_back(current_event);
              }
              continue;
          }

          if (recording) {
              // Append line and a newline to preserve the LHE structure for the stringstream
              current_event += line + "\n";
          }
      }

      file.close();
      return event_strings;
    }

    /* 
    Right, here's the plan. 

    1. Data structure at event.h & event.cpp is in place. 
        - make a constructor which takes string argument. 
    2. Here will be the logic for dividing up lhs file into multiple discreet events (from <event> to </event>)
    3. In tests, I will make a iterator which calls this, iterates through each event string, instnantiates an event, and tests it against spinas using the same values
        - Note that I will also need logic for getting N, IDs, and alpha values
        - Also note that I will want to iterate through all "events#.lhs" files 
    */
}

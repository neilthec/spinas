
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2025 Gabe Minney, Neil Christensen

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

#include "event.h"

namespace spinas {

    event::event(std::string event_string) {
        
        std::stringstream ss(event_string);

        // 1. Parse the first line (Global event info)
        if (!(ss >> n >> id >> w >> q >> a_em >> a_s)) {
            throw std::runtime_error("Malformed event string (header)");
        }

        // Pre-reserve memory for efficiency
        pdg.reserve(n);
        stat.reserve(n);
        mother.reserve(n);
        color.reserve(n);
        p.reserve(n);
        m.reserve(n);
        lt.reserve(n);
        spin.reserve(n);

        // 2. Parse n particles
        for (int i = 0; i < n; ++i) {
            int p_pdg, p_stat, m1, m2, c1, c2;
            ldouble px, py, pz, e, mass, lifetime, s;

            // The stream operator >> ignores all whitespace/newlines automatically
            if (ss >> p_pdg >> p_stat >> m1 >> m2 >> c1 >> c2 
                   >> px >> py >> pz >> e >> mass >> lifetime >> s) {
                
                pdg.push_back(p_pdg);
                stat.push_back(p_stat);
                mother.push_back({m1, m2});
                color.push_back({c1, c2});
                p.push_back({px, py, pz, e});
                m.push_back(mass);
                lt.push_back(lifetime);
                spin.push_back(s);
            }
        }
    }

    // Add these missing implementations:

    // Scalar getters
    int event::get_n() const {
        return n;
    }

    int event::get_id() const {
        return id;
    }

    ldouble event::get_w() const {
        return w;
    }

    ldouble event::get_q() const {
        return q;
    }

    ldouble event::get_a_em() const {
        return a_em;
    }

    ldouble event::get_a_s() const {
        return a_s;
    }

    // Vector getters with proper signatures matching header
    int event::get_pdg(const int& index) const {
        try {
            if (index < 0 || index >= static_cast<int>(pdg.size())) {
                throw std::out_of_range("Index out of range");
            }
            return pdg[index];
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 0;
        }
    }

    std::pair<int, int> event::get_mother(const int& index) const {
        try {
            if (index < 0 || index >= static_cast<int>(mother.size())) {
                throw std::out_of_range("Index out of range");
            }
            return mother[index];
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return std::make_pair(-1, -1);
        }
    }

    std::pair<int, int> event::get_color(const int& index) const {
        try {
            if (index < 0 || index >= static_cast<int>(color.size())) {
                throw std::out_of_range("Index out of range");
            }
            return color[index];
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return std::make_pair(-1, -1);
        }
    }

    std::array<ldouble, 4> event::get_p(const int& index) const {
        try {
            if (index < 0 || index >= static_cast<int>(p.size())) {
                throw std::out_of_range("Index out of range");
            }
            return p[index];
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return {0.0, 0.0, 0.0, 0.0};
        }
    }

    ldouble event::get_spin(const int& index) const {
        try {
            if (index < 0 || index >= static_cast<int>(spin.size())) {
                throw std::out_of_range("Index out of range");
            }
            return spin[index];
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 0.0;
        }
    }

    // Scalar setters (fix parameter types)
    void event::set_n(int& n_val) {
        n = n_val;
    }

    void event::set_id(int& id_val) {
        id = id_val;
    }

    void event::set_w(ldouble& w_val) {
        w = w_val;
    }

    void event::set_q(ldouble& q_val) {
        q = q_val;
    }

    void event::set_a_em(ldouble& a_em_val) {
        a_em = a_em_val;
    }

    void event::set_a_s(ldouble& a_s_val) {
        a_s = a_s_val;
    }

    void event::add_pdg(const int& index) {
        pdg.push_back(index);
    }

    void event::add_stat(const int& stat_val) {
        stat.push_back(stat_val);
    }

    // Vector adders
    void event::add_mother(const int& mother1, int& mother2) {
        mother.push_back(std::make_pair(mother1, mother2));
    }

    void event::add_color(const int& color1, int& color2) {
        color.push_back(std::make_pair(color1, color2));
    }

    void event::add_p(const ldouble& p0, const ldouble& p1, const ldouble& p2, const ldouble& p3) {
        p.push_back({p0, p1, p2, p3});
    }

    void event::add_m(const ldouble& m_val) {
        m.push_back(m_val);
    }

    void event::add_lt(const ldouble& lt_val) {
        lt.push_back(lt_val);
    }

    void event::add_spin(const ldouble& spin_val) {
        spin.push_back(spin_val);
    }

    // Size getters
    int event::get_pdg_size() const {
        return static_cast<int>(pdg.size());
    }

    int event::get_stat_size() const {
        return static_cast<int>(stat.size());
    }

    int event::get_mother_size() const {
        return static_cast<int>(mother.size());
    }

    int event::get_color_size() const {
        return static_cast<int>(color.size());
    }

    int event::get_p_size() const {
        return static_cast<int>(p.size());
    }

    int event::get_m_size() const {
        return static_cast<int>(m.size());
    }

    int event::get_lt_size() const {
        return static_cast<int>(lt.size());
    }

    /*
    
    Here's an example of how to implement an import function:

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <array>

typedef double ldouble; // Using double for standard precision

struct Event {
    int n, id;
    ldouble w, q, a_em, a_s;

    std::vector<int> pdg;
    std::vector<int> stat;
    std::vector<std::pair<int, int>> mother;
    std::vector<std::pair<int, int>> color;
    std::vector<std::array<ldouble, 4>> p;
    std::vector<ldouble> m;
    std::vector<ldouble> lt;
    std::vector<ldouble> spin;

    void clear() {
        pdg.clear(); stat.clear(); mother.clear();
        color.clear(); p.clear(); m.clear(); lt.clear(); spin.clear();
    }
};

int main() {
    std::ifstream file("events.lhe");
    std::string line;
    Event ev;

    while (std::getline(file, line)) {
        // Look for the start of an event
        if (line.find("<event>") != std::string::npos) {
            ev.clear();
            
            // 1. Read the first line of the event (Global data)
            if (!(file >> ev.n >> ev.id >> ev.w >> ev.q >> ev.a_em >> ev.a_s)) break;

            // 2. Loop through each particle
            for (int i = 0; i < ev.n; ++i) {
                int p_pdg, p_stat;
                int m1, m2, c1, c2;
                ldouble px, py, pz, e, mass, lifetime, s;

                file >> p_pdg >> p_stat >> m1 >> m2 >> c1 >> c2 
                     >> px >> py >> pz >> e >> mass >> lifetime >> s;

                ev.pdg.push_back(p_pdg);
                ev.stat.push_back(p_stat);
                ev.mother.push_back({m1, m2});
                ev.color.push_back({c1, c2});
                ev.p.push_back({px, py, pz, e});
                ev.m.push_back(mass);
                ev.lt.push_back(lifetime);
                ev.spin.push_back(s);
            }
            
            // Process the event (Example: print particle count)
            std::cout << "Processed event with ID " << ev.id << " and " << ev.n << " particles." << std::endl;
        }
    }
    return 0;
}*/
    
}
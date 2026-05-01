
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

//File:  SPINAS/source/events.h

#pragma once

#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <array>
#include "types.h" // Assuming ldouble is defined here

namespace spinas {
    // Format follows the Les Houches Event (LHE) version 3.0 convention
    class event {
    private:
        // --- Event Global Information ---
        int n;              // Number of particles in this event
        int id;             // Process ID
        ldouble w;          // Event weight
        ldouble q;          // Scale (Factorization scale in GeV)
        ldouble a_em;       // QED coupling alpha
        ldouble a_s;        // QCD coupling alpha_s

        // --- Particle Specific Information ---
        std::vector<int> pdg;               // PDG ID (e.g., 22 for photon, 11 for electron)
        std::vector<int> stat;              // Status (-1: incoming, 1: outgoing)
        std::vector<std::pair<int, int>> mother; // Indices of first and last mother particles
        std::vector<std::pair<int, int>> color;  // Color and anti-color flow tags
        std::vector<std::array<ldouble, 4>> p;   // Four-momentum (E, px, py, pz) in GeV
        std::vector<ldouble> m;             // Generated mass in GeV
        std::vector<ldouble> lt;            // Proper lifetime (c*tau) in mm
        std::vector<ldouble> spin;          // Helicity or spin projection (the 13th column)

    public:
        event(); //Default constructor
        event(std::string event_string); //Constructor from string

        int get_n() const;
        int get_id() const;
        ldouble get_w() const;
        ldouble get_q() const;
        ldouble get_a_em() const;
        ldouble get_a_s() const;

        int get_pdg(const int& index) const;
        int get_stat(const int& index) const;
        std::pair<int, int> get_mother(const int& index) const;
        std::pair<int, int> get_color(const int& index) const;
        std::array<ldouble, 4> get_p(const int& index) const;
        ldouble get_m(const int& index) const;
        ldouble get_lt(const int& index) const;
        ldouble get_spin(const int& index) const;

        void set_n(const int& n_val);
        void set_id(const int& id_val);
        void set_w(const ldouble& w_val);
        void set_q(const ldouble& q_val);
        void set_a_em(const ldouble& a_em_val);
        void set_a_s(const ldouble& a_s_val);

        void add_pdg(const int& pdg_val);
        void add_stat(const int& stat_val);
        void add_mother(const int& mother1, int& mother2);
        void add_color(const int& color1, int& color2);
        void add_p(const ldouble& p0, const ldouble& p1, const ldouble& p2, const ldouble& p3);
        void add_m(const ldouble& m_val);
        void add_lt(const ldouble& lt_val);
        void add_spin(const ldouble& spin_val);

        // Vector size getters
        int get_pdg_size() const;
        int get_stat_size() const;
        int get_mother_size() const;
        int get_color_size() const;
        int get_p_size() const;
        int get_m_size() const;
        int get_lt_size() const;
    };
}
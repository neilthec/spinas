//File:  SPINAS/boost/events.cpp

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <string>
#include "types.h"
#include "utilities.h"
//#include "spinas.h"
#include "cmatrix.h"
#include "cvector.h"
#include "event.h"
#include"particle.h"
#include "propagator.h"
#include "propagator.h"
#include "sproduct.h"


using namespace spinas;

void check_complex_relative(cdouble actual, cdouble expected, ldouble tol_percent) {
    // Threshold of 1e-11
    const ldouble tiny = static_cast<ldouble>(1e-11);

    if (std::abs(expected) < static_cast<ldouble>(1e-18)) {
        BOOST_CHECK_SMALL(actual.real(), tiny);
        BOOST_CHECK_SMALL(actual.imag(), tiny);
    } else {
        BOOST_CHECK_CLOSE(actual.real(), expected.real(), tol_percent);
        BOOST_CHECK_CLOSE(actual.imag(), expected.imag(), tol_percent);
    }
}

BOOST_AUTO_TEST_SUITE(event_tests)

/* Create a new BOOST_AUTO_TEST_CASE. Do basic tests for set and get with hardcoded numbers. Repeat with random momenta (see function somewhere).
Create a constructor that doesn't take a string and creates an empty event. Test a hardcoded event string using a constructor to read it in and give an event.
Add all the different parts of the event and test them.*/ 
BOOST_AUTO_TEST_CASE(test_event_manual_build) {

    event ev;

    int n = 2;
    ev.set_n(n);

    ev.add_pdg(11);
    ev.add_pdg(-11);

    ev.add_stat(-1);
    ev.add_stat(1);

    int zero = 0;
    ev.add_mother(zero, zero);
    ev.add_mother(zero, zero);

    ev.add_color(zero, zero);
    ev.add_color(zero, zero);

    // p = (px, py, pz, E)
    ev.add_p(0.0, 0.0, 10.0, 10.0);
    ev.add_p(0.0, 0.0, -10.0, 10.0);

    ev.add_m(0.0005);
    ev.add_m(0.0005);

    ev.add_lt(0.0);
    ev.add_lt(0.0);

    ev.add_spin(1.0);
    ev.add_spin(-1.0);

    BOOST_CHECK_EQUAL(ev.get_n(), 2);
    BOOST_CHECK_EQUAL(ev.get_pdg(0), 11);
    BOOST_CHECK_EQUAL(ev.get_pdg(1), -11);

    std::array<ldouble, 4> p0 = ev.get_p(0);
    BOOST_CHECK_CLOSE(p0[3], 10.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(test_event_random_momenta) {

    event ev;
    int n = 5;
    ev.set_n(n);

    for (int i = 0; i < n; ++i) {
        ldouble pvec[4];

        ldouble mass = choose_random_momentum(pvec, -10.0, 10.0);

        ev.add_p(pvec[0], pvec[1], pvec[2], pvec[3]);
        ev.add_m(mass);

        // Add minimal required structure
        ev.add_pdg(choose_random_int(-20, 20));
        ev.add_stat(1);

        int zero = 0;
        ev.add_mother(zero, zero);
        ev.add_color(zero, zero);
        ev.add_lt(0.0);
        ev.add_spin(0.0);
    }

    // Check on-shell condition: E^2 - |p|^2 = m^2
    for (int i = 0; i < n; ++i) {
        std::array<ldouble, 4> p = ev.get_p(i);
        ldouble m = ev.get_m(i);

        ldouble mass_sq =
            p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];

        BOOST_CHECK_CLOSE(mass_sq, m*m, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(test_event_default_constructor) {

    event ev;  // requires you implemented this

    BOOST_CHECK_EQUAL(ev.get_n(), 0);
    BOOST_CHECK_EQUAL(ev.get_pdg_size(), 0);
    BOOST_CHECK_EQUAL(ev.get_p_size(), 0);
    BOOST_CHECK_EQUAL(ev.get_m_size(), 0);
    BOOST_CHECK_EQUAL(ev.get_stat_size(), 0);
}

BOOST_AUTO_TEST_CASE(test_event_from_string) {

    std::string lhe_event =
        "2 1 1.0 1.0 0.007297 0.118\n"
        "11 -1 0 0 0 0 0.0 0.0 10.0 10.0 0.0005 0.0 1.0\n"
        "-11 1 0 0 0 0 0.0 0.0 -10.0 10.0 0.0005 0.0 -1.0\n";

    event ev(lhe_event);

    BOOST_CHECK_EQUAL(ev.get_n(), 2);
    BOOST_CHECK_EQUAL(ev.get_id(), 1);

    BOOST_CHECK_CLOSE(ev.get_w(), 1.0, 1e-9);
    BOOST_CHECK_CLOSE(ev.get_a_em(), 0.007297, 1e-6);

    BOOST_CHECK_EQUAL(ev.get_pdg(0), 11);
    BOOST_CHECK_EQUAL(ev.get_pdg(1), -11);

    std::array<ldouble, 4> p1 = ev.get_p(1);
    BOOST_CHECK_CLOSE(p1[2], -10.0, 1e-9);

    BOOST_CHECK_CLOSE(ev.get_m(0), 0.0005, 1e-9);
}

BOOST_AUTO_TEST_CASE(test_lhe_parsing_and_instantiation) {
    // 1. Setup pathing (Adjust directory as needed for your build system)
    std::string test_dir = "../event_files"; 
    int test_file_num = 1;

    // 2. Extract strings
    std::vector<std::string> event_strings = event_string_processor(test_dir, test_file_num);

    // Ensure we actually found events before proceeding
    BOOST_REQUIRE_MESSAGE(!event_strings.empty(), "No events found in file " << test_file_num);

    int processed_count = 0;

    for (const auto& raw_str : event_strings) {
        try {
            // 3. Instantiate Event
            // This populates the vectors and converts scientific notation to ldouble
            event ev(raw_str);

            // 4. Basic Validation
            // Verify that the number of particles parsed matches 'n' from the global line
            BOOST_CHECK_EQUAL(ev.get_pdg_size(), static_cast<size_t>(ev.get_n()));
            BOOST_CHECK_EQUAL(ev.get_p_size(), static_cast<size_t>(ev.get_n()));

            processed_count++;
        }
        catch (const std::exception& e) {
            BOOST_ERROR("Exception caught during Event instantiation: " << e.what());
        }
        
    }

    std::cout << "Successfully validated " << processed_count << " event objects." << std::endl;
}

BOOST_AUTO_TEST_CASE(test_spinas_matrix_elements) {
    const ldouble me = static_cast<ldouble>(0.0005);
    const ldouble echarge = static_cast<ldouble>(0.31333);
    const ldouble tol = static_cast<ldouble>(0.01); // 0.01 percent

    std::string test_dir = "../event_files"; 
    std::vector<std::string> event_strings = event_string_processor(test_dir, 1);
    BOOST_REQUIRE_MESSAGE(!event_strings.empty(), "Check path: " << test_dir);

    particle p1(me), p2(me);           
    particle p3(static_cast<ldouble>(0.0)), p4(static_cast<ldouble>(0.0)), p5(static_cast<ldouble>(0.0)); 

    sproduct s12s(SQUARE, &p1, &p2);
    sproduct a34a(ANGLE, &p3, &p4);
    sproduct s132a(SQUARE, &p1, &p3, &p2);

    propagator electron_prop(me, static_cast<ldouble>(0.0));

    for (const std::string& raw_str : event_strings) {
        event ev(raw_str); 

        ldouble moms[5][4];
        for (int i = 0; i < ev.get_n(); ++i) {
            auto m_vec = ev.get_p(i);
            int pdg_id = std::abs(ev.get_pdg(i));
            ldouble m = (pdg_id == 22) ? static_cast<ldouble>(0.0) : me;
            
            moms[i][1] = static_cast<ldouble>(m_vec[1]);
            moms[i][2] = static_cast<ldouble>(m_vec[2]);
            moms[i][3] = static_cast<ldouble>(m_vec[3]);
            moms[i][0] = std::sqrt(moms[i][1]*moms[i][1] + moms[i][2]*moms[i][2] + moms[i][3]*moms[i][3] + m*m);
            
            if (i == 0)      p1.set_momentum(moms[i]);
            else if (i == 1) p2.set_momentum(moms[i]);
            else if (i == 2) p3.set_momentum(moms[i]);
            else if (i == 3) p4.set_momentum(moms[i]);
            else if (i == 4) p5.set_momentum(moms[i]);
        }

        s12s.update();
        a34a.update();
        s132a.update();

        int s1 = -1, s2 = 1;
        cdouble val_electrons = s12s.v(s1, s2);
        cdouble val_photons   = a34a.v(); 

        ldouble prop_p[4];
        for(int j=0; j<4; j++) prop_p[j] = moms[0][j] - moms[2][j]; 
        cdouble pDen = electron_prop.denominator(prop_p);

        cdouble M = (static_cast<ldouble>(2.0) * echarge * echarge * me * val_electrons * val_photons) / pDen;

        ldouble calculated_w = std::norm(M);
        ldouble reference_w = static_cast<ldouble>(ev.get_w());

        if (std::abs(reference_w) < static_cast<ldouble>(1e-18)) {
            BOOST_CHECK_SMALL(calculated_w, static_cast<ldouble>(1e-11));
        } else {
            BOOST_CHECK_CLOSE(calculated_w, reference_w, tol);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
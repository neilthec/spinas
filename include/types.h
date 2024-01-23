
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

//File:  SPINAS/include/types.h

#pragma once



#if defined WITH_LONG_DOUBLE
typedef long double ldouble;
#else
typedef double ldouble;
#endif

#include <complex>
typedef std::complex<ldouble> cdouble;


constexpr bool ANGLE = true;
constexpr bool SQUARE = false;
constexpr bool UPPER = true;
constexpr bool LOWER = false;




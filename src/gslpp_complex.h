/***************************************************************************
 *   Copyright (C) 2007 by SUSYfit Collaboration                           *
 *   authors\susyfit.org                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/** @file   gsl_complex.h
    @author SUSYfit Collaboration
    @date   Thu Nov 21 13:16:42 2007
    @brief  Declarations of gslpp::complex
*/

#include <iostream>
#ifndef GSLPP_COMPLEX
#define GSLPP_COMPLEX
#ifndef __GSL_COMPLEX_H__
# include <gsl/gsl_complex.h>
#endif
#ifndef __GSL_COMPLEX_MATH_H__
# include <gsl/gsl_complex_math.h>
#endif

namespace gslpp
{
  /** @defgroup complex complex
   */
  /** @class gslpp::complex gslpp_complex.h <gslpp/gslpp_complex.h>
    *  Complex numbers.
    * @ingroup complex
    */
  class complex
  {
      gsl_complex _complex;
    public:
      static const complex& i();
      /** Constructor */
      complex();
      complex(const double& real, const double& imag, bool polar = false);
      /** Copy constructor */
      complex(const complex& z);
      complex(const double& a);
      /** Conversion constructor */
      complex(const gsl_complex* z);
      complex(const gsl_complex& z);
      /** Destructor */
      virtual ~complex();
      /** Check if this is purely real  */
      bool is_real() const;
      /** Check if this is purely imaginary  */
      bool is_imag() const;
      /** Assign */
      void assign(const double& real, const double& imag, bool polar);
      /** Set the real part */
      const double& real() const;
      /** Set imaginary part */
      const double& imag() const;
      /** Get the real part */
      double& real();
      /** Get the imaginary part */
      double& imag();
      /** Get arg */
      double arg() const;
      /** Get abs */
      double abs() const;
      /** Get the abs squared */
      double abs2() const;
      /** Get the log of the abs */
      double log_of_abs() const;
      /** Assign */
      complex& operator=(const complex& z);
      complex& operator=(const double& x);
      /** Compare */
      bool operator==(const complex& z1) const;
      bool operator!=(const complex& z1) const;
      /** Unary minus */
      complex operator-() const;
      /** Addition operator */
      complex operator+(const complex& z1);
      /** Subtraction operator */
      complex operator-(const complex& z1);
      /** Multiplication operator */
      complex operator*(const complex& z1);
      /** Division operator */
      complex operator/(const complex& z1);
      /** Addition assignment */
      complex& operator+=(const complex& z1);
      /** Subtraction assignment */
      complex& operator-=(const complex& z1);
      /** Multiplication assignment */
      complex& operator*=(const complex& z1);
      /** Division assignment */
      complex& operator/=(const complex& z1);
      /** Addition operator */
      complex operator+(const double& a);
      /** Subtraction assignment */
      complex operator-(const double& a);
      /** Multiplication operator */
      complex operator*(const double& a);
      /** Division operator */
      complex operator/(const double& a);
      /** Addition assignment  */
      complex& operator+=(const double& a);
      /** Subtraction assignment */
      complex& operator-=(const double& a);
      /** Multiplication assignment */
      complex& operator*=(const double& a);
      /** Division assignment */
      complex& operator/=(const double& a);
      /** Return complex conjugate  */
      complex conjugate() const;
      /** Return 1 / z */
      complex inverse() const;
      /** Conversion  */
      gsl_complex* as_gsl_type_ptr() const;
      gsl_complex& as_gsl_type();
      const gsl_complex& as_gsl_type() const;
      operator gsl_complex& ();
      operator const gsl_complex& () const;

      /** Friend functions */
      /**
        * @ingroup complex
        * @param output output stream
        * @param z1 Complex number
        * @return formatted output for complex
        */
      friend std::ostream& operator<<(std::ostream& output, const complex& z1);
      /** @{
       * @name Operations on complex numbers
       */
      /** Add a real and complex numbers
       * @ingroup complex
       * @param x1 Real number
       * @param z2 Complex number
       * @return @f$ x_1 + z_2 @f$
       */
      friend complex operator+(const double& x1, const complex& z2);

      /** Subtract a real and complex numbers
       * @ingroup complex
       * @param x1 Real number
       * @param z2 Complex number
       * @return @f$ x_1 - z_2 @f$
       */
      friend complex operator-(const double& x1, const complex& z2);

      /** Multiply a real and complex numbers
       * @ingroup complex
       * @param x1 Real number
       * @param z2 Complex number
       * @return @f$ x_1 z_2 @f$
       */
      friend complex operator*(const double& x1, const complex& z2);

      /** Divide a real and complex numbers
       * @ingroup complex
       * @param x1 Real number
       * @param z2 Complex number
       * @return @f$ x_1 / z_2 @f$
       */
      friend complex operator/(const double& x1,const complex& z2);
      /** @}
       */

      friend complex exp(const complex& z);
      friend complex log(const complex& z);
      friend complex log10(const complex& z);
      friend complex log(const complex& z,
                          const complex& b);
      friend complex sqrt(const complex& z);
      friend complex pow(const complex& z1,
                          const complex& z2);
      friend complex pow(const complex& z,
                          const double  x);
      friend complex sin(const complex& z);
      friend complex cos(const complex& z);
      friend complex tan(const complex& z);
      friend complex sec(const complex& z);
      friend complex csc(const complex& z);
      friend complex cot(const complex& z);
      friend complex arcsin(const complex& z);
      friend complex arccos(const complex& z);
      friend complex arctan(const complex& z);
      friend complex arcsec(const complex& z);
      friend complex arccsc(const complex& z);
      friend complex arccot(const complex& z);
      friend complex sinh(const complex& z);
      friend complex cosh(const complex& z);
      friend complex tanh(const complex& z);
      friend complex sech(const complex& z);
      friend complex csch(const complex& z);
      friend complex coth(const complex& z);
      friend complex arcsinh(const complex& z);
      friend complex arccosh(const complex& z);
      friend complex arctanh(const complex& z);
      friend complex arcsech(const complex& z);
      friend complex arccsch(const complex& z);
      friend complex arccoth(const complex& z);
  };
}
#endif
//
// EOF
//

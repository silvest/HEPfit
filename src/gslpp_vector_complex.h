/***************************************************************************
 *   Copyright (C) 2007 by SUSYfit Collaboration                           *
 *   authors@susyfit.org                                                   *
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
#ifndef GSLPP_VECTOR_COMPLEX_H
#define GSLPP_VECTOR_COMPLEX_H
#include <iostream>
#ifndef __GSL_VECTOR_COMPLEX_DOUBLE_H__
#include <gsl/gsl_vector_complex_double.h>
#endif
#ifndef GSLPP_COMPLEX_H
#include "gslpp_complex.h"
#endif
#ifndef GSLPP_VECTOR_BASE_H
#include "gslpp_vector_base.h"
#endif

namespace gslpp
{
  /**
   * @author SUSYfit Collaboration <authors@susyfit.org>
   */
  template <>
  class vector<complex>
  {
      gsl_vector_complex *_vector;

    public:
      /** Constructor  */
      vector(const size_t& size, const complex& z);
      vector(const size_t& size, const double& a);
      /** Copy constructor */
      vector(const vector<complex>& v);
      vector(const vector<double>& v);
      vector(const gsl_vector_complex& v);
      vector(const gsl_vector_complex* v);
      /** Destructor */
      ~vector();
      /** Get i-th element */
      const complex operator()(const size_t& i) const;
      /** Set i-th element */
//      complex& operator()(const size_t& i);
      /** Assign */
      vector<complex>& operator=(const vector<complex>& v);
      vector<complex>& operator=(double a);
      /** Assign element */
      void assign(const size_t& i, const complex& z);
      void assign(const size_t& i, const double& a);
      /** Get vector size */
      size_t size() const;
      /** Get Euclidean norm */
      double mod() const;
      /** Conversion */
      gsl_vector_complex* as_gsl_type_ptr() const;
      gsl_vector_complex& as_gsl_type();
      const gsl_vector_complex& as_gsl_type() const;
      /** Unary minus */
      vector<complex> operator-() const;
      /** Addition operator */
      vector<complex> operator+(const vector<complex>& v);
      /** Subtraction operator */
      vector<complex> operator-(const vector<complex>& v);
      /** Scalar product operator */
      complex operator*(const vector<complex>& v);
      /** Vector product operator */
//       vector<complex> operator^(const vector<complex>& v);
      /** Addition assignment */
      vector<complex>& operator+=(const vector<complex>& v);
      /** Subtraction assignment */
      vector<complex>& operator-=(const vector<complex>& v);
      /** Addition operator */
      vector<complex> operator+(const complex& z);
      /** Subtraction assignment */
      vector<complex> operator-(const complex& z);
      /** Multiplication operator */
      vector<complex> operator*(const complex& z);
      /** Division operator */
      vector<complex> operator/(const complex& z);
      /** Addition assignment  */
      vector<complex>& operator+=(const complex& z);
      /** Subtraction assignment */
      vector<complex>& operator-=(const complex& z);
      /** Multiplication assignment */
      vector<complex>& operator*=(const complex& z);
      /** Division assignment */
      vector<complex>& operator/=(const complex& z);
      /** Addition operator  */
      vector<complex> operator+(const double& a);
      /** Subtraction assignment */
      vector<complex> operator-(const double& a);
      /** Multiplication operator */
      vector<complex> operator*(const double& a);
      /** Division operator */
      vector<complex> operator/(const double& a);
      /** Addition assignment  */
      vector<complex>& operator+=(const double& a);
      /** Subtraction assignment */
      vector<complex>& operator-=(const double& a);
      /** Multiplication assignment */
      vector<complex>& operator*=(const double& a);
      /** Division assignment */
      vector<complex>& operator/=(const double& a);

      /** friend functions */
      friend std::ostream& operator<<(std::ostream& output, const vector<complex>& v);

      /** @{
       * @name Operations on vector<complex>
       */
      /** Add a complex number to a complex vector
       * @ingroup vector
       * @param z Complex number
       * @param v Complex vector
       * @return @f$ z + v @f$
       */
      friend vector<complex> operator+(const complex& z, vector<complex> v);

      /** Subtract a complex number from a complex vector
       * @ingroup vector
       * @param z Complex number
       * @param v Complex vector
       * @return @f$ z - v @f$
       */
      friend vector<complex> operator-(const complex& z, vector<complex> v);

      /** Multiply a complex number by complex vector
       * @ingroup vector
       * @param z Complex number
       * @param v Complex vector
       * @return @f$ z*v @f$
       */
      friend vector<complex> operator*(const complex& z, vector<complex> v);

      /** Add a real number to a complex vector
       * @ingroup vector
       * @param a Real number
       * @param v Complex vector
       * @return @f$ a + v @f$
       */
      friend vector<complex> operator+(const double& a, vector<complex> v);

      /** Subtract a complex vector from a real number
       * @ingroup vector
       * @param a Real number
       * @param v Complex vector
       * @return @f$ a - v @f$
       */
      friend vector<complex> operator-(const double& a, vector<complex> v);

      /** Multiply a real number by a complex vector
       * @ingroup vector
       * @param a Real number
       * @param v Complex vector
       * @return @f$ z*v @f$
       */
      friend vector<complex> operator*(const double& a, vector<complex> v);
      /** @}
       */
  };
}
#endif

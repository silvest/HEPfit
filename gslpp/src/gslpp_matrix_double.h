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
#ifndef GSLPP_MATRIX_DOUBLE_H
#define GSLPP_MATRIX_DOUBLE_H
#include <iostream>
#ifndef __GSL_MATRIX_DOUBLE_H__
#include <gsl/gsl_matrix_double.h>
#endif
#ifndef GSLPP_COMPLEX_H
#ifndef __GSL_MATRIX_DOUBLE_H__
#include <gsl/gsl_matrix_double.h>
#endif
#include "gslpp_complex.h"
#endif
#ifndef GSLPP_VECTOR_BASE_H
#include "gslpp_vector_base.h"
#endif
#ifndef GSLPP_VECTOR_DOUBLE_H
#include "gslpp_vector_double.h"
#endif
#ifndef GSLPP_VECTOR_COMPLEX_H
#include "gslpp_vector_complex.h"
#endif
#ifndef GSLPP_MATRIX_BASE_H
#include "gslpp_matrix_base.h"
#endif
#ifndef GSLPP_MATRIX_COMPLEX_H
#include "gslpp_matrix_complex.h"
#endif

namespace gslpp
{
  /**
   * @author SUSYfit Collaboration <authors@susyfit.org>
   */
  template <>
  class matrix<double>
    {
      gsl_matrix *_matrix;

    public:
      /** Constructor */
      matrix(const size_t& size_i, const size_t& size_j, const double& a);
      matrix(const size_t& size_i, const double& a);
      /** Copy constructor */
      matrix(const matrix<double>& m);
      matrix(const gsl_matrix& m);
      matrix(const gsl_matrix* m);
      matrix(const vector<double>& v);
      /** Destructor */
      ~matrix();
      /** Get element (i,j)*/
      const double operator()(const size_t& i, const size_t& j) const;
      /** Set i-th element */
      double& operator()(const size_t& i, const size_t& j);
      /** Assign */
      matrix<double>& operator=(const matrix<double>& m);
      matrix<double>& operator=(double a);
      void assign(const size_t& i, const size_t& j, const matrix<double>& a);
     /** Get matrix size */
      size_t size_i() const;
      size_t size_j() const;
      /** Identity matrix */
      static matrix<double> Id(size_t size);
      /** Transpose matrix */
      matrix<double> transpose();
      /** Inverse matrix */
      matrix<double> inverse();
      /** Eigenvalues and eigenvectors */
      void eigensystem(matrix<complex> &U, vector<complex> &S);
      /** Conversion */
      gsl_matrix* as_gsl_type_ptr() const;
      gsl_matrix& as_gsl_type();
      const gsl_matrix& as_gsl_type() const;
      /**
       * check whether two matrices are equal
       * @param m1 the first matrix
       * @param m2 the second matrix
       * @return true if equal, false otherwise
       */
      bool is_equal(const matrix<double>& m1, const matrix<double>& m2);
      /** Unary minus (matrix) */
      matrix<double> operator-() const;
      /** Addition operator (matrix) */
      matrix<double> operator+(const matrix<double>& m);
      /** Subtraction operator (matrix) */
      matrix<double> operator-(const matrix<double>& m);
      /** Product (matrix) */
      matrix<double> operator*(const matrix<double>& m);
      /** Multiplication (vector double) */
      vector<double> operator*(const vector<double>& v);
      /** Multiplication (vector complex) */
      vector<complex> operator*(const vector<complex>& v);
      /** Addition assignment (matrix) */
      matrix<double>& operator+=(const matrix<double>& m);
      /** Subtraction assignment (matrix) */
      matrix<double>& operator-=(const matrix<double>& m);
      /** Multiplication assignment (matrix) */
      matrix<double>& operator*=(const matrix<double>& m);
      /** Addition operator (double) */
      matrix<double> operator+(const double& a);
      /** Subtraction assignment (double) */
      matrix<double> operator-(const double& a);
      /** Multiplication operator (double) */
      matrix<double> operator*(const double& a);
      /** Division operator (double) */
      matrix<double> operator/(const double& a);
      /** Addition assignment (double) */
      matrix<double>& operator+=(const double& a);
      /** Subtraction assignment (double) */
      matrix<double>& operator-=(const double& a);
      /** Multiplication assignment (double) */
      matrix<double>& operator*=(const double& a);
      /** Division assignment (double) */
      matrix<double>& operator/=(const double& a);
      /** Addition operator (complex) */
      matrix<complex> operator+(const complex& z);
      /** Subtraction assignment (complex) */
      matrix<complex> operator-(const complex& z);
      /** Multiplication operator (complex) */
      matrix<complex> operator*(const complex& z);
      /** Division operator (complex) */
      matrix<complex> operator/(const complex& z);
      /** friend functions */
      friend std::ostream& operator<<(std::ostream& output, const matrix<double>& m);
      /** @{
       * @name Operations on matrix<double>
       */
      /** Add a real number to a real vector
       * @ingroup vector
       * @param a Real number
       * @param m Real vector
       * @return @f$ a + m @f$
       */
      friend matrix<double> operator+(const double& a, matrix<double> m);

      /** Subtract a real number from a real matrix
       * @ingroup matrix
       * @param a Real number
       * @param m Real matrix
       * @return @f$ a - m @f$
       */
      friend matrix<double> operator-(const double& a, matrix<double>  m);

      /** Multiply a real number by real matrix
       * @ingroup matrix
       * @param a Real number
       * @param m Real matrix
       * @return @f$ a*m @f$
       */
      friend matrix<double> operator*(const double& a, matrix<double> m);

      /** Multiply a real vector by a real matrix
       * @ingroup matrix
       * @param v Real vector
       * @param m Real matrix
       * @return @f$ v*m @f$
       */
      friend vector<double> operator*(const vector<double>& v, matrix<double> m);

      /** Multiply a complex vector by a real matrix
       * @ingroup matrix
       * @param v Complex vector
       * @param m Real matrix
       * @return @f$ v*m @f$
       */
      friend vector<complex> operator*(const vector<complex>& v, matrix<double> m);

      /** Add a complex number to a real matrix
       * @ingroup matrix
       * @param z Complex number
       * @param v Real matrix
       * @return @f$ z + v @f$
       */
      friend matrix<complex> operator+(const complex& z, matrix<double> v);

      /** Subtract a complex number from a real matrix
       * @ingroup matrix
       * @param z Complex number
       * @param v Real matrix
       * @return @f$ z - v @f$
       */
      friend matrix<complex> operator-(const complex& z, matrix<double> v);

      /** Multiply a complex number by a real matrix
       * @ingroup matrix
       * @param z Complex number
       * @param v Real matrix
       * @return @f$ a*v @f$
       */
      friend matrix<complex> operator*(const complex& z, matrix<double> v);
      /** @}
       */
  };

}

#endif

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
#ifndef GSLPP_MATRIX_COMPLEX_H
#define GSLPP_MATRIX_COMPLEX_H
#include <iostream>
#ifndef __GSL_MATRIX_COMPLEX_DOUBLE_H__
#include <gsl/gsl_matrix_complex_double.h>
#endif
#ifndef GSLPP_COMPLEX_H
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
#ifndef GSLPP_MATRIX_DOUBLE_H
#include "gslpp_matrix_double.h"
#endif

namespace gslpp
{
  /**
   * @author SUSYfit Collaboration <authors@susyfit.org>
   */
  template <>
	/**
	 * \class matrix<complex>
	 * \author M.C.
	 * \date 09/22/2010
	 * \file gslpp_matrix_complex.h
	 * \brief Complex matrices
	 */
      class matrix<complex>
  {
    gsl_matrix_complex *_matrix;

    public:
	/**
	 * \brief Constructor
	 * \param size_i number of rows
	 * \param size_j number of columns
	 * \param z initial value of entries (complex)
	 */
      matrix(const size_t& size_i, const size_t& size_j, const complex& z);
	/**
	 * \brief Constructor
	 * \param size_i number of rows
	 * \param size_j number of columns
	 * \param a initial value of entries (double)
	 */
      matrix(const size_t& size_i, const size_t& size_j, const double& a);
      /** Copy constructor */
      matrix(const matrix<complex>& m);
      matrix(const matrix<double>& m);
      matrix(const vector<complex>& v);
      matrix(const vector<double>& v);
      matrix(const gsl_matrix_complex& m);
      matrix(const gsl_matrix_complex* m);
      /** Destructor */
      ~matrix();
      /** Get i-th element */
      complex operator()(const size_t& i, const size_t& j) const;
      /** Set i-th element */
//      complex& operator()(const size_t& i, const size_t& j);
      /** Assign */
      matrix<complex>& operator=(const matrix<complex>& m);
      /** Assign element */
      void assign(const size_t& i, const size_t& j, const complex& z);
      void assign(const size_t& i, const size_t& j, const double& a);
      void assign(const size_t& i, const size_t& j, const matrix<complex>& z);
      void assign(const size_t& i, const size_t& j, const matrix<double>& a);
      /** Get matrix size */
      size_t size_i() const;
      size_t size_j() const;
      /** Identity matrix */
      static matrix<complex> Id(size_t size);
      /** Transpose matrix */
      matrix<complex> transpose();
      /** Hermitean conjugate matrix */
      matrix<complex> hconjugate();
      matrix<complex> inverse();
      /**
       * Eigenvalues and eigenvectors
       * @param U matrix<complex>& eigenvectors 
       * @param S vector<double>& eigenvalues
       */
      void eigensystem(matrix<complex> &U, vector<double> &S);
      /**
       * Singular Value Decomposition as U diagonalmatrix(S) V^+
       * @param U matrix<complex>&
       * @param V matrix<complex>&
       * @param S vector<double>&
       */
      void singularvalue(matrix<complex> &U, matrix<complex> &V, vector<double> &S);
      /** Conversion */
      gsl_matrix_complex* as_gsl_type_ptr() const;
      gsl_matrix_complex& as_gsl_type();
      const gsl_matrix_complex& as_gsl_type() const;
      /** Unary minus */
      matrix<complex> operator-() const;
      /** Addition operator */
      matrix<complex> operator+(const matrix<complex>& m);
      /** Subtraction operator */
      matrix<complex> operator-(const matrix<complex>& m);
      /** Multiplication operator */
      matrix<complex> operator*(const matrix<complex>& m);
      /** Addition operator */
      matrix<complex> operator+(const matrix<double>& m);
      /** Subtraction operator */
      matrix<complex> operator-(const matrix<double>& m);
      /** Multiplication operator */
      matrix<complex> operator*(const matrix<double>& m);
      /** Multiplication operator */
      vector<complex> operator*(const vector<complex>& v);
      vector<complex> operator*(const vector<double>& v);
      /** Addition assignment  */
      matrix<complex>& operator+=(const matrix<complex>& m);
      /** Subtraction assignment */
      matrix<complex>& operator-=(const matrix<complex>& m);
      /** Multiplication assignment */
      matrix<complex>& operator*=(const matrix<complex>& m);
      /** Addition operator  */
      matrix<complex> operator+(const complex& z);
      /** Subtraction assignment */
      matrix<complex> operator-(const complex& z);
      /** Multiplication operator */
      matrix<complex> operator*(const complex& z);
      /** Division operator */
      matrix<complex> operator/(const complex& z);
      /** Addition assignment  */
      matrix<complex>& operator+=(const complex& z);
      /** Subtraction assignment */
      matrix<complex>& operator-=(const complex& z);
      /** Multiplication assignment */
      matrix<complex>& operator*=(const complex& z);
      /** Division assignment */
      matrix<complex>& operator/=(const complex& z);
      /** Addition operator  */
      matrix<complex> operator+(const double& a);
      /** Subtraction assignment */
      matrix<complex> operator-(const double& a);
      /** Multiplication operator */
      matrix<complex> operator*(const double& a);
      /** Division operator */
      matrix<complex> operator/(const double& a);
      /** Addition assignment  */
      matrix<complex>& operator+=(const double& a);
      /** Subtraction assignment */
      matrix<complex>& operator-=(const double& a);
      /** Multiplication assignment */
      matrix<complex>& operator*=(const double& a);
      /** Division assignment */
      matrix<complex>& operator/=(const double& a);

      /** friend functions */
      friend std::ostream& operator<<(std::ostream& output, const matrix<complex>& v);

     /** @{
       * @name Operations on matrix<complex>
       */
      /** Add a double matrix to a complex matrix
       * @ingroup matrix
       * @param m1 Double matrix
       * @param m2 Complex matrix
       * @return @f$ m2 + m1 @f$
       */
      friend matrix<complex> operator+(matrix<double> m1, matrix<complex> m2);

     /** @{
       * @name Operations on matrix<complex>
       */
      /** Subtract a double matrix to a complex matrix
       * @ingroup matrix
       * @param m1 Double matrix
       * @param m2 Complex matrix
       * @return @f$ -m2 + m1 @f$
       */
      friend matrix<complex> operator-(matrix<double> m1, matrix<complex> m2);

      /** @{
       * @name Operations on matrix<complex>
       */
      /** Add a complex number to a complex matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Complex matrix
       * @return @f$ z + m @f$
       */
      friend matrix<complex> operator+(const complex& z, matrix<complex> m);

      /** Subtract a complex number from a complex matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Complex matrix
       * @return @f$ z - m @f$
       */
      friend matrix<complex> operator-(const complex& z, matrix<complex> m);

      /** Multiply a complex number by complex matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Complex matrix
       * @return @f$ z*m @f$
       */
      friend matrix<complex> operator*(const complex& z, matrix<complex> m);

      /** Multiply a complex vector by a complex matrix
       * @ingroup matrix
       * @param v Complex vector
       * @param m Complex matrix
       * @return @f$ v*m @f$
       */
      friend vector<complex> operator*(const vector<complex>& v, matrix<complex> m);

      /** Multiply a real vector by a complex matrix
       * @ingroup matrix
       * @param v Real vector
       * @param m Complex matrix
       * @return @f$ v*m @f$
       */
      friend vector<complex> operator*(const vector<double>& v, matrix<complex> m);

      /** Add a real number to a complex matrix
       * @ingroup matrix
       * @param a Real number
       * @param m Complex matrix
       * @return @f$ a + m @f$
       */
      friend matrix<complex> operator+(const double& a, matrix<complex> m);

      /** Subtract a complex matrix from a real number
       * @ingroup matrix
       * @param a Real number
       * @param m Complex matrix
       * @return @f$ a - m @f$
       */
      friend matrix<complex> operator-(const double& a, matrix<complex> m);

      /** Multiply a real number by a complex matrix
       * @ingroup matrix
       * @param a Real number
       * @param m Complex matrix
       * @return @f$ z*m @f$
       */
      friend matrix<complex> operator*(const double& a, matrix<complex> m);
      /** @}
       */
  };
}

#endif

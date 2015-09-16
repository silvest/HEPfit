/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

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
     * @class matrix<double>
     * @ingroup gslpp
     * @brief A class for constructing and defining operations on real matrices.
     * @author HEPfit Collaboration
     * @copyright GNU General Public License
     * @details This class defines some common operations on real matrices
     * using the <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
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
      double operator()(const size_t& i, const size_t& j) const;
      /** Set i-th element */
      double& operator()(const size_t& i, const size_t& j);
      /** Assign */
      matrix<double>& operator=(const matrix<double>& m);
      matrix<double>& operator=(double a);
      void assign(const size_t& i, const size_t& j, const double& a);
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
      /** Determinant matrix */
      double determinant();
      /** Eigenvalues and eigenvectors */
      void eigensystem(matrix<complex> &U, vector<complex> &S);
      /** Conversion */
      gsl_matrix* as_gsl_type_ptr() const;
      gsl_matrix& as_gsl_type();
      const gsl_matrix& as_gsl_type() const;
      /** check whether two matrices are equal */
      //bool is_equal(const matrix<double>& m1, const matrix<double>& m2);

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
       * @param m Real matrix
       * @return @f$ z + m @f$
       */
      friend matrix<complex> operator+(const complex& z, matrix<double> m);

      /** Subtract a complex number from a real matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Real matrix
       * @return @f$ z - m @f$
       */
      friend matrix<complex> operator-(const complex& z, matrix<double> m);

      /** Multiply a complex number by a real matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Real matrix
       * @return @f$ a*m @f$
       */
      friend matrix<complex> operator*(const complex& z, matrix<double> m);
      /** @}
       */
  };

}

#endif

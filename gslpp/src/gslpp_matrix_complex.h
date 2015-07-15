/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

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
     * @class matrix<complex>
     * @ingroup gslpp
     * @brief A class for constructing and defining operations on complex matrices.
     * @author HEPfit Collaboration
     * @copyright GNU General Public License
     * @details This class defines some common operations on complex matrices
     * using the <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
     */
  template <>
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
      matrix(const size_t& size_i, const complex& z);
	/**
	 * \brief Constructor
	 * \param size_i number of rows
	 * \param size_j number of columns
	 * \param a initial value of entries (double)
	 */
      matrix(const size_t& size_i, const size_t& size_j, const double& a);
      matrix(const size_t& size_i, const double& a);
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
      const complex operator()(const size_t& i, const size_t& j) const;
      /** Set i-th element */
//      complex& operator()(const size_t& i, const size_t& j);
      /** Assign */
      matrix<complex>& operator=(const matrix<complex>& m);
      /** Assign element */
      void assign(const size_t& i, const size_t& j, const complex& z);
      void assign(const size_t& i, const size_t& j, const double& a);
      void assign(const size_t& i, const size_t& j, const matrix<complex>& z);
      void assign(const size_t& i, const size_t& j, const matrix<double>& a);
      void assignre(const size_t& i, const size_t& j,const double& a);
      void assignim(const size_t& i, const size_t& j,const double& a);
      /** Get matrix size */
      size_t size_i() const;
      size_t size_j() const;
      /** Identity matrix */
      static matrix<complex> Id(size_t size);
      /** Transpose matrix */
      matrix<complex> transpose() const;
      /** Hermitean conjugate matrix */
      matrix<complex> hconjugate() const;
      matrix<complex> inverse() const;
      /** Get matrix of real parts */
      matrix<double> real() const;
      /** Get matrix of imaginary parts */
      matrix<double> imag() const;
      /**
       * Eigenvalues and eigenvectors
       * @param U matrix<complex>& eigenvectors 
       * @param S vector<double>& eigenvalues
       */
      void eigensystem(matrix<complex> &U, vector<double> &S) const;
      /**
       * Singular Value Decomposition as U diagonalmatrix(S) V^+
       * @param U matrix<complex>&
       * @param V matrix<complex>&
       * @param S vector<double>&
       */
      void singularvalue(matrix<complex> &U, matrix<complex> &V, vector<double> &S) const;
      /** Conversion */
      gsl_matrix_complex* as_gsl_type_ptr() const;
      gsl_matrix_complex& as_gsl_type();
      const gsl_matrix_complex& as_gsl_type() const;
      /** check whether two matrices are equal */
      //bool is_equal(const matrix<complex>& m1, const matrix<complex>& m2);
      /** Unary minus */
      matrix<complex> operator-() const;
      /** Addition operator */
      matrix<complex> operator+(const matrix<complex>& m) const;
      /** Subtraction operator */
      matrix<complex> operator-(const matrix<complex>& m) const;
      /** Multiplication operator */
      matrix<complex> operator*(const matrix<complex>& m) const;
      /** Addition operator */
      matrix<complex> operator+(const matrix<double>& m) const;
      /** Subtraction operator */
      matrix<complex> operator-(const matrix<double>& m) const;
      /** Multiplication operator */
      matrix<complex> operator*(const matrix<double>& m) const;
      /** Multiplication operator */
      vector<complex> operator*(const vector<complex>& v) const;
      vector<complex> operator*(const vector<double>& v) const;
      /** Addition assignment  */
      matrix<complex>& operator+=(const matrix<complex>& m);
      /** Subtraction assignment */
      matrix<complex>& operator-=(const matrix<complex>& m);
      /** Multiplication assignment */
      matrix<complex>& operator*=(const matrix<complex>& m);
      /** Addition operator  */
      matrix<complex> operator+(const complex& z) const;
      /** Subtraction operator */
      matrix<complex> operator-(const complex& z) const;
      /** Multiplication operator */
      matrix<complex> operator*(const complex& z) const;
      /** Division operator */
      matrix<complex> operator/(const complex& z) const;
      /** Addition assignment  */
      matrix<complex>& operator+=(const complex& z);
      /** Subtraction assignment */
      matrix<complex>& operator-=(const complex& z);
      /** Multiplication assignment */
      matrix<complex>& operator*=(const complex& z);
      /** Division assignment */
      matrix<complex>& operator/=(const complex& z);
      /** Addition operator  */
      matrix<complex> operator+(const double& a) const;
      /** Subtraction operator */
      matrix<complex> operator-(const double& a) const;
      /** Multiplication operator */
      matrix<complex> operator*(const double& a) const;
      /** Division operator */
      matrix<complex> operator/(const double& a) const;
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
      friend matrix<complex> operator+(matrix<double> m1, const matrix<complex> m2);

     /** @{
       * @name Operations on matrix<complex>
       */
      /** Subtract a double matrix to a complex matrix
       * @ingroup matrix
       * @param m1 Double matrix
       * @param m2 Complex matrix
       * @return @f$ -m2 + m1 @f$
       */
      friend matrix<complex> operator-(matrix<double> m1, const matrix<complex> m2);

      /** @{
       * @name Operations on matrix<complex>
       */
      /** Add a complex number to a complex matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Complex matrix
       * @return @f$ z + m @f$
       */
      friend matrix<complex> operator+(const complex& z, const matrix<complex> m);

      /** Subtract a complex number from a complex matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Complex matrix
       * @return @f$ z - m @f$
       */
      friend matrix<complex> operator-(const complex& z, const matrix<complex> m);

      /** Multiply a complex number by complex matrix
       * @ingroup matrix
       * @param z Complex number
       * @param m Complex matrix
       * @return @f$ z*m @f$
       */
      friend matrix<complex> operator*(const complex& z, const matrix<complex> m);

      /** Multiply a complex vector by a complex matrix
       * @ingroup matrix
       * @param v Complex vector
       * @param m Complex matrix
       * @return @f$ v*m @f$
       */
      friend vector<complex> operator*(const vector<complex>& v, const matrix<complex> m);

      /** Multiply a real vector by a complex matrix
       * @ingroup matrix
       * @param v Real vector
       * @param m Complex matrix
       * @return @f$ v*m @f$
       */
      friend vector<complex> operator*(const vector<double>& v, const matrix<complex> m);

      /** Add a real number to a complex matrix
       * @ingroup matrix
       * @param a Real number
       * @param m Complex matrix
       * @return @f$ a + m @f$
       */
      friend matrix<complex> operator+(const double& a, const matrix<complex> m);

      /** Subtract a complex matrix from a real number
       * @ingroup matrix
       * @param a Real number
       * @param m Complex matrix
       * @return @f$ a - m @f$
       */
      friend matrix<complex> operator-(const double& a, const matrix<complex> m);

      /** Multiply a real number by a complex matrix
       * @ingroup matrix
       * @param a Real number
       * @param m Complex matrix
       * @return @f$ z*m @f$
       */
      friend matrix<complex> operator*(const double& a, const matrix<complex> m);
      /** @}
       */
  };
}

#endif

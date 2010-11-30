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
#ifndef __GSL_BLAS_H__
#include <gsl/gsl_blas.h>
#endif
#ifndef __GSL_LINALG_H__
#include <gsl/gsl_linalg.h>
#endif
#ifndef GSLPP_MATRIX_DOUBLE_H
#include "gslpp_matrix_double.h"
#endif

namespace gslpp
{
  /** Constructor  */
  matrix<double>::matrix(const size_t& size_i, const size_t& size_j, const double &a)
  {
    _matrix = gsl_matrix_alloc(size_i, size_j);
    gsl_matrix_set_all(_matrix, a);
  }

  /** Copy constructor  */
  matrix<double>::matrix(const matrix<double>& m)
  {
    _matrix = gsl_matrix_alloc(m.size_i(), m.size_j());
    gsl_matrix_memcpy(_matrix, m.as_gsl_type_ptr());
  }

  matrix<double>::matrix(const gsl_matrix& m)
  {
    _matrix = gsl_matrix_alloc(m.size1,m.size2);
    gsl_matrix_memcpy(_matrix, &m);
  }

  matrix<double>::matrix(const gsl_matrix* m)
  {
    _matrix = gsl_matrix_alloc(m->size1,m->size2);
    gsl_matrix_memcpy(_matrix, m);
  }

  /** Destructor */
  matrix<double>::~matrix()
  {
    gsl_matrix_free(_matrix);
  }

  /** Get element (i,j) */
  const double matrix<double>::operator()(const size_t& i, const size_t& j) const
  {
    const double *x = gsl_matrix_const_ptr(_matrix, i, j);
    return *x;
  }

  /** Set element (i,j) */
  double& matrix<double>::operator()(const size_t& i, const size_t& j)
  {
    double *x = gsl_matrix_ptr(_matrix, i, j);
    return *x;
  }

  /** Assign */
  matrix<double>& matrix<double>::operator=(const matrix<double>& m)
  {
    gsl_matrix_memcpy(_matrix, m.as_gsl_type_ptr());
    return *this;
  }

  /** Get matrx size  */
  size_t matrix<double>::size_i() const
  {
   return _matrix->size1;
  }

  size_t matrix<double>::size_j() const
  {
    return _matrix->size2;
  }

  /** Identity matrix */
  matrix<double> matrix<double>::Id(size_t size)
  {
    matrix<double> m1(size, size, 0.);
    if (gsl_matrix_add_diagonal(m1.as_gsl_type_ptr(), 1.))
    {
      std::cout << "\n Error in matrix<double> Id" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Transpose matrix */
  matrix<double> matrix<double>::transpose()
  {
    matrix<double> m1(size_j(), size_i(), 0.);
    if (gsl_matrix_transpose_memcpy(m1.as_gsl_type_ptr(), _matrix))
    {
      std::cout << "\n Error in matrix<double> transpose" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Inverse matrix */
  matrix<double> matrix<double>::inverse()
  {
    matrix<double> m1(_matrix);
    matrix<double> m2(_matrix);
    int signum;
    gsl_permutation *p;

    if (size_j() != size_i())
    {
      std::cout << "\n ** Size mismatch in matrix<double> inverse" << std::endl;
      exit(EXIT_FAILURE);
    }

    if ((p = gsl_permutation_alloc(size_i())) == NULL)
    {
      std::cout << "\n ** Error in matrix<double> inverse" << std::endl;
      exit(EXIT_FAILURE);
    }

    if(gsl_linalg_LU_decomp (m1.as_gsl_type_ptr(), p, &signum))
    {
      std::cout << "\n ** Error in matrix<double> inverse" << std::endl;
      gsl_permutation_free(p);
      exit(EXIT_FAILURE);
    }

    if(gsl_linalg_LU_invert(m1.as_gsl_type_ptr(), p, m2.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<double> inverse" << std::endl;
      gsl_permutation_free(p);
      exit(EXIT_FAILURE);
    }
    gsl_permutation_free(p);
    return m2;
  }

  /** Conversion */
  gsl_matrix* matrix<double>::as_gsl_type_ptr() const
  {
    return _matrix;
  }

  gsl_matrix& matrix<double>::as_gsl_type()
  {
    return const_cast<gsl_matrix&>(*_matrix);
  }

  const gsl_matrix& matrix<double>::as_gsl_type() const
  {
    return const_cast<gsl_matrix&>(*_matrix);
  }

  /** Unary minus */
  matrix<double> matrix<double>::operator-() const
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_scale(m1.as_gsl_type_ptr(), -1.))
    {
      std::cout << "\n Error in matrix<double> unary -" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Addition operator (matrix) */
  matrix<double> matrix<double>::operator+(const matrix<double>& m)
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_add(m1.as_gsl_type_ptr(), m.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in matrix<double> +" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Subtraction operator (matrix) */
  matrix<double> matrix<double>::operator-(const matrix<double>& m)
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_sub(m1.as_gsl_type_ptr(), m.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in matrix<double> -" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Multiplication operator (matrix) */
  matrix<double> matrix<double>::operator*(const matrix<double>& m)
  {
    int i,j,k;
    matrix<double> m1(size_i(),m.size_i(),0.);

    if (size_j() != m.size_i())
    {
      std::cout << "\n Error in matrix<double> *" << std::endl;
      exit(EXIT_FAILURE);
    }

    for(i=0; i<size_i(); i++)
      for(j=0; j<m.size_j(); j++)
        for(k=0; k<m.size_i(); k++)
          m1(i,j)+=(*this)(i,k)*m(k,j);

    return m1;
  }

  /** Multiplication assignment (vector double) */
  vector<double> matrix<double>::operator*(const vector<double>& v)
  {
    vector<double> v1(size_i(),0.);

    if (size_j() != v.size())
    {
      std::cout << "\n ** Size mismatch in matrix<double> * (vector double)"
          << std::endl;
      exit(EXIT_FAILURE);
    }
    if(gsl_blas_dgemv(CblasNoTrans, 1., _matrix, v.as_gsl_type_ptr(),
       0., v1.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<double> * (vector double)"
          << std::endl;
      exit(EXIT_FAILURE);
    }
    return v1;
  }

  /** Multiplication assignment (vector complex) */
  vector<complex> matrix<double>::operator*(const vector<complex>& v)
  {
    matrix<complex> m1(*this);
    return m1 * v;
  }

  /** Addition assignment (matrix) */
  matrix<double>& matrix<double>::operator+=(const matrix<double>& m)
  {
    *this = *this + m;
    return *this;
  }
  /** Subtraction assignment (matrix) */
  matrix<double>& matrix<double>::operator-=(const matrix<double>& m)
  {
    *this = *this - m;
    return *this;
  }

  /** Multiplication assignment (matrix) */
  matrix<double>& matrix<double>::operator*=(const matrix<double>& m)
  {
    if (!((size_i() == size_j()) && (size_i() == m.size_i()) && (size_j() == m.size_j())))
    {
      std::cout << "\n Error in matrix<double> *= (matrix)" << std::endl;
      exit(EXIT_FAILURE);
    }
    *this = *this * m;
    return *this;
  }

  /** Addition operator (double) */
  matrix<double> matrix<double>::operator+(const double& a)
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_add_constant(m1.as_gsl_type_ptr(), a))
    {
      std::cout << "\n Error in matrix<double> + (double)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Subtraction operator (double) */
  matrix<double> matrix<double>::operator-(const double& a)
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_add_constant(m1.as_gsl_type_ptr(), -a))
    {
      std::cout << "\n Error in matrix<double> - (double)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Multiplication operator (double) */
  matrix<double> matrix<double>::operator*(const double& a)
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_scale(m1.as_gsl_type_ptr(), a))
    {
      std::cout << "\n Error in matrix<double> * (double)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Division operator (double) */
  matrix<double> matrix<double>::operator/(const double& a)
  {
    matrix<double> m1(_matrix);
    if (gsl_matrix_scale(m1.as_gsl_type_ptr(), 1./a))
    {
      std::cout << "\n Error in matrix<double> / (double)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Addition assignment (double) */
  matrix<double>& matrix<double>::operator+=(const double& a)
  {
    *this = *this + a;
    return *this;
  }
  /** Subtraction assignment (double) */
  matrix<double>& matrix<double>::operator-=(const double& a)
  {
    *this = *this - a;
    return *this;
  }
  /** Multiplication assignment (double) */
  matrix<double>& matrix<double>::operator*=(const double& a)
  {
    *this = *this * a;
    return *this;
  }
  /** Division assignment (double) */
  matrix<double>& matrix<double>::operator/=(const double& a)
  {
    *this = *this / a;
    return *this;
  }
  /** Addition operator (complex) */
  matrix<complex> matrix<double>::operator+(const complex& z)
  {
    matrix<complex> v1(*this);
    return v1+z;
  }
  /** Subtraction operator (complex) */
  matrix<complex> matrix<double>::operator-(const complex& z)
  {
    matrix<complex> v1(*this);
    return v1-z;
  }
  /** Multiplication operator (complex) */
  matrix<complex> matrix<double>::operator*(const complex& z)
  {
    matrix<complex> v1(*this);
    return v1*z;
  }
  /** Division operator (complex) */
  matrix<complex> matrix<double>::operator/(const complex& z)
  {
    matrix<complex> v1(*this);
    return v1/z;
  }

  /** friend functions */
  std::ostream& operator<<(std::ostream& output, const matrix<double>& m)
  {
    size_t i,j;
    for (i=0; i<m.size_i(); i++)
      {
        output << std::endl;
        output << "\t(";
        for (j=0; j<m.size_j()-1; j++)
          output << m(i,j) << ",";
        output << m(i,j) << ")";
      }
    return output;
  }
  /** @{
   * @name Operations on matrix<double>
   */
  /** Add a real number to a real matrix
   * @ingroup matrix
   * @param a Real number
   * @param m Real matrix
   * @return @f$ a + m @f$
   */
  matrix<double> operator+(const double& a, matrix<double> m)
  {
    return m+a;
  }

  /** Subtract a real number from a real matrix
   * @ingroup matrix
   * @param a Real number
   * @param m Real matrix
   * @return @f$ a - m @f$
   */
  matrix<double> operator-(const double& a, matrix<double> m)
  {
    return -m+a;
  }

  /** Multiply a real number by a real matrix
   * @ingroup matrix
   * @param a Real number
   * @param m Real matrix
   * @return @f$ a*m @f$
   */
  matrix<double> operator*(const double& a, matrix<double> m)
  {
    return m*a;
  }

  /** Multiply a real vector by a real matrix
   * @ingroup matrix
   * @param v Real vector
   * @param m Real matrix
   * @return @f$ v*m @f$
   */
  vector<double> operator*(const vector<double>& v, matrix<double> m)
  {
    return m.transpose() * v;
  }

  /** Multiply a complex vector by a real matrix
   * @ingroup matrix
   * @param v Complex vector
   * @param m Real matrix
   * @return @f$ v*m @f$
   */
  vector<complex> operator*(const vector<complex>& v, matrix<double> m)
  {
    return m.transpose() * v;
  }

  /** Add a complex number to a real matrix
   * @ingroup matrix
   * @param z Complex number
   * @param m Real matrix
   * @return @f$ z + v @f$
   */
  matrix<complex> operator+(const complex& z, matrix<double> m)
  {
    return m+z;
  }

  /** Subtract a complex number from a real matrix
   * @ingroup matrix
   * @param z Complex number
   * @param m Real matrix
   * @return @f$ z - m @f$
   */
  matrix<complex> operator-(const complex& z, matrix<double> m)
  {
    return -m+z;
  }

  /** Multiply a complex number by a real matrix
   * @ingroup matrix
   * @param z Complex number
   * @param m Real matrix
   * @return @f$ z*m @f$
   */
  matrix<complex> operator*(const complex& z, matrix<double> m)
  {
    return m*z;
  }
  /** @}
   */
}

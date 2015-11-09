/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef __GSL_BLAS_H__
#include <gsl/gsl_blas.h>
#endif
#ifndef __GSL_LINALG_H__
#include <gsl/gsl_linalg.h>
#endif
#ifndef __GSL_EIGEN_H__
#include <gsl/gsl_eigen.h>
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

  matrix<double>::matrix(const size_t& size_i, const double &a)
  {
    _matrix = gsl_matrix_alloc(size_i, size_i);
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

  matrix<double>::matrix(const vector<double>& v)
  {
    size_t i, size = v.size();
    _matrix = gsl_matrix_alloc(size, size);
    gsl_matrix_set_all(_matrix, 0.);
    for(i=0; i<size; ++i)
        gsl_matrix_set(_matrix, i, i, v(i));
  }

  /** Destructor */
  matrix<double>::~matrix()
  {
    gsl_matrix_free(_matrix);
  }

  /** Get element (i,j) */
  double matrix<double>::operator()(const size_t& i, const size_t& j) const
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
      if(size_i()==m.size_i() && size_j()==m.size_j())
           gsl_matrix_memcpy(_matrix, m.as_gsl_type_ptr());
      else
      {
          std::cout << "\n ** Wrong size assign in matrix<complex> operator =" << std::endl;
          exit(EXIT_FAILURE);
      }
    return *this;
  }
  
  matrix<double>& matrix<double>::operator=(double a)
  {
    gsl_matrix_set_all(_matrix, a);
    return *this;
  }

  void matrix<double>::assign(const size_t& i, const size_t& j, const double& a)
  {
    gsl_matrix_set(_matrix, i, j, a);
  }

  /** Assign submatrix */
  void matrix<double>::assign(const size_t& i, const size_t& j, const matrix<double>& a)
  {
      size_t ki,kj;
      double *x;
      if(i+a.size_i() <= size_i() && j+a.size_j() <= size_j())
          for(ki=i;ki<i+a.size_i();ki++)
              for(kj=j;kj<j+a.size_j();kj++)
              {
                  x = gsl_matrix_ptr(_matrix, ki, kj);
                  *x = a(ki-i,kj-j);
              }
      else
      {
          std::cout << "\n ** Wrong size assign in matrix<double> assign submatrix" << std::endl;
          exit(EXIT_FAILURE);
      }
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
  
  /** Determinant of matrix */
  double matrix<double>::determinant(){
      
    matrix<double> m1(_matrix);  
    int signum;
    gsl_permutation *p;
    
    if (size_j() != size_i())
    {
      std::cout << "\n ** Size mismatch in matrix<double> determinant" << std::endl;
      exit(EXIT_FAILURE);
    }
   
    if ((p = gsl_permutation_alloc(size_i())) == NULL)
    {
      std::cout << "\n ** Error in matrix<double> determinant" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if(gsl_linalg_LU_decomp (m1.as_gsl_type_ptr(), p, &signum))
    {
      std::cout << "\n ** Error in matrix<double> determinant" << std::endl;
      gsl_permutation_free(p);
      exit(EXIT_FAILURE);
    }
    gsl_permutation_free(p);
    return gsl_linalg_LU_det(m1.as_gsl_type_ptr() , signum);
  }
  
  void matrix<double>::eigensystem(matrix<complex> &U, vector<complex> &S) {
      matrix<double> m(*this);

      gsl_eigen_nonsymmv_workspace *ws;

      ws = gsl_eigen_nonsymmv_alloc(size_i());

      gsl_eigen_nonsymmv(m.as_gsl_type_ptr(), S.as_gsl_type_ptr(),
              U.as_gsl_type_ptr(), ws);

      gsl_eigen_nonsymmv_free(ws);
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

//  bool matrix<double>::is_equal(const matrix<double>& m1, const matrix<double>& m2){
//      return gsl_matrix_equal (m1.as_gsl_type_ptr(),m2.as_gsl_type_ptr());
//  }
  
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
    unsigned int i,j,k;
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

  matrix<double> operator+(const double& a, matrix<double> m)
  {
    return m+a;
  }

  matrix<double> operator-(const double& a, matrix<double> m)
  {
    return -m+a;
  }

  matrix<double> operator*(const double& a, matrix<double> m)
  {
    return m*a;
  }

  vector<double> operator*(const vector<double>& v, matrix<double> m)
  {
    return m.transpose() * v;
  }

  vector<complex> operator*(const vector<complex>& v, matrix<double> m)
  {
    return m.transpose() * v;
  }

  matrix<complex> operator+(const complex& z, matrix<double> m)
  {
    return m+z;
  }

  matrix<complex> operator-(const complex& z, matrix<double> m)
  {
    return -m+z;
  }

  matrix<complex> operator*(const complex& z, matrix<double> m)
  {
    return m*z;
  }
}

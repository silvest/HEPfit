/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <math.h>
#ifndef __GSL_BLAS_H__
#include <gsl/gsl_blas.h>
#endif
#ifndef __GSL_LINALG_H__
#include <gsl/gsl_linalg.h>
#endif
#ifndef __GSL_EIGEN_H__
#include <gsl/gsl_eigen.h>
#endif
#ifndef GSLPP_MATRIX_COMPLEX_H
#include "gslpp_matrix_complex.h"
#endif

namespace gslpp
{
  /** Constructor  */
  matrix<complex>::matrix(const size_t& size_i, const size_t& size_j, const complex& z)
  {
    _matrix = gsl_matrix_complex_alloc(size_i, size_j);
    gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
  }

  matrix<complex>::matrix(const size_t& size_i, const complex& z)
  {
    _matrix = gsl_matrix_complex_alloc(size_i, size_i);
    gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
  }

  matrix<complex>::matrix(const size_t& size_i, const size_t& size_j, const double& a)
  {
    complex z(a);
    _matrix = gsl_matrix_complex_alloc(size_i, size_j);
    gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
  }

  matrix<complex>::matrix(const size_t& size_i, const double& a)
  {
    complex z(a);
    _matrix = gsl_matrix_complex_alloc(size_i, size_i);
    gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
  }
  /** Copy constructor  */
  matrix<complex>::matrix(const matrix<complex>& m)
  {
    _matrix = gsl_matrix_complex_alloc(m.size_i(), m.size_j());
    gsl_matrix_complex_memcpy(_matrix, m.as_gsl_type_ptr());
  }

  matrix<complex>::matrix(const matrix<double>& m)
  {
    size_t i, j, size_i, size_j;
    size_i=m.size_i();
    size_j=m.size_j();
    _matrix = gsl_matrix_complex_alloc(size_i, size_j);
    for(i=0; i<size_i; ++i)
      for(j=0; j<size_j; ++j)
        gsl_matrix_complex_set(_matrix, i, j, m(i,j)+0.*complex::i());
  }

  matrix<complex>::matrix(const vector<complex>& v)
  {
    size_t i, size = v.size();
    complex z(0,0,false);
    _matrix = gsl_matrix_complex_alloc(size, size);
    gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
    for(i=0; i<size; ++i)
        gsl_matrix_complex_set(_matrix, i, i, v(i));
  }

  matrix<complex>::matrix(const vector<double>& v)
  {
    size_t i, size = v.size();
    complex z(0,0,false);
    _matrix = gsl_matrix_complex_alloc(size, size);
    gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
    for(i=0; i<size; ++i)
        gsl_matrix_complex_set(_matrix, i, i, v(i)+0.*complex::i());
  }

  matrix<complex>::matrix(const gsl_matrix_complex& m)
  {
    _matrix = gsl_matrix_complex_alloc(m.size1, m.size2);
    gsl_matrix_complex_memcpy(_matrix, &m);
  }

  matrix<complex>::matrix(const gsl_matrix_complex* m)
  {
    _matrix = gsl_matrix_complex_alloc(m->size1, m->size2);
    gsl_matrix_complex_memcpy(_matrix, m);
  }

  /** Destructor */
  matrix<complex>::~matrix()
  {
    gsl_matrix_complex_free(_matrix);
  }

  /** Get i-th element */
  const complex matrix<complex>::operator()(const size_t& i, const size_t& j) const
  {
    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
    return complex(x);
  }

  /** Set i-th element */
//  complex& matrix<complex>::operator()(const size_t& i, const size_t& j)
//  {
//    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
//    return complex(x);
//  }

  /** Assign */
  matrix<complex>& matrix<complex>::operator=(const matrix<complex>& m)
  {
      if(size_i()==m.size_i() && size_j()==m.size_j())
          gsl_matrix_complex_memcpy(_matrix, m.as_gsl_type_ptr());
      else
      {
          std::cout << "\n ** Wrong size assign in matrix<complex> operator =" << std::endl;
          exit(EXIT_FAILURE);
      }
    return *this;
  }

  /** Assign element */
  void matrix<complex>::assign(const size_t& i, const size_t& j, const complex& z)
  {
    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
    *x = z.as_gsl_type();
  }

  void matrix<complex>::assign(const size_t& i, const size_t& j,const double& a)
  {
    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
    *x = complex(a).as_gsl_type();
  }
  
  void matrix<complex>::assignre(const size_t& i, const size_t& j,const double& a)
  {
    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
    GSL_SET_REAL(x,a);
  }
  
  void matrix<complex>::assignim(const size_t& i, const size_t& j,const double& a)
  {
    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
    GSL_SET_IMAG(x,a);
  }
  
  /** Assign submatrix */
  void matrix<complex>::assign(const size_t& i, const size_t& j, const matrix<complex>& z)
  {
      size_t ki,kj;
      gsl_complex *x;
      if(i+z.size_i() <= size_i() && j+z.size_j() <= size_j())
          for(ki=i;ki<i+z.size_i();ki++)
              for(kj=j;kj<j+z.size_j();kj++)
              {
                  x = gsl_matrix_complex_ptr(_matrix, ki, kj);
                  *x = z(ki-i,kj-j).as_gsl_type();
              }
      else
      {
          std::cout << "\n ** Wrong size assign in matrix<complex> assign submatrix" << std::endl;
          exit(EXIT_FAILURE);
      }
  }

  void matrix<complex>::assign(const size_t& i, const size_t& j, const matrix<double>& a)
  {
      matrix<complex> z(a);
      assign(i,j,z);
  }

  /** Get matrix size */
  size_t matrix<complex>::size_i() const
  {
    return _matrix->size1;
  }

  size_t matrix<complex>::size_j() const
  {
    return _matrix->size2;
  }

  /** Identity matrix */
  matrix<complex> matrix<complex>::Id(size_t size)
  {
    return matrix<double>::Id(size) * (1.+0.*complex::i());
  }

  /** Transpose matrix */
  matrix<complex> matrix<complex>::transpose() const
  {
    matrix<complex> m1(size_j(), size_i(), 0.);
    if (gsl_matrix_complex_transpose_memcpy(m1.as_gsl_type_ptr(), _matrix))
    {
      std::cout << "\n ** Error in matrix<complex> transpose" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Hermitean conjugate matrix */
  matrix<complex> matrix<complex>::hconjugate() const
  {
    matrix<complex> m1(*this);
    gsl_complex z1, z2;

    GSL_SET_COMPLEX(&z1,1.,0.);
    GSL_SET_COMPLEX(&z2,0.,0.);

    if(gsl_blas_zgemm(CblasConjTrans,CblasNoTrans, z1, _matrix,
       Id(size_j()).as_gsl_type_ptr(), z2, m1.as_gsl_type_ptr()) != 0)
    {
      std::cout << "\n ** Error in matrix<complex> conjugate" << std::endl;
      exit(EXIT_FAILURE);
    }

    return m1;
  }

  /** Inverse matrix */
  matrix<complex> matrix<complex>::inverse() const
  {
    matrix<complex> m1(_matrix);
    matrix<complex> m2(_matrix);
    int signum;
    gsl_permutation *p;

    if (size_j() != size_i())
    {
      std::cout << "\n ** Size mismatch in matrix<complex> inverse" << std::endl;
      exit(EXIT_FAILURE);
    }

    if ((p = gsl_permutation_alloc(size_i())) == NULL)
    {
      std::cout << "\n ** Error in matrix<complex> inverse" << std::endl;
      exit(EXIT_FAILURE);
    }

    if(gsl_linalg_complex_LU_decomp (m1.as_gsl_type_ptr(), p, &signum))
    {
      std::cout << "\n ** Error in matrix<complex> inverse" << std::endl;
      gsl_permutation_free(p);
      exit(EXIT_FAILURE);
    }

    if(gsl_linalg_complex_LU_invert(m1.as_gsl_type_ptr(), p,
       m2.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> inverse" << std::endl;
      gsl_permutation_free(p);
      exit(EXIT_FAILURE);
    }
    gsl_permutation_free(p);
    return m2;
  }
  
      /** Get matrix of real parts */
    matrix<double> matrix<complex>::real() const 
    {
        matrix<double> res(size_i(), size_j());
        for(size_t i = 0; i < size_i(); i++)  {
            gsl_vector_complex_view tmp = gsl_matrix_complex_row(_matrix, i);
            gsl_vector_view tmpd = gsl_vector_complex_real(&tmp.vector);
            gsl_matrix_set_row(res.as_gsl_type_ptr(),i,&tmpd.vector);
        }
       return res;
    }

      /** Get matrix of imaginary parts */
    matrix<double> matrix<complex>::imag() const 
    {
        matrix<double> res(size_i(), size_j());
        for(size_t i = 0; i < size_i(); i++)  {
            gsl_vector_complex_view tmp = gsl_matrix_complex_row(_matrix, i);
            gsl_vector_view tmpd = gsl_vector_complex_imag(&tmp.vector);
            gsl_matrix_set_row(res.as_gsl_type_ptr(),i,&tmpd.vector);
        }
       return res;
    }
  
  void matrix<complex>::eigensystem(matrix<complex> &U, vector<double> &S) const
  {
      matrix<complex> m(*this);

      gsl_eigen_hermv_workspace *ws;

      ws = gsl_eigen_hermv_alloc(size_i());

      gsl_eigen_hermv(m.as_gsl_type_ptr(), S.as_gsl_type_ptr(), U.as_gsl_type_ptr(), ws);

      gsl_eigen_hermv_sort(S.as_gsl_type_ptr(), U.as_gsl_type_ptr(),
             GSL_EIGEN_SORT_VAL_ASC);

      gsl_eigen_hermv_free(ws);
  }

  void matrix<complex>::singularvalue(matrix<complex> &U, matrix<complex> &V,
          vector<double> &S) const
  {
     size_t i;
     matrix<complex> m(*this);
     vector<complex> v1(size_i(),0.);

     m = (*this)*hconjugate();
     m.eigensystem(U, S);
     m = hconjugate()*(*this);
     m.eigensystem(V, S);
     m = U.hconjugate()*(*this)*V;
     for(i=0; i<m.size_i(); i++)
     {
         v1.assign(i,complex(1.,-m(i,i).arg(),true));
         S(i)=m(i,i).abs();
     }
     V=V*matrix<complex>(v1);
  }

  /** Conversion  */
  gsl_matrix_complex* matrix<complex>::as_gsl_type_ptr() const
  {
    return _matrix;
  }

  gsl_matrix_complex& matrix<complex>::as_gsl_type()
  {
    return const_cast<gsl_matrix_complex&>(*_matrix);
  }

  const gsl_matrix_complex& matrix<complex>::as_gsl_type() const
  {
    return const_cast<gsl_matrix_complex&>(*_matrix);
  }

//  bool matrix<complex>::is_equal(const matrix<complex>& m1, const matrix<complex>& m2){
//      return gsl_matrix_complex_equal(m1.as_gsl_type_ptr(),m2.as_gsl_type_ptr());
//  }
  
  /** Unary minus */
  matrix<complex> matrix<complex>::operator-() const
  {
    matrix<complex> m1(_matrix);
    gsl_complex z1;
    GSL_SET_COMPLEX(&z1,-1.,0.);
    if (gsl_matrix_complex_scale(m1.as_gsl_type_ptr(), z1))
    {
      std::cout << "\n ** Error in matrix<complex> unary -" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Addition operator (matrix) */
  matrix<complex> matrix<complex>::operator+(const matrix<complex>& m) const
  {
    matrix<complex> m1(_matrix);
    if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> +" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Subtraction operator (matrix) */
  matrix<complex> matrix<complex>::operator-(const matrix<complex>& m) const 
  {
    matrix<complex> m1(_matrix);
    matrix<complex> m2 = -m;
    if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m2.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> +" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Multiplication operator (matrix complex) */
  matrix<complex> matrix<complex>::operator*(const matrix<complex>& m) const 
  {
    matrix<complex> m1(_matrix);
    gsl_complex z1, z2;

    if (size_j() != m.size_i())
    {
      std::cout << "\n ** Size mismatch in matrix<complex> *" << std::endl;
      exit(EXIT_FAILURE);
    }
    GSL_SET_COMPLEX(&z1,1.,0.);
    GSL_SET_COMPLEX(&z2,0.,0.);
    if(gsl_blas_zgemm(CblasNoTrans,CblasNoTrans, z1, _matrix,
       m.as_gsl_type_ptr(), z2, m1.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> *" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

 /** Addition operator (matrix) */
  matrix<complex> matrix<complex>::operator+(const matrix<double>& m) const 
  {
    matrix<complex> m1(_matrix);
    matrix<complex> m2(m);
    if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m2.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> +" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  /** Subtraction operator (matrix) */
  matrix<complex> matrix<complex>::operator-(const matrix<double>& m) const
  {
    matrix<complex> m1(_matrix);
    matrix<complex> m2(-m);
    if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m2.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> +" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Multiplication operator (matrix double) */
  matrix<complex> matrix<complex>::operator*(const matrix<double>& m) const
  {
    matrix<complex> m1(m);
    return *this * m1;
  }

  /** Multiplication operator (vector complex) */
  vector<complex> matrix<complex>::operator*(const vector<complex>& v) const
  {
    gsl_complex z1, z2;
    vector<complex> v1(size_i(),0.);

    if (size_j() != v.size())
    {
      std::cout << "\n ** Size mismatch in matrix<complex> * (vector complex)"
          << std::endl;
      exit(EXIT_FAILURE);
    }
    GSL_SET_COMPLEX(&z1, 1., 0.);
    GSL_SET_COMPLEX(&z2, 0., 0.);
    if(gsl_blas_zgemv(CblasNoTrans, z1, _matrix, v.as_gsl_type_ptr(),
       z2, v1.as_gsl_type_ptr()))
    {
      std::cout << "\n ** Error in matrix<complex> * (vector complex)"
          << std::endl;
      exit(EXIT_FAILURE);
    }
    return v1;
  }

  /** Multiplication operator (vector double) */
  vector<complex> matrix<complex>::operator*(const vector<double>& v) const
  {
    vector<complex> v1(v);
    return *this * v1;
  }

  /** Addition assignment (matrix) */
  matrix<complex>& matrix<complex>::operator+=(const matrix<complex>& m)
  {
    *this = *this + m;
    return *this;
  }
  /** Subtraction assignment (matrix) */
  matrix<complex>& matrix<complex>::operator-=(const matrix<complex>& m)
  {
    *this = *this - m;
    return *this;
  }

  /** Multiplication assignment (matrix) */
  matrix<complex>& matrix<complex>::operator*=(const matrix<complex>& m)
  {
    if (!((size_i() == size_j()) && (size_i() == m.size_i()) &&
           (size_j() == m.size_j())))
    {
      std::cout << "\n ** Size mismatch in matrix<complex> *= (matrix)"
          << std::endl;
      exit(EXIT_FAILURE);
    }
    *this = *this * m;
    return *this;
  }

  /** Addition operator (complex) */
  matrix<complex> matrix<complex>::operator+(const complex& z) const
  {
    matrix<complex> m1(_matrix);
    if (gsl_matrix_complex_add_constant(_matrix, z.as_gsl_type()))
    {
      std::cout << "\n ** Error in matrix<complex> + (complex)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Subtraction operator (complex) */
  matrix<complex> matrix<complex>::operator-(const complex& z) const
  {
    matrix<complex> m1(_matrix);
    if (gsl_matrix_complex_add_constant(_matrix, (-z).as_gsl_type()))
    {
      std::cout << "\n ** Error in matrix<complex> - (complex)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Multiplication operator (complex) */
  matrix<complex> matrix<complex>::operator*(const complex& z) const
  {
    matrix<complex> m1(_matrix);
    if (gsl_matrix_complex_scale(m1.as_gsl_type_ptr(), z.as_gsl_type()))
    {
      std::cout << "\n ** Error in matrix<complex> * (complex)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }

  /** Division operator (complex) */
  matrix<complex> matrix<complex>::operator/(const complex& z) const
  {
    matrix<complex> m1(_matrix);
    if (gsl_matrix_complex_scale(m1.as_gsl_type_ptr(), z.inverse().as_gsl_type()))
    {
      std::cout << "\n ** Error in matrix<complex> / (complex)" << std::endl;
      exit(EXIT_FAILURE);
    }
    return m1;
  }
  
  /** Addition assignment (complex) */
  matrix<complex>& matrix<complex>::operator+=(const complex& z)
  {
    *this = *this + z;
    return *this;
  }
  /** Subtraction assignment (complex) */
  matrix<complex>& matrix<complex>::operator-=(const complex& z)
  {
    *this = *this - z;
    return *this;
  }
  /** Multiplication assignment (complex) */
  matrix<complex>& matrix<complex>::operator*=(const complex& z)
  {
    *this = *this * z;
    return *this;
  }
  /** Division assignment (complex) */
  matrix<complex>& matrix<complex>::operator/=(const complex& z)
  {
    *this = *this / z;
    return *this;
  }

  /** Addition operator (double) */
  matrix<complex> matrix<complex>::operator+(const double& a) const
  {
    complex z(a);
    return *this + z;
  }

  /** Subtraction assignment (double) */
  matrix<complex> matrix<complex>::operator-(const double& a) const
  {
    complex z(a);
    return *this - z;
  }

  /** Multiplication operator (double) */
  matrix<complex> matrix<complex>::operator*(const double& a) const
  {
    complex z(a);
    return *this * z;
  }

  /** Division operator (double) */
  matrix<complex> matrix<complex>::operator/(const double& a) const
  {
    complex z(a);
    return *this / z;
  }
  /** Addition assignment (double) */
  matrix<complex>& matrix<complex>::operator+=(const double& a)
  {
    *this = *this + a;
    return *this;
  }
  /** Subtraction assignment (double) */
  matrix<complex>& matrix<complex>::operator-=(const double& a)
  {
    *this = *this - a;
    return *this;
  }
  /** Multiplication assignment (double) */
  matrix<complex>& matrix<complex>::operator*=(const double& a)
  {
    *this = *this * a;
    return *this;
  }
  /** Division assignment (double) */
  matrix<complex>& matrix<complex>::operator/=(const double& a)
  {
    *this = *this / a;
    return *this;
  }
  /** friend functions */
  std::ostream& operator<<(std::ostream& output, const matrix<complex>& m)
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

  matrix<complex> operator+(const matrix<double> m1, const matrix<complex> m2) {
      return m2 + m1;
  }

  matrix<complex> operator-(const matrix<double> m1, const matrix<complex> m2) {
      return -m2 + m1;
  }

  matrix<complex> operator+(const complex& z, const matrix<complex> m)
  {
    return m + z;
  }

  matrix<complex> operator-(const complex& z, const matrix<complex> m)
  {
    return -m + z;
  }

  matrix<complex> operator*(const complex& z, const matrix<complex> m)
  {
    return m * z;
  }

  vector<complex> operator*(const vector<complex>& v, const matrix<complex> m) 
  {
    return m.transpose() * v;
  }

  vector<complex> operator*(const vector<double>& v, const matrix<complex> m)
  {
    return m.transpose() * v;
  }

  matrix<complex> operator+(const double& a, const matrix<complex> m)
  {
    return m + a;
  }

  matrix<complex> operator-(const double& a, const matrix<complex> m)
  {
    return -m + a;
  }

  matrix<complex> operator*(const double& a, const matrix<complex> m)
  {
    return m * a;
  }
}

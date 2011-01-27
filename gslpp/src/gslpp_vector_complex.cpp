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
#include <stdlib.h>
#include <iostream>
#ifndef __GSL_BLAS_H__
#include <gsl/gsl_blas.h>
#endif
#ifndef GSLPP_VECTOR_COMPLEX_H
#include "gslpp_vector_complex.h"
#endif
#ifndef GSLPP_VECTOR_DOUBLE_H
#include "gslpp_vector_double.h"
#endif

namespace gslpp
{
  /** Constructor  */
  vector<complex>::vector(const size_t& size, const complex& z)
  {
    _vector = gsl_vector_complex_alloc(size);
    gsl_vector_complex_set_all(_vector, z.as_gsl_type());
  }

  vector<complex>::vector(const size_t& size, const double& a)
  {
    complex z(a);
    _vector = gsl_vector_complex_alloc(size);
    gsl_vector_complex_set_all(_vector, z.as_gsl_type());
  }
  /** Copy constructor  */
  vector<complex>::vector(const vector<complex>& v)
  {
    _vector = gsl_vector_complex_alloc(v.size());
    gsl_vector_complex_memcpy(_vector, v.as_gsl_type_ptr());
  }

  vector<complex>::vector(const vector<double>& v)
  {
    size_t i,n;
    n=v.size();
    _vector = gsl_vector_complex_alloc(n);
    for(i=0; i<n; i++)
      gsl_vector_complex_set(_vector, i, v(i)+0.*complex::i());
  }

  vector<complex>::vector(const gsl_vector_complex& v)
  {
    _vector = gsl_vector_complex_alloc(v.size);
    gsl_vector_complex_memcpy(_vector, &v);
  }

  vector<complex>::vector(const gsl_vector_complex* v)
  {
    _vector = gsl_vector_complex_alloc(v->size);
    gsl_vector_complex_memcpy(_vector, v);
  }

  /** Destructor */
  vector<complex>::~vector()
  {
    gsl_vector_complex_free(_vector);
  }

  /** Get i-th element */
  const complex vector<complex>::operator()(const size_t& i) const
  {
    const gsl_complex *x = gsl_vector_complex_const_ptr(_vector, i);
    return complex(x);
  }

  /** Set i-th element */
//  complex& vector<complex>::operator()(const size_t& i)
//  {
//    gsl_complex *x = gsl_vector_complex_ptr(_vector, i);
//    return complex(x);
//  }

  /** Assign */
  vector<complex>& vector<complex>::operator=(const vector<complex>& v)
  {
    gsl_vector_complex_memcpy(_vector, v.as_gsl_type_ptr());
    return *this;
  }

  /** Assign element */
  void vector<complex>::assign(const size_t& i, const complex& z)
  {
    gsl_complex *x = gsl_vector_complex_ptr(_vector, i);
    *x = z.as_gsl_type();
  }

  void vector<complex>::assign(const size_t& i, const double& a)
  {
    gsl_complex *x = gsl_vector_complex_ptr(_vector, i);
    *x = complex(a).as_gsl_type();
  }

  /** Get vector size */
  size_t vector<complex>::size() const
  {
    return _vector->size;
  }

  /** Get Euclidean norm */
  double vector<complex>::mod() const
  {
    return gsl_blas_dznrm2(_vector);
  }

  /** Conversion  */
  gsl_vector_complex* vector<complex>::as_gsl_type_ptr() const
  {
    return _vector;
  }

  gsl_vector_complex& vector<complex>::as_gsl_type()
  {
    return const_cast<gsl_vector_complex&>(*_vector);
  }

  const gsl_vector_complex& vector<complex>::as_gsl_type() const
  {
    return const_cast<gsl_vector_complex&>(*_vector);
  }

  /** Unary minus */
  vector<complex> vector<complex>::operator-() const
  {
    vector<complex> v1(_vector);
    gsl_blas_zdscal(-1.,v1.as_gsl_type_ptr());
    return v1;
  }

  /** Addition operator (vector) */
  vector<complex> vector<complex>::operator+(const vector<complex>& v)
  {
    vector<complex> v1(_vector);
    gsl_complex z1;
    GSL_SET_COMPLEX(&z1,1.,0.);
    if (gsl_blas_zaxpy(z1, v.as_gsl_type_ptr(),v1.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in vector<complex> +" << std::endl;
      exit(EXIT_FAILURE);
    }
    return v1;
  }
  /** Subtraction operator (vector) */
  vector<complex> vector<complex>::operator-(const vector<complex>& v)
  {
    vector<complex> v1(_vector);
    gsl_complex z1;
    GSL_SET_COMPLEX(&z1,-1.,0.);
    if (gsl_blas_zaxpy(z1, v.as_gsl_type_ptr(), v1.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in vector<complex> -" << std::endl;
      exit(EXIT_FAILURE);
    }
    return v1;
  }
  /** Scalar product operator (vector) */
  complex vector<complex>::operator*(const vector<complex>& v)
  {
    complex z1(0.,0.,false);
    vector<complex> v1(_vector);
    if(gsl_blas_zdotu(v1.as_gsl_type_ptr(), v.as_gsl_type_ptr(), z1.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in vector<complex> *" << std::endl;
      exit(EXIT_FAILURE);
    }
    return z1;
  }
  /** Vector product operator (vector) */
//   vector<complex> vector<complex>::operator^(const vector<complex>& v)
//   {
//     std::cout << "\n To be implemented" << std::endl;
//     exit(EXIT_FAILURE);
//   }

  /** Addition assignment (vector) */
  vector<complex>& vector<complex>::operator+=(const vector<complex>& v)
  {
/*    gsl_complex z1;
    GSL_SET_COMPLEX(&z1,1.,0.);
    if (gsl_blas_zaxpy(z1, v.as_gsl_type_ptr(), _vector))
    {
      std::cout << "\n Error in vector<complex> +=" << std::endl;
      exit(EXIT_FAILURE);
    }
    return *this;*/
    *this = *this + v;
    return *this;
  }
  /** Subtraction assignment (vector) */
  vector<complex>& vector<complex>::operator-=(const vector<complex>& v)
  {
    *this = *this - v;
    return *this;
  }

  /** Addition operator (complex) */
  vector<complex> vector<complex>::operator+(const complex& z)
  {
    vector<complex> v1(_vector);
    gsl_complex z1;
    gsl_vector_complex *v2;
    GSL_SET_COMPLEX(&z1,1.,0.);
    v2 = gsl_vector_complex_alloc(v1.size());
    gsl_vector_complex_set_all(v2,z1);
    if (gsl_blas_zaxpy(z, v2, v1.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in vector<complex> + (double)" << std::endl;
      exit(EXIT_FAILURE);
    }
    gsl_vector_complex_free(v2);
    return v1;
  }

  /** Subtraction assignment (complex) */
  vector<complex> vector<complex>::operator-(const complex& z)
  {
    vector<complex> v1(_vector);
    gsl_complex z1;
    gsl_vector_complex *v2;
    GSL_SET_COMPLEX(&z1,-1.,0.);
    v2 = gsl_vector_complex_alloc(v1.size());
    gsl_vector_complex_set_all(v2,z1);
    if (gsl_blas_zaxpy(z, v2, v1.as_gsl_type_ptr()))
    {
      std::cout << "\n Error in vector<complex> - (double)" << std::endl;
      exit(EXIT_FAILURE);
    }
    gsl_vector_complex_free(v2);
    return v1;
  }

  /** Multiplication operator (complex) */
  vector<complex> vector<complex>::operator*(const complex& z)
  {
    vector<complex> v1(_vector);
    gsl_blas_zscal(z.as_gsl_type(), v1.as_gsl_type_ptr());
    return v1;
  }

  /** Division operator (complex) */
  vector<complex> vector<complex>::operator/(const complex& z)
  {
    vector<complex> v1(_vector);
    gsl_blas_zscal(z.inverse().as_gsl_type(), v1.as_gsl_type_ptr());
    return v1;
  }
  /** Addition assignment (complex) */
  vector<complex>& vector<complex>::operator+=(const complex& z)
  {
    *this = *this + z;
    return *this;
  }
  /** Subtraction assignment (complex) */
  vector<complex>& vector<complex>::operator-=(const complex& z)
  {
    *this = *this - z;
    return *this;
  }
  /** Multiplication assignment (complex) */
  vector<complex>& vector<complex>::operator*=(const complex& z)
  {
    *this = *this * z;
    return *this;
  }
  /** Division assignment (complex) */
  vector<complex>& vector<complex>::operator/=(const complex& z)
  {
    *this = *this / z;
    return *this;
  }

  /** Addition operator (double) */
  vector<complex> vector<complex>::operator+(const double& a)
  {
    complex z(a);
    return *this + z;
  }

  /** Subtraction assignment (double) */
  vector<complex> vector<complex>::operator-(const double& a)
  {
    complex z(a);
    return *this - z;
  }

  /** Multiplication operator (double) */
  vector<complex> vector<complex>::operator*(const double& a)
  {
    complex z(a);
    return *this * z;
  }

  /** Division operator (double) */
  vector<complex> vector<complex>::operator/(const double& a)
  {
    complex z(a);
    return *this / z;
  }
  /** Addition assignment (double) */
  vector<complex>& vector<complex>::operator+=(const double& a)
  {
    *this = *this + a;
    return *this;
  }
  /** Subtraction assignment (double) */
  vector<complex>& vector<complex>::operator-=(const double& a)
  {
    *this = *this - a;
    return *this;
  }
  /** Multiplication assignment (double) */
  vector<complex>& vector<complex>::operator*=(const double& a)
  {
    *this = *this * a;
    return *this;
  }
  /** Division assignment (double) */
  vector<complex>& vector<complex>::operator/=(const double& a)
  {
    *this = *this / a;
    return *this;
  }
  /** friend functions */
  std::ostream& operator<<(std::ostream& output, const vector<complex>& v)
  {
    size_t i;
    output << "(";
    for (i=0; i<v.size()-1; i++)
      output << v(i) << ",";
    output << v(i) << ")";
    return output;
  }
  /** @{
   * @name Operations on vector<complex>
   */
  /** Add a complex number to a complex vector
   * @ingroup vector
   * @param z Complex number
   * @param v Complex vector
   * @return @f$ z + v @f$
   */
  vector<complex> operator+(const complex& z, vector<complex> v)
  {
    return v+z;
  }

  /** Subtract a complex number from a complex vector
   * @ingroup vector
   * @param z Complex number
   * @param v Complex vector
   * @return @f$ z - v @f$
   */
  vector<complex> operator-(const complex& z, vector<complex> v)
  {
    return -v+z;
  }

  /** Multiply a complex number by a complex vector
   * @ingroup vector
   * @param z Complex number
   * @param v Complex vector
   * @return @f$ z*v @f$
   */
  vector<complex> operator*(const complex& z, vector<complex> v)
  {
    return v*z;
  }

  /** Add a real number to a complex vector
   * @ingroup vector
   * @param a Real number
   * @param v Complex vector
   * @return @f$ a + v @f$
   */
  vector<complex> operator+(const double& a, vector<complex> v)
  {
    return v+a;
  }

  /** Subtract a complex vector from a real number
   * @ingroup vector
   * @param a Real number
   * @param v Complex vector
   * @return @f$ a - v @f$
   */
  vector<complex> operator-(const double& a, vector<complex> v)
  {
    return -v+a;
  }

  /** Multiply a real number by a complex vector
   * @ingroup vector
   * @param a Real number
   * @param v Complex vector
   * @return @f$ z*v @f$
   */
  vector<complex> operator*(const double& a, vector<complex> v)
  {
    return v*a;
  }
  /** @}
   */
}

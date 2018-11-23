/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#ifndef __GSL_BLAS_H__
#include <gsl/gsl_blas.h>
#endif
#ifndef GSLPP_VECTOR_DOUBLE_H
#include "gslpp_vector_double.h"
#endif
#ifndef GSLPP_VECTOR_COMPLEX_H
#include "gslpp_vector_complex.h"
#endif

namespace gslpp
{
  /** Constructor  */
  vector<double>::vector(const size_t& size, const double& a)
  {
    _vector = gsl_vector_alloc(size);
    gsl_vector_set_all(_vector, a);
  }

  vector<double>::vector(const size_t& size)
  {
    _vector = gsl_vector_alloc(size);
    gsl_vector_set_all(_vector, 0.);
  }

  /** Copy constructor */
  vector<double>::vector(const vector<double>& v)
  {
    _vector = gsl_vector_alloc(v.size());
    gsl_vector_memcpy(_vector, v.as_gsl_type_ptr());
  }

  vector<double>::vector(const gsl_vector& v)
  {
    _vector = gsl_vector_alloc(v.size);
    gsl_vector_memcpy(_vector, &v);
  }

  vector<double>::vector(const gsl_vector* v)
  {
    _vector = gsl_vector_alloc(v->size);
    gsl_vector_memcpy(_vector, v);
  }

  /** Destructor */
  vector<double>::~vector()
  {
    gsl_vector_free(_vector);
  }

  /** Get i-th element */
  const double& vector<double>::operator()(const size_t& i) const
    {
      const double *x = gsl_vector_const_ptr(_vector, i);
      return *x;
    }

  /** Set i-th element */
  double& vector<double>::operator()(const size_t& i)
  {
    double *x = gsl_vector_ptr(_vector, i);
    return *x;
  }

  /** Assign */
  void vector<double>::reset()
    {
        gsl_vector_set_zero(_vector);
    }
    

  vector<double>& vector<double>::operator=(const vector<double>& v)
  {
    gsl_vector_memcpy(_vector, v.as_gsl_type_ptr());
    return *this;
  }

  /** Get vector size */
  size_t vector<double>::size() const
  {
    return _vector->size;
  }

  /** Get Euclidean norm */
  double vector<double>::mod() const
  {
    return gsl_blas_dnrm2(_vector);
  }
  
  /** Get max element */
  double vector<double>::max() const
  {
    return gsl_vector_max(_vector);
  }
  
  /** Get min element */
  double vector<double>::min() const
  {
    return gsl_vector_min(_vector);
  }

  /** Conversion */
  gsl_vector* vector<double>::as_gsl_type_ptr() const
  {
    return _vector;
  }

  gsl_vector& vector<double>::as_gsl_type()
  {
    return const_cast<gsl_vector&>(*_vector);
  }

  const gsl_vector& vector<double>::as_gsl_type() const
  {
    return const_cast<gsl_vector&>(*_vector);
  }

  /** Unary minus */
  vector<double> vector<double>::operator-() const
  {
    vector<double> v1(_vector);
    if (gsl_vector_scale(v1.as_gsl_type_ptr(), -1.))
      {
        std::cout << "\n Error in vector<double> unary -" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }

  /** Addition operator (vector) */
  vector<double> vector<double>::operator+(const vector<double>& v) const
  {
    vector<double> v1(_vector);
    if (gsl_vector_add(v1.as_gsl_type_ptr(), v.as_gsl_type_ptr()))
      {
        std::cout << "\n Error in vector<double> +" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }
  vector<complex> vector<double>::operator+(const vector<complex>& v) const
  {
    vector<complex> v1(_vector);
    return v1 + v;
  }
  /** Subtraction operator (vector) */
  vector<double> vector<double>::operator-(const vector<double>& v) const
  {
    vector<double> v1(_vector);
    if (gsl_vector_sub(v1.as_gsl_type_ptr(), v.as_gsl_type_ptr()))
      {
        std::cout << "\n Error in vector<double> -" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }
  vector<complex> vector<double>::operator-(const vector<complex>& v) const
  {
      vector<complex> v1(_vector);
      return v1-v;
  }
  /** Scalar product operator (vector) */
  double vector<double>::operator*(const vector<double>& v) const
  {
    double a=0.;
    vector<double> v1(_vector);
    if(gsl_blas_ddot(v1.as_gsl_type_ptr(), v.as_gsl_type_ptr(), &a))
      {
        std::cout << "\n Error in vector<double> *" << std::endl;
        exit(EXIT_FAILURE);
      }
    return a;
  }
  complex vector<double>::operator*(const vector<complex>& v) const
  {
      vector<complex> v1(_vector);
      return v*(*this);
  }

  /** Vector product operator (vector) */
//   vector<double> vector<double>::operator^(const vector<double>& v)
//   {
//     std::cout << "\n To be implemented" << std::endl;
//     exit(EXIT_FAILURE);
//   }

  /** Addition assignment (vector) */
  vector<double>& vector<double>::operator+=(const vector<double>& v)
  {
    *this = *this + v;
    return *this;
  }
  /** Subtraction assignment (vector) */
  vector<double>& vector<double>::operator-=(const vector<double>& v)
  {
    *this = *this - v;
    return *this;
  }

  /** Addition operator (double) */
  vector<double> vector<double>::operator+(const double& a) const
  {
    vector<double> v1(_vector);
    if (gsl_vector_add_constant(v1.as_gsl_type_ptr(), a))
      {
        std::cout << "\n Error in vector<double> + (double)" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }
  /** Subtraction operator (double) */
  vector<double> vector<double>::operator-(const double& a) const
  {
    vector<double> v1(_vector);
    if (gsl_vector_add_constant(v1.as_gsl_type_ptr(), -a))
      {
        std::cout << "\n Error in vector<double> - (double)" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }
  /** Multiplication operator (double) */
  vector<double> vector<double>::operator*(const double& a) const
  {
    vector<double> v1(_vector);
    if (gsl_vector_scale(v1.as_gsl_type_ptr(), a))
      {
        std::cout << "\n Error in vector<double> * (double)" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }
  /** Division operator (double) */
  vector<double> vector<double>::operator/(const double& a) const
  {
    vector<double> v1(_vector);
    if (gsl_vector_scale(v1.as_gsl_type_ptr(), 1./a))
      {
        std::cout << "\n Error in vector<double> / (double)" << std::endl;
        exit(EXIT_FAILURE);
      }
    return v1;
  }
  /** Addition assignment (double) */
  vector<double>& vector<double>::operator+=(const double& a)
  {
    *this = *this + a;
    return *this;
  }
  /** Subtraction assignment (double) */
  vector<double>& vector<double>::operator-=(const double& a)
  {
    *this = *this - a;
    return *this;
  }
  /** Multiplication assignment (double) */
  vector<double>& vector<double>::operator*=(const double& a)
  {
    *this = (*this) * a;
    return *this;
  }
  /** Division assignment (double) */
  vector<double>& vector<double>::operator/=(const double& a)
  {
    *this = *this / a;
    return *this;
  }
  /** Addition operator (complex) */
  vector<complex> vector<double>::operator+(const complex& z) const
  {
    vector<complex> v1(*this);
    return v1+z;
  }
  /** Subtraction operator (complex) */
  vector<complex> vector<double>::operator-(const complex& z) const
  {
    vector<complex> v1(*this);
    return v1-z;
  }
  /** Multiplication operator (complex) */
  vector<complex> vector<double>::operator*(const complex& z) const
  {
    vector<complex> v1(*this);
    return v1*z;
  }
  /** Division operator (complex) */
  vector<complex> vector<double>::operator/(const complex& z) const
  {
    vector<complex> v1(*this);
    return v1/z;
  }
  /** Comparison == (vector) */
  bool vector<double>::operator==(const vector<double>& a) const
  {
    if(a.size() != size())
    {
        std::cout << "\n Error in vector<double>::operator== (vector): cannot compare vectors of different size" << std::endl;
        exit(EXIT_FAILURE);
    }
    for(size_t i = 0; i < size(); i++)
        if(a(i) != (*this)(i)) return(false);
    return(true);
  }

  /** friend functions */
  std::ostream& operator<<(std::ostream& output, const vector<double>& v)
  {
    size_t i;
    output << "(";
    for (i=0; i<v.size()-1; i++)
      output << v(i) << ",";
    output << v(i) << ")";
    return output;
  }

  vector<double> operator+(const double& a, vector<double> v)
  {
    return v+a;
  }

  vector<double> operator-(const double& a, vector<double> v)
  {
    return -v+a;
  }

  vector<double> operator*(const double& a, vector<double> v)
  {
    return v*a;
  }

  vector<complex> operator+(const complex& z, vector<double> v)
  {
    return v+z;
  }

  vector<complex> operator-(const complex& z, vector<double> v)
  {
    return -v+z;
  }

  vector<complex> operator*(const complex& z, vector<double> v)
  {
    return v*z;
  }
}

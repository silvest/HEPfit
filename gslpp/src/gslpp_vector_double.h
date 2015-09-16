/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GSLPP_VECTOR_DOUBLE_H
#define GSLPP_VECTOR_DOUBLE_H
#include <iostream>
#ifndef __GSL_VECTOR_DOUBLE_H__
# include <gsl/gsl_vector_double.h>
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
     * @class vector<double>
     * @ingroup gslpp
     * @brief A class for constructing and defining operations on real vectors.
     * @author HEPfit Collaboration
     * @copyright GNU General Public License
     * @details This class defines some common operations on real vectors
     * using the <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
     */
  template <>
  class vector<double>
  {
    gsl_vector *_vector;

    public:
      /** Constructor */
      vector(const size_t& size);
      vector(const size_t& size, const double& a);
      /** Copy constructor */
      vector(const vector<double>& v);
      vector(const gsl_vector& v);
      vector(const gsl_vector* v);
      /** Destructor */
      ~vector();
      /** Get i-th element */
      double operator()(const size_t& i) const;
      /** Set i-th element */
      double& operator()(const size_t& i);
      /** Assign */
      vector<double>& operator=(const vector<double>& v);
      /** Get vector size */
      size_t size() const;
      /** Get Euclidean norm */
      double mod() const;
      /** Conversion */
      gsl_vector* as_gsl_type_ptr() const;
      gsl_vector& as_gsl_type();
      const gsl_vector& as_gsl_type() const;
      /** Unary minus (vector) */
      vector<double> operator-() const;
      /** Addition operator (vector) */
      vector<double> operator+(const vector<double>& v);
      /** Subtraction operator (vector) */
      vector<double> operator-(const vector<double>& v);
      /** Scalar product operator (vector) */
      double operator*(const vector<double>& v);
      /** Vector product operator */
//       vector<double> operator^(const vector<double>& v);
      /** Addition assignment (vector) */
      vector<double>& operator+=(const vector<double>& v);
      /** Subtraction assignment (vector) */
      vector<double>& operator-=(const vector<double>& v);
      /** Addition operator (double) */
      vector<double> operator+(const double& a);
      /** Subtraction assignment (double) */
      vector<double> operator-(const double& a);
      /** Multiplication operator (double) */
      vector<double> operator*(const double& a);
      /** Division operator (double) */
      vector<double> operator/(const double& a);
      /** Addition assignment (double) */
      vector<double>& operator+=(const double& a);
      /** Subtraction assignment (double) */
      vector<double>& operator-=(const double& a);
      /** Multiplication assignment (double) */
      vector<double>& operator*=(const double& a);
      /** Division assignment (double) */
      vector<double>& operator/=(const double& a);
      /** Addition operator (complex) */
      vector<complex> operator+(const complex& z);
      /** Subtraction assignment (complex) */
      vector<complex> operator-(const complex& z);
      /** Multiplication operator (complex) */
      vector<complex> operator*(const complex& z);
      /** Division operator (complex) */
      vector<complex> operator/(const complex& z);

      /** friend functions */
      friend std::ostream& operator<<(std::ostream& output, const vector<double>& v);
      /** @{
       * @name Operations on vector<double>
       */
      /** Add a real number to a real vector
       * @ingroup vector
       * @param a Real number
       * @param v Real vector
       * @return @f$ a + v @f$
       */
      friend vector<double> operator+(const double& a, vector<double> v);

      /** Subtract a real number from a real vector
       * @ingroup vector
       * @param a Real number
       * @param v Real vector
       * @return @f$ a - v @f$
       */
      friend vector<double> operator-(const double& a, vector<double>  v);

      /** Multiply a real number by real vector
       * @ingroup vector
       * @param a Real number
       * @param v Real vector
       * @return @f$ a*v @f$
       */
      friend vector<double> operator*(const double& a, vector<double> v);

      /** Add a complex number to a real vector
       * @ingroup vector
       * @param z Complex number
       * @param v Real vector
       * @return @f$ z + v @f$
       */
      friend vector<complex> operator+(const complex& z, vector<double> v);

      /** Subtract a complex number from a real vector
       * @ingroup vector
       * @param z Complex number
       * @param v Real vector
       * @return @f$ z - v @f$
       */
      friend vector<complex> operator-(const complex& z, vector<double> v);

      /** Multiply a complex number by a real vector
       * @ingroup vector
       * @param z Complex number
       * @param v Real vector
       * @return @f$ z*v @f$
       */
      friend vector<complex> operator*(const complex& z, vector<double> v);
      /** @}
       */
  };
}

#endif

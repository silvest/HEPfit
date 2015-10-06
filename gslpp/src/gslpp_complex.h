/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#ifndef GSLPP_COMPLEX
#define GSLPP_COMPLEX
#ifndef __GSL_COMPLEX_H__
# include <gsl/gsl_complex.h>
#endif
#ifndef __GSL_COMPLEX_MATH_H__
# include <gsl/gsl_complex_math.h>
#endif

/**
 * @namespace gslpp 
 * @brief Complex number, vector and matrix manipulation
 * using <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
 * @{
 */
namespace gslpp
{
    /**
     * @class complex
     * @ingroup gslpp
     * @brief A class for defining operations on and functions of complex numbers.
     * @author HEPfit Collaboration
     * @copyright GNU General Public License
     * @details  This class defines some common operations on complex variables
     * using the <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
     */
  class complex
  {
      gsl_complex _complex;
    public:
      static const complex& i();
      /**
       * @brief Default constructor for the complex class.
       */
      complex();
      
      /**
       * @brief Overleaded constructor for the complex class.
       * @param[in] real real part of the comlex number
       * @param[in] imag imaginary part of the complex number
       * @param[in] polar boolean switch for specifying polar form of
       * of the complex number (default: false)
       */
      complex(const double& real, const double& imag, bool polar = false);
      /**
       * @brief Copy constructor for the complex class.
       */
      complex(const complex& z);
      /**
       * @brief Conversion constructor for the complex class.
       * Converts a real double to a complex type.
       */
      complex(const double& a);
      /**
       * @brief Conversion constructor for the complex class.
       * Converts a gsl_complex type pointer to a complex type.
       */
      complex(const gsl_complex* z);
      /**
       * @brief Conversion constructor for the complex class.
       * Converts a gsl_complex type reference to a complex type.
       */
      complex(const gsl_complex& z);
      /**
       * @brief Default destructor for the complex class.
       */
      virtual ~complex();
      /**
       * @brief Check if complex number is purely real
       * @return boolean true/false
       */
      bool is_real() const;
      /**
       * @brief Check if complex number is purely imaginary
       * @return boolean true/false
       */
      bool is_imag() const;
      /** Assign */
      void assign(const double& real, const double& imag, bool polar);
      /** Set the real part */
      const double& real() const;
      /** Set imaginary part */
      const double& imag() const;
      /** Get the real part */
      double& real();
      /** Get the imaginary part */
      double& imag();
      /**
       * @return The argument of a complex number
       */
      double arg() const;
      /**
       * @return The absolute value of a complex number
       */
      double abs() const;
      /**
       * @return The square of the absolute value of a complex number
       */
      double abs2() const;
      /**
       * @return The logarithm of the absolute value of a complex number
       */
      double log_of_abs() const;
      /**
       * @brief Assignment operator for a complex variable of complex type
       */
      complex& operator=(const complex& z);
      /**
       * @brief Assignment operator for a double variable to complex type
       */
      complex& operator=(const double& x);
      /**
       * @brief Equivalence operator between two complex variables
       */
      bool operator==(const complex& z1) const;
      /**
       * @brief Inequivalence operator between two complex variables
       */
      bool operator!=(const complex& z1) const;
      /**
       * @brief Unary minus operator for a complex number
       */
      complex operator-() const;
      /**
       * @brief Addition operator for a complex number
       */
      complex operator+(const complex& z1) const;
      /**
       * @brief Subtraction operator for a complex number
       */
      complex operator-(const complex& z1) const;
      /**
       * @brief Multiplication operator for a complex number
       */
      complex operator*(const complex& z1) const;
      /**
       * @brief Division operator for a complex number
       */
      complex operator/(const complex& z1) const;
      /**
       * @brief Addition assignment operator for a complex number
       */
      complex& operator+=(const complex& z1);
      /**
       * @brief Subtraction assignment operator for a complex number
       */
      complex& operator-=(const complex& z1);
      /**
       * @brief Muliplication assignment operator for a complex number
       */
      complex& operator*=(const complex& z1);
      /**
       * @brief Division assignment operator for a complex number
       */
      complex& operator/=(const complex& z1);
      /**
       * @brief Addition operator for adding a real number to a complex number
       */
      complex operator+(const double& a) const;
      /**
       * @brief Subtraction operator for subtracting a real number from a complex number
       */
      complex operator-(const double& a) const;
      /**
       * @brief Multiplication operator for multiplying a real number to a complex number
       */
      complex operator*(const double& a) const;
      /**
       * @brief Divsion operator for dividing a complex number by a real number
       */
      complex operator/(const double& a) const;
      /**
       * @brief Addition assignment operator for adding a real number to a complex number
       */
      complex& operator+=(const double& a);
      /**
       * @brief Subtraction assignment operator for subtracting a real number from a complex number
       */
      complex& operator-=(const double& a);
      /**
       * @brief Multiplication assignment operator for mulitplying a real number to a complex number
       */
      complex& operator*=(const double& a);
      /**
       * @brief Division assignment operator for dividing a complex number by a real number
       */
      complex& operator/=(const double& a);
      /**
       * @return The conjugate of a complex number
       */
      complex conjugate() const;
      /**
       * @return The inverse of a complex number
       */
      complex inverse() const;
      /** Conversion  */
      gsl_complex* as_gsl_type_ptr() const;
      gsl_complex& as_gsl_type();
      const gsl_complex& as_gsl_type() const;
      operator gsl_complex& ();
      operator const gsl_complex& () const;

      /** Friend functions */
      /**
        * @ingroup complex
        * @param[in] output output stream
        * @param[in] z Complex number
        * @return formatted output for complex
        */
      friend std::ostream& operator<<(std::ostream& output, const complex& z);
      /** @{
       * @name Operations on complex numbers
       */
      /** Add a real and complex numbers
       * @ingroup complex
       * @param[in] x1 Real number
       * @param[in] z2 Complex number
       * @return @f$ x_1 + z_2 @f$
       */
      friend complex operator+(const double& x1, const complex& z2);

      /** Subtract a real and complex numbers
       * @ingroup complex
       * @param[in] x1 Real number
       * @param[in] z2 Complex number
       * @return @f$ x_1 - z_2 @f$
       */
      friend complex operator-(const double& x1, const complex& z2);

      /** Multiply a real and complex numbers
       * @ingroup complex
       * @param[in] x1 Real number
       * @param[in] z2 Complex number
       * @return @f$ x_1 z_2 @f$
       */
      friend complex operator*(const double& x1, const complex& z2);

      /** Divide a real and complex numbers
       * @ingroup complex
       * @param[in] x1 Real number
       * @param[in] z2 Complex number
       * @return @f$ x_1 / z_2 @f$
       */
      friend complex operator/(const double& x1,const complex& z2);
      /** @}
       */

      friend complex exp(const complex& z);
      friend complex log(const complex& z);
      friend complex log10(const complex& z);
      friend complex log(const complex& z,
                          const complex& b);
      friend complex dilog(const complex& z);
      friend complex sqrt(const complex& z);
      friend complex pow(const complex& z1,
                          const complex& z2);
      friend complex pow(const complex& z,
                          const double  x);
      friend complex sin(const complex& z);
      friend complex cos(const complex& z);
      friend complex tan(const complex& z);
      friend complex sec(const complex& z);
      friend complex csc(const complex& z);
      friend complex cot(const complex& z);
      friend complex arcsin(const complex& z);
      friend complex arccos(const complex& z);
      friend complex arctan(const complex& z);
      friend complex arcsec(const complex& z);
      friend complex arccsc(const complex& z);
      friend complex arccot(const complex& z);
      friend complex sinh(const complex& z);
      friend complex cosh(const complex& z);
      friend complex tanh(const complex& z);
      friend complex sech(const complex& z);
      friend complex csch(const complex& z);
      friend complex coth(const complex& z);
      friend complex arcsinh(const complex& z);
      friend complex arccosh(const complex& z);
      friend complex arctanh(const complex& z);
      friend complex arcsech(const complex& z);
      friend complex arccsch(const complex& z);
      friend complex arccoth(const complex& z);
  };
}
#endif
/** @} */

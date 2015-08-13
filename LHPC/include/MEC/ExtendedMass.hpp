/*
 * ExtendedMass.hpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef EXTENDEDMASS_HPP_
#define EXTENDEDMASS_HPP_

#include <string>
#include "BOLlib/include/BOLlib.hpp"
#include "RunningConstant.hpp"
#include "RunningConstantError.hpp"

namespace LHPC
{
  // this is a class to hold the information about the mass of a particle in
  // the FLHA format.
  class ExtendedMass
  {
  public:
    ExtendedMass();
    ExtendedMass( ExtendedMass const& copySource );
    ~ExtendedMass();

    double
    getMass() const;
    double
    getMinusUncertainty() const;
    double
    getPlusUncertainty() const;
    int
    getScheme() const;
    double
    getScale() const;
    void
    setValues( double const massValue,
               double const minusUncertainty,
               double const plusUncertainty,
               int const schemeType,
               double const evaluationScale );
    void
    setUncertainties( double const minusUncertainty,
                      double const plusUncertainty );


  protected:
    double massValue;
    double minusUncertainty;
    double plusUncertainty;
    int schemeType;
    double evaluationScale;
  };





  inline double
  ExtendedMass::getMass() const
  {
    return massValue;
  }

  inline double
  ExtendedMass::getMinusUncertainty() const
  {
    return minusUncertainty;
  }

  inline double
  ExtendedMass::getPlusUncertainty() const
  {
    return plusUncertainty;
  }

  inline int
  ExtendedMass::getScheme() const
  {
    return schemeType;
  }

  inline double
  ExtendedMass::getScale() const
  {
    return evaluationScale;
  }

  inline void
  ExtendedMass::setValues( double const massValue,
                           double const minusUncertainty,
                           double const plusUncertainty,
                           int const schemeType,
                           double const evaluationScale )
  {
    this->massValue = massValue;
    this->minusUncertainty = minusUncertainty;
    this->plusUncertainty = plusUncertainty;
    this->schemeType = schemeType;
    this->evaluationScale = evaluationScale;
  }

  inline void
  ExtendedMass::setUncertainties( double const minusUncertainty,
                                  double const plusUncertainty )
  {
    this->minusUncertainty = minusUncertainty;
    this->plusUncertainty = plusUncertainty;
  }

}

#endif /* EXTENDEDMASS_HPP_ */

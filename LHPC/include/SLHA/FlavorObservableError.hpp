/*
 * FlavorObservableError.hpp
 *
 *  Created on: Apr 1, 2012 (really!)
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef FLAVOROBSERVABLEERROR_HPP_
#define FLAVOROBSERVABLEERROR_HPP_

#include <string>
#include "BOLlib/include/BOLlib.hpp"
#include "FlavorObservable.hpp"

namespace LHPC
{
  // this is a class to hold the information about the mass of a particle in
  // the FLHA format.
  class FlavorObservableError
  {
  public:
    FlavorObservableError();
    FlavorObservableError( FlavorObservableError const& copySource );
    ~FlavorObservableError();

    double
    getMinusUncertainty() const;
    double
    getPlusUncertainty() const;
    double
    getScale() const;
    int
    getNumberOfDaughterParticles() const;
    std::list< int >&
    getDaughterParticleList();
    std::list< int > const&
    getDaughterParticleList() const;
    void
    setValues( double const minusUncertainty,
               double const plusUncertainty,
               double const evaluationScale,
               std::list< int > const& daughterParticleCodes );
    void
    setFromString( std::string const& valuesString );
    std::string
    getAsString() const;


  protected:
    double minusUncertainty;
    double plusUncertainty;
    double evaluationScale;
    std::list< int > daughterParticleCodes;
  };





  inline double
  FlavorObservableError::getMinusUncertainty() const
  {
    return minusUncertainty;
  }

  inline double
  FlavorObservableError::getPlusUncertainty() const
  {
    return plusUncertainty;
  }

  inline double
  FlavorObservableError::getScale() const
  {
    return evaluationScale;
  }

  inline int
  FlavorObservableError::getNumberOfDaughterParticles() const
  {
    return (int)(daughterParticleCodes.size());
  }

  inline std::list< int >&
  FlavorObservableError::getDaughterParticleList()
  {
    return daughterParticleCodes;
  }

  inline std::list< int > const&
  FlavorObservableError::getDaughterParticleList() const
  {
    return daughterParticleCodes;
  }

  inline void
  FlavorObservableError::setValues(double const minusUncertainty,
                                   double const plusUncertainty,
                                   double const evaluationScale,
                                std::list< int > const& daughterParticleCodes )
  {
    this->minusUncertainty = minusUncertainty;
    this->plusUncertainty = plusUncertainty;
    this->evaluationScale = evaluationScale;
    this->daughterParticleCodes = daughterParticleCodes;
  }

}

#endif /* FLAVOROBSERVABLEERROR_HPP_ */

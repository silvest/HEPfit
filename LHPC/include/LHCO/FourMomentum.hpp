/*
 * FourMomentum.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef FOURMOMENTUM_HPP_
#define FOURMOMENTUM_HPP_

#include <list>
#include <vector>
#include <sstream>
#include "ObjectLine.hpp"
#include "../LHEF/ParticleLine.hpp"

namespace LHPC
{
  // this is a class to just provide some easy conversion to 4-momenta, with
  // some useful functions.
  class FourMomentum
  {
  public:
    typedef LHEF::ParticleLine const* LhefPointer;
    typedef LHCO::ObjectLine const* LhcoPointer;
    enum VectorComponent
    {
      tComponent = 0,
      xComponent = 1,
      yComponent = 2,
      zComponent = 3
    };

    FourMomentum();
    FourMomentum( double const initialEnergy,
                  double const initialXMomentum,
                  double const initialYMomentum,
                  double const initialZMomentum );
    FourMomentum( FourMomentum const& copySource );
    FourMomentum( LHEF::ParticleLine const& copySource );
    FourMomentum( LHCO::ObjectLine const& copySource );
    FourMomentum( std::pair< LhefPointer, LhefPointer > const& copySource );
    FourMomentum( std::vector< LhefPointer > const& copySource );
    FourMomentum( std::list< LhefPointer > const& copySource );
    FourMomentum( std::pair< LhcoPointer, LhcoPointer > const& copySource );
    FourMomentum( std::vector< LhcoPointer > const& copySource );
    FourMomentum( std::list< LhcoPointer > const& copySource );
    ~FourMomentum();

    void
    assignFrom( FourMomentum const& copySource );
    void
    operator=( FourMomentum const& copySource ){ assignFrom( copySource ); }
    void
    assignFrom( LHEF::ParticleLine const& copySource );
    void
    operator=( LHEF::ParticleLine const& copySource ){ assignFrom(
                                                                copySource ); }
    void
    assignFrom( LHCO::ObjectLine const& copySource );
    void
    operator=( LHCO::ObjectLine const& copySource ){ assignFrom(
                                                                copySource ); }
    double&
    operator[]( int const whichComponent );
    double const&
    operator[]( int const whichComponent ) const;
    void
    operator+=( FourMomentum const& sourceFourMomentum );
    void
    operator-=( FourMomentum const& sourceFourMomentum );
    void
    operator*=( double const scalingFactor );
    void
    operator/=( double const scalingFactor );
    FourMomentum
    operator+( FourMomentum const& sourceFourMomentum );
    FourMomentum
    operator-( FourMomentum const& sourceFourMomentum );
    double
    operator*( FourMomentum const& sourceFourMomentum );
    double
    getT() const{ return (*this)[ (int)tComponent ]; }
    void
    setT( double const inputValue ){ (*this)[ (int)tComponent ] = inputValue; }
    double
    getX() const{ return (*this)[ (int)xComponent ]; }
    void
    setX( double const inputValue ){ (*this)[ (int)xComponent ] = inputValue; }
    double
    getY() const{ return (*this)[ (int)yComponent ]; }
    void
    setY( double const inputValue ){ (*this)[ (int)yComponent ] = inputValue; }
    double
    getZ() const{ return (*this)[ (int)zComponent ]; }
    void
    setZ( double const inputValue ){ (*this)[ (int)zComponent ] = inputValue; }
    double
    getTransverseMagnitudeSquared() const;
    double
    getTransverseMagnitude() const;
    double
    getSpatialMagnitudeSquared() const;
    double
    getSpatialMagnitude() const;
    double
    getInvariantMassSquared() const;
    double
    getInvariantMass() const;
    std::string
    toString() const;


  protected:
    std::vector< double > momentumComponents;
  };





  inline void
  FourMomentum::assignFrom( FourMomentum const& copySource )
  {
    momentumComponents = copySource.momentumComponents;
  }

  inline void
  FourMomentum::assignFrom( LHEF::ParticleLine const& copySource )
  {
    momentumComponents[ (int)tComponent ] = copySource.getEnergy();
    momentumComponents[ (int)xComponent ] = copySource.getXMomentum();
    momentumComponents[ (int)yComponent ] = copySource.getYMomentum();
    momentumComponents[ (int)zComponent ] = copySource.getZMomentum();
  }

  inline void
  FourMomentum::assignFrom( LHCO::ObjectLine const& copySource )
  {
    momentumComponents[ (int)xComponent ]
    = ( copySource.getTransverseMomentum()
        * cos( copySource.getAzimuthalAngle() ) );
    momentumComponents[ (int)yComponent ]
    = ( copySource.getTransverseMomentum()
        * sin( copySource.getAzimuthalAngle() ) );
    momentumComponents[ (int)zComponent ]
    = ( copySource.getTransverseMomentum()
        * sinh( copySource.getPseudorapidity() ) );
    momentumComponents[ (int)tComponent ]
    = sqrt( ( copySource.getInvariantMass() * copySource.getInvariantMass() )
            + getSpatialMagnitudeSquared() );
  }

  inline double&
  FourMomentum::operator[]( int const whichComponent )
  {
    return momentumComponents[ whichComponent ];
  }

  inline double const&
  FourMomentum::operator[]( int const whichComponent ) const
  {
    return momentumComponents[ whichComponent ];
  }

  inline void
  FourMomentum::operator+=( FourMomentum const& sourceFourMomentum )
  {
    momentumComponents[ (int)tComponent ]
    += sourceFourMomentum[ (int)tComponent ];
    momentumComponents[ (int)xComponent ]
    += sourceFourMomentum[ (int)xComponent ];
    momentumComponents[ (int)yComponent ]
    += sourceFourMomentum[ (int)yComponent ];
    momentumComponents[ (int)zComponent ]
    += sourceFourMomentum[ (int)zComponent ];
  }

  inline void
  FourMomentum::operator-=( FourMomentum const& sourceFourMomentum )
  {
    momentumComponents[ (int)tComponent ]
    -= sourceFourMomentum[ (int)tComponent ];
    momentumComponents[ (int)xComponent ]
    -= sourceFourMomentum[ (int)xComponent ];
    momentumComponents[ (int)yComponent ]
    -= sourceFourMomentum[ (int)yComponent ];
    momentumComponents[ (int)zComponent ]
    -= sourceFourMomentum[ (int)zComponent ];
  }

  inline void
  FourMomentum::operator*=( double const scalingFactor )
  {
    momentumComponents[ (int)tComponent ] *= scalingFactor;
    momentumComponents[ (int)xComponent ] *= scalingFactor;
    momentumComponents[ (int)yComponent ] *= scalingFactor;
    momentumComponents[ (int)zComponent ] *= scalingFactor;
  }

  inline void
  FourMomentum::operator/=( double const scalingFactor )
  {
    momentumComponents[ (int)tComponent ] /= scalingFactor;
    momentumComponents[ (int)xComponent ] /= scalingFactor;
    momentumComponents[ (int)yComponent ] /= scalingFactor;
    momentumComponents[ (int)zComponent ] /= scalingFactor;
  }

  inline FourMomentum
  FourMomentum::operator+( FourMomentum const& sourceFourMomentum )
  {
    FourMomentum returnFourMomentum( *this );
    returnFourMomentum += sourceFourMomentum;
    return returnFourMomentum;
  }

  inline FourMomentum
  FourMomentum::operator-( FourMomentum const& sourceFourMomentum )
  {
    FourMomentum returnFourMomentum( *this );
    returnFourMomentum -= sourceFourMomentum;
    return returnFourMomentum;
  }

  inline double
  FourMomentum::operator*( FourMomentum const& sourceFourMomentum )
  {
    return ( ( momentumComponents[ (int)tComponent ]
               * sourceFourMomentum[ (int)tComponent ] )
             - ( momentumComponents[ (int)xComponent ]
                 * sourceFourMomentum[ (int)xComponent ] )
             - ( momentumComponents[ (int)yComponent ]
                 * sourceFourMomentum[ (int)yComponent ] )
             - ( momentumComponents[ (int)zComponent ]
                 * sourceFourMomentum[ (int)zComponent ] ) );
  }

  inline double
  FourMomentum::getTransverseMagnitudeSquared() const
  {
    return ( ( momentumComponents[ (int)xComponent ]
               * momentumComponents[ (int)xComponent ] )
             + ( momentumComponents[ (int)yComponent ]
                 * momentumComponents[ (int)yComponent ] ) );
  }

  inline double
  FourMomentum::getTransverseMagnitude() const
  {
    return sqrt( getTransverseMagnitudeSquared() );
  }

  inline double
  FourMomentum::getSpatialMagnitudeSquared() const
  {
    return ( getTransverseMagnitudeSquared()
             + ( momentumComponents[ (int)zComponent ]
                 * momentumComponents[ (int)zComponent ] ) );
  }

  inline double
  FourMomentum::getSpatialMagnitude() const
  {
    return sqrt( getSpatialMagnitudeSquared() );
  }

  inline double
  FourMomentum::getInvariantMassSquared() const
  {
    return ( ( momentumComponents[ (int)tComponent ]
               * momentumComponents[ (int)tComponent ] )
             - getSpatialMagnitudeSquared() );
  }

  inline double
  FourMomentum::getInvariantMass() const
  {
    return sqrt( getInvariantMassSquared() );
  }

  inline std::string
  FourMomentum::toString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "( " << getT() << ", " << getX() << ", " << getY() << ", " << getZ()
    << " )";
    return stringBuilder.str();
  }

} /* namespace LHPC */
#endif /* FOURMOMENTUM_HPP_ */

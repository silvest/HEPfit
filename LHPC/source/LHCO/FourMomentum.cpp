/*
 * FourMomentum.cpp
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

#include "LHCO.hpp"

namespace LHPC
{
  FourMomentum::FourMomentum() :
      momentumComponents( 4,
                          BOL::UsefulStuff::notANumber )
  {
    // just an initialization list.
  }

  FourMomentum::FourMomentum( double const initialEnergy,
                              double const initialXMomentum,
                              double const initialYMomentum,
                              double const initialZMomentum ) :
      momentumComponents( 4,
                          initialEnergy )
  {
    momentumComponents[ (int)xComponent ] = initialXMomentum;
    momentumComponents[ (int)yComponent ] = initialYMomentum;
    momentumComponents[ (int)zComponent ] = initialZMomentum;
  }

  FourMomentum::FourMomentum( FourMomentum const& copySource ) :
    momentumComponents( copySource.momentumComponents )
  {
    // just an initialization list.
  }

  FourMomentum::FourMomentum( LHEF::ParticleLine const& copySource ) :
      momentumComponents( 4,
                          BOL::UsefulStuff::notANumber )
  {
    assignFrom( copySource );
  }

  FourMomentum::FourMomentum( LHCO::ObjectLine const& copySource ) :
      momentumComponents( 4,
                          BOL::UsefulStuff::notANumber )
  {
    assignFrom( copySource );
  }

  FourMomentum::FourMomentum(
                   std::pair< LhefPointer, LhefPointer > const& copySource ) :
      momentumComponents( 4,
                          0.0 )
  {
    assignFrom( *(copySource.first) );
    momentumComponents[ (int)tComponent ] += copySource.second->getEnergy();
    momentumComponents[ (int)xComponent ] += copySource.second->getXMomentum();
    momentumComponents[ (int)yComponent ] += copySource.second->getYMomentum();
    momentumComponents[ (int)zComponent ] += copySource.second->getZMomentum();
  }

  FourMomentum::FourMomentum( std::vector< LhefPointer > const& copySource ) :
      momentumComponents( 4,
                          0.0 )
  {
    for( int whichElement( copySource.size() - 1 );
         0 <= whichElement;
         --whichElement )
    {
      momentumComponents[ (int)tComponent ]
      += copySource[ whichElement ]->getEnergy();
      momentumComponents[ (int)xComponent ]
      += copySource[ whichElement ]->getXMomentum();
      momentumComponents[ (int)yComponent ]
      += copySource[ whichElement ]->getYMomentum();
      momentumComponents[ (int)zComponent ]
      += copySource[ whichElement ]->getZMomentum();
    }
  }

  FourMomentum::FourMomentum( std::list< LhefPointer > const& copySource ) :
      momentumComponents( 4,
                          0.0 )
  {
    for( std::list< LhefPointer >::const_iterator
         whichLine( copySource.begin() );
         copySource.end() != whichLine;
         ++whichLine )
    {
      momentumComponents[ (int)tComponent ] += (*whichLine)->getEnergy();
      momentumComponents[ (int)xComponent ] += (*whichLine)->getXMomentum();
      momentumComponents[ (int)yComponent ] += (*whichLine)->getYMomentum();
      momentumComponents[ (int)zComponent ] += (*whichLine)->getZMomentum();
    }
  }

  FourMomentum::FourMomentum(
                    std::pair< LhcoPointer, LhcoPointer > const& copySource ) :
      momentumComponents( 4,
                          0.0 )
  {
    assignFrom( *(copySource.first) );
    FourMomentum copyMomentum( *(copySource.second) );
    momentumComponents[ (int)tComponent ] += copyMomentum.getT();
    momentumComponents[ (int)xComponent ] += copyMomentum.getX();
    momentumComponents[ (int)yComponent ] += copyMomentum.getY();
    momentumComponents[ (int)zComponent ] += copyMomentum.getZ();
  }

  FourMomentum::FourMomentum( std::vector< LhcoPointer > const& copySource ) :
      momentumComponents( 4,
                          0.0 )
  {
    FourMomentum copyMomentum;
    for( int whichElement( copySource.size() - 1 );
         0 <= whichElement;
         --whichElement )
    {
      copyMomentum = *(copySource[ whichElement ]);
      momentumComponents[ (int)tComponent ] += copyMomentum.getT();
      momentumComponents[ (int)xComponent ] += copyMomentum.getX();
      momentumComponents[ (int)yComponent ] += copyMomentum.getY();
      momentumComponents[ (int)zComponent ] += copyMomentum.getZ();
    }
  }

  FourMomentum::FourMomentum( std::list< LhcoPointer > const& copySource ) :
      momentumComponents( 4,
                          0.0 )
  {
    FourMomentum copyMomentum;
    for( std::list< LhcoPointer >::const_iterator
         whichLine( copySource.begin() );
         copySource.end() != whichLine;
         ++whichLine )
    {
      copyMomentum = *(*whichLine);
      momentumComponents[ (int)tComponent ] += copyMomentum.getT();
      momentumComponents[ (int)xComponent ] += copyMomentum.getX();
      momentumComponents[ (int)yComponent ] += copyMomentum.getY();
      momentumComponents[ (int)zComponent ] += copyMomentum.getZ();
    }
  }

  FourMomentum::~FourMomentum()
  {
    // does nothing.
  }

} /* namespace LHPC */

/*
 * StandardPreselector.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2013 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "LHEF.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      StandardPreselector::StandardPreselector( LhefParser& lhefParser,
                                                int const soughtParticleCode,
                                            double const transverseMomentumCut,
                                             double const pseudorapidityCut ) :
          filterRules(),
          lineFilter(),
          filteredLines( lineFilter.getFilteredLines() )
      // transverseMomentumCut is a lower bound on the transverse momentum of a
      // ParticleLine to be selected.
      // pseudorapidityCut is an upper bound on the absolute value of the
      // pseudorapidity of a ParticleLine to be selected, and negative a value
      // indicates that no cut on pseudorapidity should be used.
      {
        setUpFilter( lhefParser,
                     std::vector< int >( 1,
                                         soughtParticleCode ),
                     transverseMomentumCut,
                     pseudorapidityCut );
      }

      StandardPreselector::StandardPreselector( LhefParser& lhefParser,
                                 std::vector< int > const& soughtParticleCodes,
                                            double const transverseMomentumCut,
                                             double const pseudorapidityCut ) :
          filterRules(),
          lineFilter(),
          filteredLines( lineFilter.getFilteredLines() )
      // constructor taking a vector of ints to select all particles with any
      // of the given ints as particle code.
      {
        setUpFilter( lhefParser,
                     soughtParticleCodes,
                     transverseMomentumCut,
                     pseudorapidityCut );
      }

      StandardPreselector::StandardPreselector( LhefParser& lhefParser,
                                    MassEigenstate const& soughtMassEigenstate,
                                            double const transverseMomentumCut,
                                             double const pseudorapidityCut ) :
          filterRules(),
          lineFilter(),
          filteredLines( lineFilter.getFilteredLines() )
      // constructor taking a vector of ints to select all particles with any
      // of the given ints as particle code.
      {
        setUpFilter( lhefParser,
                     soughtMassEigenstate.getAllCodes(),
                     transverseMomentumCut,
                     pseudorapidityCut );
      }

      StandardPreselector::~StandardPreselector()
      {
        for( std::vector< FilterRule* >::iterator
             whichRule( filterRules.begin() );
             filterRules.end() > whichRule;
             ++whichRule )
        {
          delete *whichRule;
        }
      }


      void
      StandardPreselector::setUpFilter( LhefParser& lhefParser,
                                 std::vector< int > const& soughtParticleCodes,
                                        double const transverseMomentumCut,
                                        double const pseudorapidityCut )
      {
        filterRules.push_back( new FilterOnState( FilterOnState::finalState,
                                                  true ) );
        filterRules.push_back( new FilterOnParticleCode( soughtParticleCodes,
                                                         true ) );
        if( 0.0 < transverseMomentumCut )
        {
          filterRules.push_back( new FilterOnTransverseMomentum(
                                                         transverseMomentumCut,
                                                                 true ) );
        }
        if( 0.0 < pseudorapidityCut )
        {
          filterRules.push_back( new FilterOnPseudorapidity( pseudorapidityCut,
                                                             true ) );
        }
        for( std::vector< FilterRule* >::iterator
             whichRule( filterRules.begin() );
             filterRules.end() > whichRule;
             ++whichRule )
        {
          lineFilter.addFilterRule( *(*whichRule) );
        }
        lhefParser.registerFilter( lineFilter );
      }

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */

/*
 * StandardPreselector.hpp
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

#ifndef STANDARDPRESELECTOR_HPP_
#define STANDARDPRESELECTOR_HPP_

#include <stdexcept>
#include "../../MEC.hpp"
#include "../FilterRuleClasses.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      // this is a class to streamline the process of creating an
      // AutomaticEventFilter & set of FilterRules to select final-state
      // particles with a given particle code, and optionally, filtering on a
      // minimum transverse momentum and on a maximum pseudorapidity magnitude.
      class StandardPreselector
      {
      public:
        StandardPreselector( LhefParser& lhefParser,
                             int const soughtParticleCode,
                             double const transverseMomentumCut = 0.0,
                             double const pseudorapidityCut = -1.0 );
        // transverseMomentumCut is a lower bound on the transverse momentum of
        // a ParticleLine to be selected.
        // pseudorapidityCut is an upper bound on the absolute value of the
        // pseudorapidity of a ParticleLine to be selected, and negative a
        // value indicates that no cut on pseudorapidity should be used.
        StandardPreselector( LhefParser& lhefParser,
                             std::vector< int > const& soughtParticleCodes,
                             double const transverseMomentumCut = 0.0,
                             double const pseudorapidityCut = -1.0 );
        // constructor taking a vector of ints to select all particles with any
        // of the given ints as particle code.
        StandardPreselector( LhefParser& lhefParser,
                             MassEigenstate const& soughtMassEigenstate,
                             double const transverseMomentumCut = 0.0,
                             double const pseudorapidityCut = -1.0 );
        // constructor taking a MassEigenstate reference instead of an int.
        virtual
        ~StandardPreselector();

        ParticleLine const&
        operator[]( int elementIndex ) const;
        // this accesses the elements starting with index 0 (like normal
        // arrays).
        bool
        isEmpty() const;
        int
        getSize() const;
        std::list< LHPC::LHEF::ParticleLine const* >&
        getList();


      protected:
        std::vector< FilterRule* > filterRules;
        AutomaticEventFilter lineFilter;
        std::list< LHPC::LHEF::ParticleLine const* >& filteredLines;

        void
        setUpFilter( LhefParser& lhefParser,
                     std::vector< int > const& soughtParticleCodes,
                     double const transverseMomentumCut,
                     double const pseudorapidityCut );
      };





      inline ParticleLine const&
      StandardPreselector::operator[]( int elementIndex ) const
      // this accesses the elements starting with index 0 (like normal arrays).
      {
        if( ( 0 > elementIndex )
            ||
            ( filteredLines.size() <= (unsigned int)elementIndex ) )
        {
          throw std::out_of_range(
                        "StandardPreselector::operator[] index out of range" );
        }
        std::list< LHPC::LHEF::ParticleLine const* >::iterator
        lineIterator( filteredLines.begin() );
        while( 0 <= (--elementIndex ) )
        {
          ++lineIterator;
        }
        return *(*lineIterator);
      }

      inline bool
      StandardPreselector::isEmpty() const
      {
        return filteredLines.empty();
      }

      inline int
      StandardPreselector::getSize() const
      {
        return (int)(filteredLines.size());
      }

      inline std::list< LHPC::LHEF::ParticleLine const* >&
      StandardPreselector::getList()
      {
        return filteredLines;
      }

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
#endif /* STANDARDPRESELECTOR_HPP_ */

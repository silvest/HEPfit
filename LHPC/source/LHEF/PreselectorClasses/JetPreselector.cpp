/*
 * JetPreselector.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LHEF.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      JetPreselector::JetPreselector( LhefParser& lhefParser,
                                      bool const includeBottoms,
                                      double const transverseMomentumCut,
                                      double const pseudorapidityCut ) :
          StandardPreselector( lhefParser,
                               jetCodes( includeBottoms ),
                               transverseMomentumCut,
                               pseudorapidityCut )
      {
        // just an initialization list.
      }

      JetPreselector::~JetPreselector()
      {
        // does nothing.
      }


      std::vector< int >
      JetPreselector::jetCodes( bool const includeBottoms )
      {
        std::vector< int > codeVector;
        codeVector.push_back( PDGIX::downOne );
        codeVector.push_back( -PDGIX::downOne );
        codeVector.push_back( PDGIX::downTwo );
        codeVector.push_back( -PDGIX::downTwo );
        codeVector.push_back( PDGIX::upOne );
        codeVector.push_back( -PDGIX::upOne );
        codeVector.push_back( PDGIX::upTwo );
        codeVector.push_back( -PDGIX::upTwo );
        codeVector.push_back( PDGVII::downOne );
        codeVector.push_back( -PDGVII::downOne );
        codeVector.push_back( PDGVII::downTwo );
        codeVector.push_back( -PDGVII::downTwo );
        codeVector.push_back( PDGVII::upOne );
        codeVector.push_back( -PDGVII::upOne );
        codeVector.push_back( PDGVII::upTwo );
        codeVector.push_back( -PDGVII::upTwo );
        if( includeBottoms )
        {
          codeVector.push_back( PDGIX::downThree );
          codeVector.push_back( -PDGIX::downThree );
          codeVector.push_back( PDGVII::downThree );
          codeVector.push_back( -PDGVII::downThree );
        }
        return codeVector;
      }

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */

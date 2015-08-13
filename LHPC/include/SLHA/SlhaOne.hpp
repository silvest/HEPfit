/*
 * SlhaOne.hpp
 *
 *  Created on: Feb 6, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SLHAONE_HPP_
#define SLHAONE_HPP_

#include "BlockTypes.hpp"
#include "SlhaParser.hpp"

namespace LHPC
{
  /* this is just a collection of instances of classes derived from SlhaBlock
   * which covers all the blocks specified in SLHA1 (JHEP 0407 (2004) 036
   * [hep-ph/0311123]).
   */
  class SlhaOne
  {
  public:
    SlhaOne( SlhaParser& fileParser,
             bool const isVerbose = false );
    virtual
    ~SlhaOne();

    SLHA::SparseSinglyIndexedBlock< double > MODSEL;
    SLHA::SparseSinglyIndexedBlock< double > SMINPUTS;
    SLHA::SparseSinglyIndexedBlock< double > MINPAR;
    SLHA::SparseSinglyIndexedBlock< double > EXTPAR;
    SLHA::SparseSinglyIndexedBlock< double > MASS;
    SLHA::DenseDoublyIndexedBlock< double > NMIX;
    SLHA::DenseDoublyIndexedBlock< double > UMIX;
    SLHA::DenseDoublyIndexedBlock< double > VMIX;
    SLHA::DenseDoublyIndexedBlock< double > STOPMIX;
    SLHA::DenseDoublyIndexedBlock< double > SBOTMIX;
    SLHA::DenseDoublyIndexedBlock< double > STAUMIX;
    SLHA::JustSingleValueBlock< double > ALPHA;
    SLHA::SparseSinglyIndexedBlock< double > HMIX;
    SLHA::DenseSinglyIndexedBlock< double > GAUGE;
    SLHA::SparseSinglyIndexedBlock< double > MSOFT;
    SLHA::DenseDoublyIndexedBlock< double > AU;
    SLHA::DenseDoublyIndexedBlock< double > AD;
    SLHA::DenseDoublyIndexedBlock< double > AE;
    SLHA::DenseDoublyIndexedBlock< double > YU;
    SLHA::DenseDoublyIndexedBlock< double > YD;
    SLHA::DenseDoublyIndexedBlock< double > YE;
    SLHA::LinesAsStringsBlock SPINFO;


  protected:
    SlhaParser& fileParser;
    bool const isVerbose;
  };

}

#endif /* SLHAONE_HPP_ */

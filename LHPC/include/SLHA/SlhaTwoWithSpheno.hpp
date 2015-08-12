/*
 * SlhaTwoWithSpheno.hpp
 *
 *  Created on: Feb 22, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SLHATWOWITHSPHENO_HPP_
#define SLHATWOWITHSPHENO_HPP_

#include "SlhaTwo.hpp"

namespace LHPC
{
  /* this is just a collection of instances of classes derived from SlhaBlock
   * which covers all the blocks specified in SLHA2 (Comput. Phys. Commun. 180
   * (2009) 8 [arXiv:0801.0045 [hep-ph]]) & those used by SPheno (Comput. Phys.
   * Commun. 153 (2003) 275 [arXiv:hep-ph/0301101]).
   */
  class SlhaTwoWithSpheno : public SlhaTwo
  {
  public:
    SlhaTwoWithSpheno( SlhaParser& fileParser,
                       bool const isVerbose = false );
    virtual
    ~SlhaTwoWithSpheno();

    SLHA::SparseSinglyIndexedBlock< double > SPHENOINFO;
    SLHA::LinesAsStringsBlock SPHENOINPUT;
    SLHA::LinesAsStringsBlock SPHENOCROSSSECTIONS;
    SLHA::SparseSinglyIndexedBlock< double > SPHENOLOWENERGY;
  };

}

#endif /* SLHATWOWITHSPHENO_HPP_ */

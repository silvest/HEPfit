/*
 * NeutralinosOneToFive.hpp
 *
 *  Created on: Jan 27, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef NEUTRALINOSONETOFIVE_HPP_
#define NEUTRALINOSONETOFIVE_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"
#include "NeutralinosOneToFour.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class NeutralinosOneToFive : public virtual MassSpectrum,
                                 public NeutralinosOneToFour
    {
    public:
      NeutralinosOneToFive( bool const isVerbose = false,
                            std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~NeutralinosOneToFive();

      MassEigenstate&
      getNeutralinoFive();
      MassEigenstate const&
      getNeutralinoFive() const;


    protected:
      MassEigenstate neutralinoFive;
    };



    inline MassEigenstate&
    NeutralinosOneToFive::getNeutralinoFive()
    {
      return neutralinoFive;
    }

    inline MassEigenstate const&
    NeutralinosOneToFive::getNeutralinoFive() const
    {
      return neutralinoFive;
    }

  }

}

#endif /* NEUTRALINOSONETOFIVE_HPP_ */

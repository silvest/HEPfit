/*
 * NeutralinosOneToFour.hpp
 *
 *  Created on: Jan 18, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef NEUTRALINOSONETOFOUR_HPP_
#define NEUTRALINOSONETOFOUR_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class NeutralinosOneToFour : public virtual MassSpectrum
    {
    public:
      NeutralinosOneToFour( bool const isVerbose = false,
                            std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~NeutralinosOneToFour();

      MassEigenstate&
      getNeutralinoOne();
      MassEigenstate const&
      getNeutralinoOne() const;
      MassEigenstate&
      getNeutralinoTwo();
      MassEigenstate const&
      getNeutralinoTwo() const;
      MassEigenstate&
      getNeutralinoThree();
      MassEigenstate const&
      getNeutralinoThree() const;
      MassEigenstate&
      getNeutralinoFour();
      MassEigenstate const&
      getNeutralinoFour() const;
      std::vector< MassEigenstate* >&
      getNeutralinos();
      std::vector< MassEigenstate* > const&
      getNeutralinos() const;


    protected:
      MassEigenstate neutralinoOne;
      MassEigenstate neutralinoTwo;
      MassEigenstate neutralinoThree;
      MassEigenstate neutralinoFour;
      std::vector< MassEigenstate* > neutralinoPointers;
    };



    inline MassEigenstate&
    NeutralinosOneToFour::getNeutralinoOne()
    {
      return neutralinoOne;
    }

    inline MassEigenstate const&
    NeutralinosOneToFour::getNeutralinoOne() const
    {
      return neutralinoOne;
    }

    inline MassEigenstate&
    NeutralinosOneToFour::getNeutralinoTwo()
    {
      return neutralinoTwo;
    }

    inline MassEigenstate const&
    NeutralinosOneToFour::getNeutralinoTwo() const
    {
      return neutralinoTwo;
    }

    inline MassEigenstate&
    NeutralinosOneToFour::getNeutralinoThree()
    {
      return neutralinoThree;
    }

    inline MassEigenstate const&
    NeutralinosOneToFour::getNeutralinoThree() const
    {
      return neutralinoThree;
    }

    inline MassEigenstate&
    NeutralinosOneToFour::getNeutralinoFour()
    {
      return neutralinoFour;
    }

    inline MassEigenstate const&
    NeutralinosOneToFour::getNeutralinoFour() const
    {
      return neutralinoFour;
    }

    inline std::vector< MassEigenstate* >&
    NeutralinosOneToFour::getNeutralinos()
    {
      return neutralinoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    NeutralinosOneToFour::getNeutralinos() const
    {
      return neutralinoPointers;
    }

  }

}

#endif /* NEUTRALINOSONETOFOUR_HPP_ */

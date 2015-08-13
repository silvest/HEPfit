/*
 * CharginosOneToTwo.hpp
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

#ifndef CHARGINOSONETOTWO_HPP_
#define CHARGINOSONETOTWO_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class CharginosOneToTwo : public virtual MassSpectrum
    {
    public:
      CharginosOneToTwo( bool const isVerbose = false,
                         std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~CharginosOneToTwo();

      MassEigenstate&
      getPositiveCharginoOne();
      MassEigenstate const&
      getPositiveCharginoOne() const;
      MassEigenstate&
      getNegativeCharginoOne();
      MassEigenstate const&
      getNegativeCharginoOne() const;
      MassEigenstate&
      getPositiveCharginoTwo();
      MassEigenstate const&
      getPositiveCharginoTwo() const;
      MassEigenstate&
      getNegativeCharginoTwo();
      MassEigenstate const&
      getNegativeCharginoTwo() const;
      std::vector< MassEigenstate* >&
      getPositiveCharginos();
      std::vector< MassEigenstate* > const&
      getPositiveCharginos() const;
      std::vector< MassEigenstate* >&
      getNegativeCharginos();
      std::vector< MassEigenstate* > const&
      getNegativeCharginos() const;


    protected:
      MassEigenstate positiveCharginoOne;
      MassEigenstate negativeCharginoOne;
      MassEigenstate positiveCharginoTwo;
      MassEigenstate negativeCharginoTwo;
      std::vector< MassEigenstate* > positiveCharginoPointers;
      std::vector< MassEigenstate* > negativeCharginoPointers;
    };



    inline MassEigenstate&
    CharginosOneToTwo::getPositiveCharginoOne()
    {
      return positiveCharginoOne;
    }

    inline MassEigenstate const&
    CharginosOneToTwo::getPositiveCharginoOne() const
    {
      return positiveCharginoOne;
    }

    inline MassEigenstate&
    CharginosOneToTwo::getNegativeCharginoOne()
    {
      return negativeCharginoOne;
    }

    inline MassEigenstate const&
    CharginosOneToTwo::getNegativeCharginoOne() const
    {
      return negativeCharginoOne;
    }

    inline MassEigenstate&
    CharginosOneToTwo::getPositiveCharginoTwo()
    {
      return positiveCharginoTwo;
    }

    inline MassEigenstate const&
    CharginosOneToTwo::getPositiveCharginoTwo() const
    {
      return positiveCharginoTwo;
    }

    inline MassEigenstate&
    CharginosOneToTwo::getNegativeCharginoTwo()
    {
      return negativeCharginoTwo;
    }

    inline MassEigenstate const&
    CharginosOneToTwo::getNegativeCharginoTwo() const
    {
      return negativeCharginoTwo;
    }

    inline std::vector< MassEigenstate* >&
    CharginosOneToTwo::getPositiveCharginos()
    {
      return positiveCharginoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    CharginosOneToTwo::getPositiveCharginos() const
    {
      return positiveCharginoPointers;
    }

    inline std::vector< MassEigenstate* >&
    CharginosOneToTwo::getNegativeCharginos()
    {
      return negativeCharginoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    CharginosOneToTwo::getNegativeCharginos() const
    {
      return negativeCharginoPointers;
    }

  }

}

#endif /* CHARGINOSONETOTWO_HPP_ */

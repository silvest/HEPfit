/*
 * SneutrinosOneToThree.hpp
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

#ifndef SNEUTRINOSONETOTHREE_HPP_
#define SNEUTRINOSONETOTHREE_HPP_

#include "../CodesAndDataForMassEigenstates.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class SneutrinosOneToThree : public virtual MassSpectrum
    {
    public:
      SneutrinosOneToThree( bool const isVerbose = false,
                            bool const flavorConserving = false,
                            std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~SneutrinosOneToThree();

      MassEigenstate&
      getAntisneutrinoOne();
      MassEigenstate const&
      getAntisneutrinoOne() const;
      MassEigenstate&
      getElectronAntisneutrinoL(){ return getAntisneutrinoOne(); }
      MassEigenstate const&
      getElectronAntisneutrinoL() const{ return getAntisneutrinoOne(); }
      MassEigenstate&
      getSneutrinoOne();
      MassEigenstate const&
      getSneutrinoOne() const;
      MassEigenstate&
      getElectronSneutrinoL(){ return getSneutrinoOne(); }
      MassEigenstate const&
      getElectronSneutrinoL() const{ return getSneutrinoOne(); }
      MassEigenstate&
      getAntisneutrinoTwo();
      MassEigenstate const&
      getAntisneutrinoTwo() const;
      MassEigenstate&
      getMuonAntisneutrinoL(){ return getAntisneutrinoTwo(); }
      MassEigenstate const&
      getMuonAntisneutrinoL() const{ return getAntisneutrinoTwo(); }
      MassEigenstate&
      getSneutrinoTwo();
      MassEigenstate const&
      getSneutrinoTwo() const;
      MassEigenstate&
      getMuonSneutrinoL(){ return getSneutrinoTwo(); }
      MassEigenstate const&
      getMuonSneutrinoL() const{ return getSneutrinoTwo(); }
      MassEigenstate&
      getAntisneutrinoThree();
      MassEigenstate const&
      getAntisneutrinoThree() const;
      MassEigenstate&
      getTauAntisneutrinoL(){ return getAntisneutrinoThree(); }
      MassEigenstate const&
      getTauAntisneutrinoL() const{ return getAntisneutrinoThree(); }
      MassEigenstate&
      getSneutrinoThree();
      MassEigenstate const&
      getSneutrinoThree() const;
      MassEigenstate&
      getTauSneutrinoL(){ return getSneutrinoThree(); }
      MassEigenstate const&
      getTauSneutrinoL() const{ return getSneutrinoThree(); }
      std::vector< MassEigenstate* >&
      getAntisneutrinos();
      std::vector< MassEigenstate* > const&
      getAntisneutrinos() const;
      std::vector< MassEigenstate* >&
      getSneutrinos();
      std::vector< MassEigenstate* > const&
      getSneutrinos() const;


    protected:
      MassEigenstate antisneutrinoOne;
      MassEigenstate sneutrinoOne;
      MassEigenstate antisneutrinoTwo;
      MassEigenstate sneutrinoTwo;
      MassEigenstate antisneutrinoThree;
      MassEigenstate sneutrinoThree;
      std::vector< MassEigenstate* > antisneutrinoPointers;
      std::vector< MassEigenstate* > sneutrinoPointers;
    };



    inline MassEigenstate&
    SneutrinosOneToThree::getAntisneutrinoOne()
    {
      return antisneutrinoOne;
    }

    inline MassEigenstate const&
    SneutrinosOneToThree::getAntisneutrinoOne() const
    {
      return antisneutrinoOne;
    }

    inline MassEigenstate&
    SneutrinosOneToThree::getSneutrinoOne()
    {
      return sneutrinoOne;
    }

    inline MassEigenstate const&
    SneutrinosOneToThree::getSneutrinoOne() const
    {
      return sneutrinoOne;
    }

    inline MassEigenstate&
    SneutrinosOneToThree::getAntisneutrinoTwo()
    {
      return antisneutrinoTwo;
    }

    inline MassEigenstate const&
    SneutrinosOneToThree::getAntisneutrinoTwo() const
    {
      return antisneutrinoTwo;
    }

    inline MassEigenstate&
    SneutrinosOneToThree::getSneutrinoTwo()
    {
      return sneutrinoTwo;
    }

    inline MassEigenstate const&
    SneutrinosOneToThree::getSneutrinoTwo() const
    {
      return sneutrinoTwo;
    }

    inline MassEigenstate&
    SneutrinosOneToThree::getAntisneutrinoThree()
    {
      return antisneutrinoThree;
    }

    inline MassEigenstate const&
    SneutrinosOneToThree::getAntisneutrinoThree() const
    {
      return antisneutrinoThree;
    }

    inline MassEigenstate&
    SneutrinosOneToThree::getSneutrinoThree()
    {
      return sneutrinoThree;
    }

    inline MassEigenstate const&
    SneutrinosOneToThree::getSneutrinoThree() const
    {
      return sneutrinoThree;
    }
    inline std::vector< MassEigenstate* >&
    SneutrinosOneToThree::getAntisneutrinos()
    {
      return antisneutrinoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SneutrinosOneToThree::getAntisneutrinos() const
    {
      return antisneutrinoPointers;
    }

    inline std::vector< MassEigenstate* >&
    SneutrinosOneToThree::getSneutrinos()
    {
      return sneutrinoPointers;
    }

    inline std::vector< MassEigenstate* > const&
    SneutrinosOneToThree::getSneutrinos() const
    {
      return sneutrinoPointers;
    }

  }

}

#endif /* SNEUTRINOSONETOTHREE_HPP_ */

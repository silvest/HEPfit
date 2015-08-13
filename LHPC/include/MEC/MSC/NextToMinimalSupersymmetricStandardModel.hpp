/*
 * NextToMinimalSupersymmetricStandardModel.hpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef NEXTTOMINIMALSUPERSYMMETRICSTANDARDMODEL_HPP_
#define NEXTTOMINIMALSUPERSYMMETRICSTANDARDMODEL_HPP_

#include "MinimalSupersymmetricStandardModel.hpp"
#include "PSC/NmssmExtraEwsbSpinZeroBosonSet.hpp"
#include "PSC/NeutralinosOneToFive.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    class NextToMinimalSupersymmetricStandardModel :
                                                   public virtual MassSpectrum,
                                                     public MSSM,
                                         public NmssmExtraEwsbSpinZeroBosonSet,
                                                    public NeutralinosOneToFive
    {
    public:
      NextToMinimalSupersymmetricStandardModel( bool const isVerbose = false,
                                       bool const neutrinosAreMajorana = false,
                              std::vector< bool >* const defaultFlags = NULL );
      virtual
      ~NextToMinimalSupersymmetricStandardModel();
    };
    typedef NextToMinimalSupersymmetricStandardModel NMSSM;

  }
  typedef MassSpectrumClass::NMSSM NmssmSpectrum;

}

#endif /* NEXTTOMINIMALSUPERSYMMETRICSTANDARDMODEL_HPP_ */

/*
 * SlhaParser.hpp
 *
 *  Created on: Jan 11, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SLHAPARSER_HPP_
#define SLHAPARSER_HPP_

#include <string>
#include <map>
#include "BOLlib/include/BOLlib.hpp"
#include "../MEC/MassSpectrum.hpp"
#include "../MEC/RunningConstant.hpp"
#include "../MEC/RunningConstantError.hpp"
#include "../MEC/SpectrumUpdater.hpp"
#include "BlockTypes.hpp"

namespace LHPC
{
  // this is a class for reading in a SLHA format file & parsing the data from
  // it.
  class SlhaParser : public BOL::PushingObserved< SpectrumUpdater >
  {
  public:
    static void
    copyWithoutBlock( std::string const& originalFilename,
                      std::string const& blockToStrip,
                      std::string const& copyFilename );

    SlhaParser( bool const shouldRecordBlocksNotRegistered = true,
                bool const isVerbose = true);
    virtual
    ~SlhaParser();

    virtual void
    registerBlock( SLHA::BaseSlhaBlock& blockToUpdate );
    // this registers blockToUpdate so that its data get updated every time a
    // new block of the appropriate name is read.
    void
    registerSpectrum( MassSpectrum& spectrumToUpdate );
    // this adds a pointer to spectrumToUpdate to spectraToUpdate so that its
    // data get updated during each readFile().
    bool
    readFile( std::string const& slhaFileName );
    /* this opens the file with name given by slhaFileName, parses its data
     * into strings, & passes each registered SlhaBlock & each BaseSlhaDecay
     * its appropriate string to interpret.
     */
    SLHA::SameNameBlockSet*
    getBlockAsStrings( std::string blockName );
    SLHA::SameNameBlockSet const*
    getBlockAsStrings( std::string blockName ) const;


  protected:
    bool const isVerbose;
    bool const shouldRecordBlocksNotRegistered;
    BOL::CommentedTextParser fileParser;
    std::map< std::string, SLHA::SameNameBlockSet* > blockMap;
    std::map< std::string, SLHA::SameNameBlockSet* >::iterator
    blockMapIterator;
    std::pair< std::string, SLHA::SameNameBlockSet* > mapInserter;
    SLHA::SameNameBlockSet* currentBlockPointer;
    std::string dataString;
    std::string commentString;
    std::string firstWordOfLine;
    BOL::VectorlikeArray< std::string > wordsOfLine;
    SpectrumUpdater observingSpectrumUpdater;
    bool ownsFmassBlock;
    SLHA::BaseSlhaBlock* fmassBlockPointer;
    std::multimap< int, RunningConstant > const* fmassMap;
    bool ownsFmasserrBlock;
    SLHA::BaseSlhaBlock* fmasserrBlockPointer;
    std::multimap< int, RunningConstantError > const* fmasserrMap;
    bool ownsMassBlock;
    SLHA::BaseSlhaBlock* massBlockPointer;
    std::map< int, double > const* massMap;
    bool successfullyRead;

    void
    clearBlocks();
    // this goes through all the blocks in blockMap & calls their
    // clearEntries() member functions.
    void
    checkForMassBlocksForSpectrum();
    // this ensures that if there is a spectrum to update, there are both
    // blocks for MASS & FMASS.
    void
    prepareToReadNewBlock();
    // this parses the block header line & sets currentBlockPointer
    // appropriately, & calls checkForBlockScaleOrDecayWidth().
    void
    prepareToReadNewDecay();
    // this parses the decay header line & sets currentBlockPointer
    // appropriately, & calls checkForBlockScaleOrDecayWidth().
    void
    recordDecayLine();
    // this interprets the current line as a decay for the spectrum.
    void
    finishUpEitherBlockOrDecay();
    // this sets up the common parsing of the line being read.
    void
    ensureSpectraRecordMasses();
    // this reads the masses from the FMASS & MASS blocks into the spectrum, if
    // necessary.
    void
    prepareForEitherBlockOrDecay();
    // this sets up the common parsing of the line being read.
  };





  inline void
  SlhaParser::registerSpectrum( MassSpectrum& spectrumToUpdate )
  // this adds a pointer to spectrumToUpdate to spectraToUpdate so that its
  // data get updated during each readFile().
  {
    registerObserver( &spectrumToUpdate );
  }

  inline SLHA::SameNameBlockSet*
  SlhaParser::getBlockAsStrings( std::string blockName )
  {
    BOL::StringParser::transformToUppercase( blockName );
    blockMapIterator = blockMap.find( blockName );
    if( blockMap.end() == blockMapIterator )
    {
      return NULL;
    }
    else
    {
      return blockMapIterator->second;
    }
  }

  inline SLHA::SameNameBlockSet const*
  SlhaParser::getBlockAsStrings( std::string blockName ) const
  {
    BOL::StringParser::transformToUppercase( blockName );
    std::map< std::string, SLHA::SameNameBlockSet* >::const_iterator
    constBlockMapIterator( blockMap.find( blockName ) );
    if( blockMap.end() == constBlockMapIterator )
    {
      return NULL;
    }
    else
    {
      return constBlockMapIterator->second;
    }
  }

  inline void
  SlhaParser::clearBlocks()
  // this goes through all the blocks in blockMap & calls their clearEntries()
  // member functions.
  {
    /* observingSpectrumUpdater is set so that only decays will be updated
     * in the observing spectra until later, when the mass map pointers will be
     * set so that the masses can be recorded.
     */
    observingSpectrumUpdater.setFmassMap( NULL );
    observingSpectrumUpdater.setMassMap( NULL );
    updateObservers();
    // the MassSpectrum class over-rides respondToObservedSignal() to call its
    // clearMassesAndDecays() function.
    blockMapIterator = blockMap.begin();
    while( blockMap.end() != blockMapIterator )
    {
      blockMapIterator->second->clearEntries();
      ++blockMapIterator;
    }
    currentBlockPointer = NULL;
  }

  inline void
  SlhaParser::prepareToReadNewDecay()
  // this parses the block header line & sets currentMassEigenstate
  // appropriately, & records its decay width.
  {
    prepareForEitherBlockOrDecay();
    if( !(observerList.empty())
        &&
        ( 3 == wordsOfLine.getSize() ) )
      // if there is a spectrum for recording decays, & if the line has the
      // right number of entries ("DECAY", particle code, decay width)...
    {
      observingSpectrumUpdater.recordDecayHeader( wordsOfLine[ 1 ],
                                                  wordsOfLine[ 2 ] );
    }
  }

  inline void
  SlhaParser::recordDecayLine()
  // this interprets the current line as a decay for the spectrum.
  {
    wordsOfLine.clearEntries();
    BOL::StringParser::parseByChar( dataString,
                                    wordsOfLine,
                                    BOL::StringParser::whitespaceChars );
    if( !(wordsOfLine.isEmpty()) )
      // if there is a decay to record...
    {
      observingSpectrumUpdater.recordDecayLine( wordsOfLine );
    }
  }

  inline void
  SlhaParser::finishUpEitherBlockOrDecay()
  // this either pushes a read decay to the observing spectra or gets the
  // read block to be pushed to its observers.
  {
    if( NULL != currentBlockPointer )
    {
      currentBlockPointer->finishRecordingLines();
      currentBlockPointer = NULL;
    }
    if( observingSpectrumUpdater.isHoldingDecay() )
    {
      updateObservers( observingSpectrumUpdater );
      observingSpectrumUpdater.clearDecay();
    }
  }

  inline void
  SlhaParser::ensureSpectraRecordMasses()
  // this reads the masses from the FMASS & MASS blocks into the spectrum, if
  // necessary.
  {
    if( !(observerList.empty()) )
    {
      observingSpectrumUpdater.setFmassMap( fmassBlockPointer->getFmassMap() );
      observingSpectrumUpdater.setFmasserrMap(
                                      fmasserrBlockPointer->getFmasserrMap() );
      observingSpectrumUpdater.setMassMap( massBlockPointer->getMassMap() );
      updateObservers( observingSpectrumUpdater );
    }
  }

  inline void
  SlhaParser::prepareForEitherBlockOrDecay()
  // this sets up the common parsing of the line being read.
  {
    finishUpEitherBlockOrDecay();
    wordsOfLine.clearEntries();
    BOL::StringParser::parseByChar( dataString,
                                    wordsOfLine,
                                    BOL::StringParser::whitespaceChars );
  }

}

#endif /* SLHAPARSER_HPP_ */

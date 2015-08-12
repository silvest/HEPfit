/*
 * SlhaParser.cpp
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

#include "SLHA.hpp"

namespace LHPC
{
  void
  SlhaParser::copyWithoutBlock( std::string const& originalFilename,
                                std::string const& blockToStrip,
                                std::string const& copyFilename )
  {
    BOL::CommentedTextParser fileParser( "#",
                                         false );
    fileParser.openFile( originalFilename );
    std::string dataString;
    std::string commentString;
    std::string firstWordOfLine;
    std::string restOfLine;
    bool copyingLine( true );
    std::stringstream copyStream;
    while( fileParser.parseNextLineOfFile( dataString,
                                           commentString ) )
    {
      firstWordOfLine.assign( BOL::StringParser::firstWordOf( dataString,
                                                              &restOfLine ) );
      if( BOL::StringParser::stringsMatchIgnoringCase( firstWordOfLine,
                   SLHA::BlockClass::BaseStringBlock::decayIdentifierString ) )
      {
        // if we find a DECAY, then we're not stripping out a block any more.
        copyingLine = true;
      }
      else if( BOL::StringParser::stringsMatchIgnoringCase( firstWordOfLine,
                   SLHA::BlockClass::BaseStringBlock::blockIdentifierString ) )
      {
        /* if we find a BLOCK, then we're not stripping any more if this BLOCK
         * is not the type that we are stripping, but we are (still) stripping
         * if it is.
         */
        copyingLine = !( BOL::StringParser::stringsMatchIgnoringCase(
                                  BOL::StringParser::firstWordOf( restOfLine ),
                                                              blockToStrip ) );
      }
      if( copyingLine
          &&
          !( dataString.empty()
             && commentString.empty() ) )
      {
        // if we're not stripping out a block, record the line again (except
        // for blank lines).
        copyStream << dataString << commentString << "\n";
      }
    }
    std::ofstream outputStream( copyFilename.c_str() );
    outputStream << copyStream.str();
    outputStream.close();
  }

  SlhaParser::SlhaParser( bool const shouldRecordBlocksNotRegistered,
                          bool const isVerbose ) :
      isVerbose( isVerbose ),
      shouldRecordBlocksNotRegistered( shouldRecordBlocksNotRegistered ),
      fileParser( "#",
                  this->isVerbose ),
      blockMap(),
      blockMapIterator(),
      mapInserter( "",
                   NULL ),
      currentBlockPointer( NULL ),
      dataString( "" ),
      commentString( "" ),
      firstWordOfLine( "" ),
      wordsOfLine( 2 ),
      observingSpectrumUpdater(),
      ownsFmassBlock( false ),
      fmassBlockPointer( NULL ),
      fmassMap(),
      ownsFmasserrBlock( false ),
      fmasserrBlockPointer( NULL ),
      fmasserrMap(),
      ownsMassBlock( false ),
      massBlockPointer( NULL ),
      massMap(),
      successfullyRead( false )
  {
    // just an initialization list.
  }

  SlhaParser::~SlhaParser()
  {
    if( ownsMassBlock )
    {
      delete massBlockPointer;
    }
    if( ownsFmassBlock )
    {
      delete fmassBlockPointer;
    }
    blockMapIterator = blockMap.begin();
    while( blockMap.end() != blockMapIterator )
    {
      delete blockMapIterator->second;
      ++blockMapIterator;
    }
  }



  void
  SlhaParser::registerBlock( SLHA::BaseSlhaBlock& blockToUpdate )
  // this registers blockToUpdate so that its data get updated every time a
  // new block of the appropriate name is read.
  {
    blockMapIterator = blockMap.find( blockToUpdate.getName() );
    if( blockMap.end() == blockMapIterator )
    {
      mapInserter.first.assign( blockToUpdate.getName() );
      currentBlockPointer = new SLHA::SameNameBlockSet( mapInserter.first );
      mapInserter.second = currentBlockPointer;
      blockMap.insert( mapInserter );
    }
    else
    {
      currentBlockPointer = blockMapIterator->second;
    }
    if( blockToUpdate.isFmassBlock() )
    {
      if( ownsFmassBlock )
      {
        delete fmassBlockPointer;
        ownsFmassBlock = false;
      }
      fmassBlockPointer = &blockToUpdate;
    }
    else if( blockToUpdate.isFmasserrBlock() )
    {
      if( ownsFmasserrBlock )
      {
        delete fmasserrBlockPointer;
        ownsFmasserrBlock = false;
      }
      fmasserrBlockPointer = &blockToUpdate;
    }
    else if( blockToUpdate.isMassBlock() )
    {
      if( ownsMassBlock )
      {
        delete massBlockPointer;
        ownsMassBlock = false;
      }
      massBlockPointer = &blockToUpdate;
    }
    currentBlockPointer->registerObserver( &blockToUpdate );
  }

  bool
  SlhaParser::readFile( std::string const& slhaFileName )
  /* this opens the file with name given by slhaFileName, parses its data
   * into strings, stores them in its BaseStringBlocks within its
   * SameNameBlockSets, & updates all the observing SlhaBlocks & MassSpectrums.
   */
  {
    checkForMassBlocksForSpectrum();
    clearBlocks();
    // new data should overwrite the old data, not append to it.
    successfullyRead = fileParser.openFile( slhaFileName );
    while( fileParser.parseNextLineOfFile( dataString,
                                           commentString ) )
    {
      firstWordOfLine.assign( BOL::StringParser::firstWordOf( dataString ) );
      // some perverts use tabs in their SLHA files.
      BOL::StringParser::transformToLowercase( firstWordOfLine );
      if( BOL::StringParser::stringsMatchIgnoringCase( firstWordOfLine,
                   SLHA::BlockClass::BaseStringBlock::blockIdentifierString ) )
      {
        prepareToReadNewBlock();
      }
      else if( BOL::StringParser::stringsMatchIgnoringCase( firstWordOfLine,
                   SLHA::BlockClass::BaseStringBlock::decayIdentifierString ) )
      {
        prepareToReadNewDecay();
      }
      else if( NULL != currentBlockPointer )
        // otherwise, if a block is being recorded...
      {
        currentBlockPointer->recordBodyLine( dataString,
                                             commentString );
      }
      else if( observingSpectrumUpdater.isHoldingDecay() )
        // otherwise, if a decay is being recorded...
      {
        recordDecayLine();
      }
      commentString.clear();
    }
    // the last block is closed after the last line has been read in:
    finishUpEitherBlockOrDecay();
    // after reading in the file, any recorded masses are passed to
    // spectrumToUpdate if it is not NULL:
    ensureSpectraRecordMasses();
    return successfullyRead;
  }

  void
  SlhaParser::checkForMassBlocksForSpectrum()
  // this ensures that if there is a spectrum to update, there are both
  // blocks for MASS & FMASS.
  {
    if( !(observerList.empty()) )
    {
      if( NULL == fmassBlockPointer )
        // if there is at least 1 spectrum to update, but no fmass block...
      {
        ownsFmassBlock = true;
        fmassBlockPointer
        = new SLHA::SinglyIndexedMultipleEntriesBlock< RunningConstant >(
                                                                       "FMASS",
                                                             RunningConstant(),
                                                                     isVerbose,
                                                                          9 );
        blockMapIterator = blockMap.find( fmassBlockPointer->getName() );
        if( blockMap.end() == blockMapIterator )
        {
          mapInserter.first.assign( fmassBlockPointer->getName() );
          currentBlockPointer
          = new SLHA::SameNameBlockSet( mapInserter.first );
          mapInserter.second = currentBlockPointer;
          blockMap.insert( mapInserter );
        }
        else
        {
          currentBlockPointer = blockMapIterator->second;
        }
        currentBlockPointer->registerObserver( fmassBlockPointer );
      }
      if( NULL == fmasserrBlockPointer )
        // if there is at least 1 spectrum to update, but no fmasserr block...
      {
        ownsFmasserrBlock = true;
        fmasserrBlockPointer
        = new SLHA::SinglyIndexedMultipleEntriesBlock< RunningConstantError >(
                                                                    "FMASSERR",
                                                        RunningConstantError(),
                                                                     isVerbose,
                                                                           9 );
        blockMapIterator = blockMap.find( fmasserrBlockPointer->getName() );
        if( blockMap.end() == blockMapIterator )
        {
          mapInserter.first.assign( fmasserrBlockPointer->getName() );
          currentBlockPointer
          = new SLHA::SameNameBlockSet( mapInserter.first );
          mapInserter.second = currentBlockPointer;
          blockMap.insert( mapInserter );
        }
        else
        {
          currentBlockPointer = blockMapIterator->second;
        }
        currentBlockPointer->registerObserver( fmasserrBlockPointer );
      }
      if( NULL == massBlockPointer )
        // if there is at least 1 spectrum to update, but no mass block...
      {
        ownsMassBlock = true;
        massBlockPointer
        = new SLHA::SparseSinglyIndexedBlock< double >( "MASS",
                                                  BOL::UsefulStuff::notANumber,
                                                        isVerbose,
                                                        9 );
        blockMapIterator = blockMap.find( massBlockPointer->getName() );
        if( blockMap.end() == blockMapIterator )
        {
          mapInserter.first.assign( massBlockPointer->getName() );
          currentBlockPointer
          = new SLHA::SameNameBlockSet( mapInserter.first );
          mapInserter.second = currentBlockPointer;
          blockMap.insert( mapInserter );
        }
        else
        {
          currentBlockPointer = blockMapIterator->second;
        }
        currentBlockPointer->registerObserver( massBlockPointer );
      }
    }
  }

  void
  SlhaParser::prepareToReadNewBlock()
  // this parses the block header line & sets currentBlockPointer
  // appropriately, & calls checkForBlockScaleOrDecayWidth().
  {
    prepareForEitherBlockOrDecay();
    if( 2 <= wordsOfLine.getSize() )
      // if the line at least has a string for the block name...
    {
      BOL::StringParser::transformToUppercase( wordsOfLine[ 1 ] );
      blockMapIterator = blockMap.find( wordsOfLine[ 1 ] );
      if( blockMap.end() != blockMapIterator )
        // if the name corresponds to a block that already exists in the map...
      {
        currentBlockPointer = blockMapIterator->second;
      }
      else if( shouldRecordBlocksNotRegistered )
        // otherwise, if it should be recorded anyway...
      {
        mapInserter.first.assign( wordsOfLine[ 1 ] );
        currentBlockPointer = new SLHA::SameNameBlockSet( mapInserter.first );
        // a new block is required.
        mapInserter.second = currentBlockPointer;
        blockMap.insert( mapInserter );
      }
      if( NULL != currentBlockPointer )
      {
        double currentBlockScale( 0.0 );
        // if no scale is given by "Q=", the block is assigned the default
        // scale of 0.0 GeV.
        if( 3 <= wordsOfLine.getSize() )
        {
          currentBlockScale =
          BOL::StringParser::stringToDouble( wordsOfLine.getBack() );
        }
        currentBlockPointer->recordHeader( dataString,
                                           commentString,
                                           currentBlockScale );
      }
    }
  }

}

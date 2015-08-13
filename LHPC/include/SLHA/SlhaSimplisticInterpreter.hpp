/*
 * SlhaSimplisticInterpreter.hpp
 *
 *  Created on: Sep 13, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHASIMPLISTICINTERPRETER_HPP_
#define SLHASIMPLISTICINTERPRETER_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include "SlhaParser.hpp"
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  /* this is a class for an object that reads in an SLHA file & then returns
   * strings interpreting given keys as block names with sets of index
   * integers. it also includes functionality to interpret the string as a
   * double or an int.
   */
  class SlhaSimplisticInterpreter
  {
  public:
    SlhaSimplisticInterpreter( std::string const& slhaFilename );
    ~SlhaSimplisticInterpreter();

    std::string
    operator()( std::string blockNameAndIndices );
    std::string
    withMap( std::string blockNameAndIndices );
    double
    getDouble( std::string blockNameAndIndices );
    int
    getInt( std::string blockNameAndIndices );
    double
    getLowestScale( std::string const& blockName ) const;
    bool
    readFile( std::string const& slhaFileName );
    // this opens the file with name given by slhaFileName with slhaParser.


  protected:
    SlhaParser slhaParser;
    std::stringstream stringParser;
    std::map< std::string, std::string > keyedResults;
    std::map< std::string, std::string >::iterator mapIterator;

    std::stringstream&
    getStringParser( std::string const& newStringForParser );
  };




  inline std::string
  SlhaSimplisticInterpreter::withMap( std::string blockNameAndIndices )
  {
    mapIterator = keyedResults.find( blockNameAndIndices );
    if( keyedResults.end() == mapIterator )
    {
      mapIterator = keyedResults.insert( keyedResults.begin(),
                                         std::pair< std::string, std::string >(
                                                           blockNameAndIndices,
                                            (*this)( blockNameAndIndices ) ) );
    }
    return mapIterator->second;
  }

  inline double
  SlhaSimplisticInterpreter::getDouble( std::string blockNameAndIndices )
  {
    double returnValue( 0.0 );
    getStringParser( (*this)( blockNameAndIndices ) ) >> returnValue;
    return returnValue;
  }

  inline int
  SlhaSimplisticInterpreter::getInt( std::string blockNameAndIndices )
  {
    int returnValue( 0 );
    getStringParser( (*this)( blockNameAndIndices ) ) >> returnValue;
    return returnValue;
  }

  inline double
  SlhaSimplisticInterpreter::getLowestScale(
                                           std::string const& blockName ) const
  {
    return (*slhaParser.getBlockAsStrings( blockName ))[ 0 ].getScale();
  }

  inline bool
  SlhaSimplisticInterpreter::readFile( std::string const& slhaFileName )
  // this opens the file with name given by slhaFileName with slhaParser.
  {
    keyedResults.clear();
    return slhaParser.readFile( slhaFileName );
  }

  inline std::stringstream&
  SlhaSimplisticInterpreter::getStringParser(
                                        std::string const& newStringForParser )
  {
    stringParser.clear();
    stringParser.str( newStringForParser );
    return stringParser;
  }

} /* namespace LHPC */
#endif /* SLHASIMPLISTICINTERPRETER_HPP_ */

/*
 * AutomaticEventFilter.hpp
 *
 *  Created on: Jan 26, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef AUTOMATICEVENTFILTER_HPP_
#define AUTOMATICEVENTFILTER_HPP_

#include <list>
#include "LhefEvent.hpp"
#include "ParticleLine.hpp"
#include "FilterRule.hpp"

namespace LHPC
{
  namespace LHEF
  {
    // this class is for automatically filtering ParticleLines from a LhefEvent
    // according to specified selections or vetoes.
    class AutomaticEventFilter
    {
    public:
      AutomaticEventFilter();
      ~AutomaticEventFilter();

      AutomaticEventFilter&
      addFilterRule( FilterRule const& ruleToAdd );
      // the filters get applied to the events in the order in which they are
      // added with this function.
      std::list< ParticleLine const* >&
      getFilteredLines();

      // functions for automatic updating:
      void
      updateForNewEvent( LhefEvent const& eventSource );
      /* this clears filteredLines then goes through eventSource & records in
       * filteredLines all pointers to ParticleLines that pass the filters.
       * this is automatically called by the LhefParser that created this
       * instance when it reads in a new LHEF event, so it shouldn't have to be
       * called by anything else.
       */

    protected:
      std::list< ParticleLine const* > filteredLines;
      std::vector< FilterRule const* > filterRules;
      int numberOfRules;
      int rulesToPass;

      bool
      passesAllFilters( ParticleLine const& lineToCheck );
    };





    inline AutomaticEventFilter&
    AutomaticEventFilter::addFilterRule( FilterRule const& ruleToAdd )
    // the filters get applied to the events in the order in which they are
    // added with this function.
    {
      filterRules.push_back( &ruleToAdd );
      ++numberOfRules;
      return *this;
    }

    inline std::list< ParticleLine const* >&
    AutomaticEventFilter::getFilteredLines()
    {
      return filteredLines;
    }

    inline void
    AutomaticEventFilter::updateForNewEvent( LhefEvent const& eventSource )
    /* this clears filteredLines then goes through eventSource & records in
     * filteredLines all pointers to ParticleLines that pass the filters.
     * this is automatically called by the LhefParser that created this
     * instance when it reads in a new LHEF event, so it shouldn't have to be
     * called by anything else.
     */
    {
      filteredLines.clear();
      for( int whichLine( eventSource.getNumberOfParticles() );
           0 < whichLine;
           --whichLine )
      {
        if( passesAllFilters( eventSource[ whichLine ] ) )
        {
          filteredLines.push_back( &(eventSource[ whichLine ]) );
        }
      }
    }

    inline bool
    AutomaticEventFilter::passesAllFilters( ParticleLine const& lineToCheck )
    /* this clears filteredLines then goes through eventSource & records in
     * filteredLines all pointers to ParticleLines that pass the filters.
     * this is automatically called by the LhefParser that created this
     * instance when it reads in a new LHEF event, so it shouldn't have to be
     * called by anything else.
     */
    {
      rulesToPass = numberOfRules;
      while( ( 0 < rulesToPass )
             &&
             ( (*(filterRules[ numberOfRules
                               - rulesToPass ]))( lineToCheck ) ) )
      {
        --rulesToPass;
      }
      return( 0 == rulesToPass );
    }

  }

}

#endif /* AUTOMATICEVENTFILTER_HPP_ */

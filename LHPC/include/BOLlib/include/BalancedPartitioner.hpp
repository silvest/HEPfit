/*
 * BalancedPartitioner.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef BALANCEDPARTITIONER_HPP_
#define BALANCEDPARTITIONER_HPP_

#include <vector>
#include <list>
#include "VectorlikeArray.hpp"
#include "BalancedPartitionCandidate.hpp"

namespace BOL
{
  // this template class takes a std::vector of pointers to instances of
  // ClassToBeSorted & a pointer to a function which takes a const reference to
  // a ClassToBeSorted instance & returns a double, & partitions the
  // std::vector according to the doubles returned from the provided
  // ClassToBeSorted references.
  // basically this takes a set of objects that have some way of providing
  // weights as doubles, & finds the set of objects such that the sum of their
  // their weights is as close to half of the total weight as possible, thus in
  // effect providing 2 sets with as equal weights as possible.
  // the output of this class is a vector of unsigned ints which are the
  // indices of the elements in the original vector which belong to either the
  // collection with lower or higher weight, depending on whether
  // getLowerWeightPartitionIndices() or getHigherWeightPartitionIndices() is
  // called.
  template< typename ClassToBeSorted >
  class BalancedPartitioner
  {
  public:
    BalancedPartitioner( double (*weightGetter)( ClassToBeSorted const& ) );
    ~BalancedPartitioner();

    void
    makeMinimumDifferencePartition(
                      std::vector< ClassToBeSorted > const& vectorToBeSorted );

    std::vector< unsigned int > const&
    getLowerWeightPartitionIndices() const;
    // this returns a reference to an internally-held std::vector of
    // unsigned ints corresponding to indices in the vector given to the last
    // call of makeMinimumDifferencePartition(...). the collection of objects
    // from that vector with these indices is the partition with GREATEST
    // weight UNDER (or equal to) half the the total weight of the set of
    // objects last given to makeMinimumDifferencePartition(...).
    // NOTE BENE: this std::vector gets OVER-WRITTEN by subsequent calls to
    // makeMinimumDifferencePartition(...)! if you want to keep the returned
    // set, you have to make your own *copy* of the std::vector referenced by
    // the return value of this function.
    double
    getLowerWeight() const;
    // this returns the sum of absolute weights of the objects with indices in
    // the return reference of getLowerWeightPartitionIndices().
    std::vector< unsigned int > const&
    getHigherWeightPartitionIndices() const;
    // this returns a reference to an internally-held std::vector of
    // unsigned ints corresponding to indices in the vector given to the last
    // call of makeMinimumDifferencePartition(...). the collection of objects
    // from that vector with these indices is the partition with LEAST weight
    // OVER (or equal to) half the the total weight of the set of objects last
    // given to makeMinimumDifferencePartition(...).
    // NOTE BENE: this std::vector gets OVER-WRITTEN by subsequent calls to
    // makeMinimumDifferencePartition(...)! if you want to keep the returned
    // set, you have to make your own *copy* of the std::vector referenced by
    // the return value of this function.
    double
    getHigherWeight() const;
    // this returns the sum of absolute weights of the objects with indices in
    // the return reference of getHigherWeightPartitionIndices().


  protected:
    typedef typename std::pair< unsigned int, double > IndexWithWeight;
    typedef typename
    std::list< IndexWithWeight >::iterator IndexWithWeightListIterator;
    typedef typename
    std::list< IndexWithWeight >::const_iterator
    IndexWithWeightListConstIterator;
    typedef typename
    std::vector< ClassToBeSorted >::iterator ClassToBeSortedIterator;
    typedef typename
    std::list< BalancedPartitionCandidate* >::iterator CandidateListIterator;
    typedef typename
    std::vector< unsigned int >::const_iterator
    CandidateVectorConstIterator;

    unsigned int sizeOfVectorToSort;
    // these are the sets of objects with the weight closest to half the total
    // weight of the given argument set:
    std::vector< unsigned int > lowerWeightPartition;
    double lowerWeight;
    std::vector< unsigned int > higherWeightPartition;
    double higherWeight;
    double (*weightGetter)( ClassToBeSorted const& );
    double weightSum;
    double positiveWeight;
    std::list< IndexWithWeight > orderedList;
    // this holds indices of the objects to be partitioned, ordered by weight.
    BalancedPartitionCandidate* currentCandidate;
    BalancedPartitionCandidate* bestCandidateSoFar;
    double currentWeight;
    VectorlikeArray< BalancedPartitionCandidate >
    vectorManagingCandidateMemory;
    // this is a vector of candidates for the partitionedVector. it is for
    // memory management: the actual set of candidates used for finding the
    // partition is:
    std::list< BalancedPartitionCandidate* > candidatesToBuildOn;
    std::list< BalancedPartitionCandidate* > extraCandidates;
    // this is used to add the candidates from building with a given object to
    // candidatesToBuildOn after the object has finished looking through
    // candidatesToBuildOn, to avoid double-counting.

    void
    clearSavedCandidates();

    static bool
    orderByDouble( IndexWithWeight const& firstPair,
                   IndexWithWeight const& secondPair );
    // this returns true if firstPair->second <= secondPair->second, false
    // otherwise.

    void
    fillVectorFromBestCandidate(
                              std::vector< unsigned int >& partitionToBeMade );
    void
    makeOtherPartition( std::vector< unsigned int >& partitionToBeMade,
                        std::vector< unsigned int >& partitionAlreadyMade );
    // this REMOVES ALL OF THE PAIRS WITH THE SAME INDICES AS IN
    // partitionAlreadyMade FROM orderedList & then fills partitionToBeMade
    // with the indices that remain.
  };




  template< typename ClassToBeSorted >
  inline
  BalancedPartitioner< ClassToBeSorted >::BalancedPartitioner(
                           double (*weightGetter)( ClassToBeSorted const& ) ) :
      sizeOfVectorToSort( 0 ),
      lowerWeightPartition(),
      lowerWeight( 0.0 ),
      higherWeightPartition(),
      higherWeight( 0.0 ),
      weightGetter( weightGetter ),
      weightSum( 0.0 ),
      positiveWeight( 0.0 ),
      orderedList(),
      currentCandidate( NULL ),
      bestCandidateSoFar( NULL ),
      currentWeight( 0.0 ),
      vectorManagingCandidateMemory(),
      candidatesToBuildOn(),
      extraCandidates()
  {
    // just an initialization list.
  }

  template< typename ClassToBeSorted >
  inline
  BalancedPartitioner< ClassToBeSorted >::~BalancedPartitioner()
  {
    // does nothing.
  }

  template< typename ClassToBeSorted >
  inline void
  BalancedPartitioner< ClassToBeSorted >::makeMinimumDifferencePartition(
                       std::vector< ClassToBeSorted > const& vectorToBeSorted )
  {
    // 1st, we remove all the stuff that might have been hanging around since
    // the last partitioning:
    clearSavedCandidates();
    weightSum = 0.0;

    sizeOfVectorToSort = vectorToBeSorted.size();
    // next we make a list of the objects with the weights (absolute values, &
    // we also ignore any objects with zero weight) & find the total weight:
    for( unsigned int whichIndex( 0 );
         sizeOfVectorToSort > whichIndex;
         ++whichIndex )
    {
      positiveWeight = (*weightGetter)( vectorToBeSorted[ whichIndex ] );
      if( 0.0 > positiveWeight )
      {
        positiveWeight = -positiveWeight;
      }
      if( 0.0 < positiveWeight )
      {
        weightSum += positiveWeight;
        orderedList.push_back( IndexWithWeight( whichIndex,
                                                positiveWeight ) );
      }
    }

    double const halfTotalWeight( 0.5 * weightSum );
    if( orderedList.empty() )
      // if the given set didn't have enough objects to sort...
    {
      // the vectors that get returned have already been cleared.
      lowerWeight = 0.0;
      higherWeight = 0.0;
    }
    else
    {
      orderedList.sort( &orderByDouble );
      // now orderedList is all the objects with non-zero weight, ordered by
      // their absolute weights.

      if( halfTotalWeight <= orderedList.back().second )
      // if the object with greatest weight takes up half or more of the
      // total weight by itself...
      {
        // it alone becomes higherWeightPartition.
        higherWeightPartition.push_back( orderedList.back().first );
        higherWeight = orderedList.back().second;
        lowerWeight = ( weightSum - higherWeight );
        makeOtherPartition( lowerWeightPartition,
                            higherWeightPartition );
      }
      else if( orderedList.size() < 4 )
        // or if the given set was too small to have any non-trivial
        // partitions...
      {
        // the object with greatest weight alone becomes lowerWeightPartition
        // (if it should become higherWeightPartition, it'll have happened in
        // the previous if statement):
        lowerWeightPartition.push_back( orderedList.back().first );
        lowerWeight = orderedList.back().second;
        higherWeight = ( weightSum - lowerWeight );
        makeOtherPartition( higherWeightPartition,
                            lowerWeightPartition );
      }
      else
      // otherwise, we need to check all the required combinations (partitions
      // of sizes up to half the number of objects, with weights less than half
      // the total weight, *plus* cases where the weight is over half the total
      // weight, if removing a single element would bring it back under half
      // the total weight):
      {
        // 1st we note that the highest-weight element on its own is still the
        // leading candidate partition:
        bestCandidateSoFar = &(vectorManagingCandidateMemory.newEnd());
        bestCandidateSoFar->getVector().assign( 1,
                                                orderedList.back().first );
        lowerWeight = orderedList.back().second;
        higherWeight = ( weightSum - lowerWeight );
        bestCandidateSoFar->setWeight( lowerWeight );
        bestCandidateSoFar->setInverseWeight( higherWeight );
        // now bestCandidate so far is the highest-weight element on its own, &
        // vectorManagingCandidateMemory is holding a pointer to it in case the
        // bestCandidateSoFar pointer gets set to pointing at a different
        // partitionCandidate.


        // now we start going through candidate partitions which consist of 2
        // or more elements.
        IndexWithWeightListConstIterator listIterator( orderedList.begin() );
        ++listIterator;
        // listIterator is now at the pair with the 2nd object by weight.
        currentCandidate = &(vectorManagingCandidateMemory.newEnd());
        // note where the new partitionCandidate is in memory, so that we can
        // delete it later.
        candidatesToBuildOn.push_back( currentCandidate );
        currentCandidate->getVector().assign( 2,
                                              orderedList.front().first );
        // currentCandidate is now the set { [lowest weight element] }.
        currentCandidate->getVector()[ 1 ] = listIterator->first;
        // currentCandidate is now the set
        // { [lowest weight element], [2nd-lowest weight element] }.
        currentWeight
        = ( orderedList.front().second + listIterator->second );
        // the only way that this combination could have half the weight or
        // over is if there was only 3 elements or less, which should have
        // been caught by earlier if statements.
        currentCandidate->setWeight( currentWeight );
        currentCandidate->setInverseWeight( ( weightSum - currentWeight ) );

        if( lowerWeight < currentWeight )
        // if this is a better candidate than the single highest-weight
        // object...
        {
          // note the weights of this new current best candidate:
          lowerWeight = currentWeight;
          higherWeight = ( weightSum - currentWeight );
          // note the actual new best current candidate:
          bestCandidateSoFar = currentCandidate;
        }

        IndexWithWeightListConstIterator lowerWeightIterator;
        // this is used for making new candidates of 2 elements.

        // now we loop over the rest of the objects in orderedList (all but the
        // 2 lowest-weight objects) to build all partitions with 2 or more
        // objects.
        // we don't bother building with partition candidates which are already
        // over half the total weight; however, we do keep those which are over
        // half the weight after being built from a candidate with less than
        // half the weight, if its weight is less that the higherWeight of the
        // best candidate so far, in which case it becomes the new best
        // candidate.
        // we also don't bother building candidates with more than half the
        // total number of elements, because their partners should already have
        // been found.

        // move on to the 3rd object by weight to begin the loop:
        ++listIterator;
        while( orderedList.end() != listIterator )
        // until we reach the end of orderedList...
        {
          for( CandidateListIterator
               candidateIterator( candidatesToBuildOn.begin() );
               candidatesToBuildOn.end() != candidateIterator;
               ++candidateIterator )
          {
            if( ( 2 * ( (*candidateIterator)->getVector().size() + 1 ) )
               <= orderedList.size() )
              // if we're not building a candidate that consists of too many
              // objects...
            {
              currentWeight = ( (*candidateIterator)->getWeight()
                               + listIterator->second );

              if( halfTotalWeight < currentWeight )
                // if this candidate is already over half the weight, it cannot
                // be used to build upon, so we only keep it if it's a better
                // candidate than bestCandidateSoFar. in either case, we remove
                // *candidateIterator from candidatesToBuildOn because it's
                // not going to be a better candidate when paired with an
                // object of even higher weight.
              {
                if( higherWeight > currentWeight )
                  // if this *is* a better candidate (including not consisting
                  // of too many elements; if a candidate has half the elements
                  // but over half the weight, the its partner has half the
                  // elements & less than half the weight, so it'll be found
                  // anyway)...
                {
                  // note the weights of this new current best candidate:
                  higherWeight = currentWeight;
                  lowerWeight = ( weightSum - currentWeight );
                  // make the new best candidate & manage the memory:
                  bestCandidateSoFar
                  = &(vectorManagingCandidateMemory.newEnd());
                  // set the candidate's properties:
                  bestCandidateSoFar->buildFrom( *(*candidateIterator),
                                                 listIterator->first,
                                                 higherWeight,
                                                 lowerWeight );
                }
                // otherwise, we don't do anything with this combination of
                // *candidateIterator & *listIterator.
                // in either case we remove *candidateIterator from
                // candidatesToBuildOn because it's not going to be a better
                // candidate when paired with an object of even higher weight.
                candidateIterator
                = candidatesToBuildOn.erase( candidateIterator );
              }
              else
                // otherwise, currentWeight <= halfTotalWeight, so we check to
                // see if it's a better candidate than bestCandidateSoFar, &
                // also if we should keep it for building on, if it does not
                // consist of too many objects already...
              {
                if( lowerWeight < currentWeight )
                  // if this *is* a better candidate...
                {
                  // note the weights of this new current best candidate:
                  lowerWeight = currentWeight;
                  higherWeight = ( weightSum - currentWeight );
                  // make the new best candidate & manage the memory:
                  bestCandidateSoFar
                  = &(vectorManagingCandidateMemory.newEnd());
                  // set the candidate's properties:
                  bestCandidateSoFar->buildFrom( *(*candidateIterator),
                                                 listIterator->first,
                                                 lowerWeight,
                                                 higherWeight );

                  if( ( 2 * ( (*candidateIterator)->getVector().size() + 2 ) )
                      <= orderedList.size() )
                    // if building on this would make a candidate that still
                    // would not consist of too many objects...
                  {
                    extraCandidates.push_back( bestCandidateSoFar );
                  }
                }
                else if( ( 2 * ( (*candidateIterator)->getVector().size()
                                 + 2 ) )
                         <= orderedList.size() )
                 // otherwise if we're building a candidate which is not the
                 // best candidate so far, but building on it would make a
                 // candidate that still would not consist of too many
                 // objects...
                {
                  // make the new candidate & manage the memory:
                  currentCandidate
                  = &(vectorManagingCandidateMemory.newEnd());
                  extraCandidates.push_back( currentCandidate );
                  // set the candidate's properties:
                  currentCandidate->buildFrom( *(*candidateIterator),
                                               listIterator->first,
                                               currentWeight,
                            // we know that ( halfTotalWeight > currentWeight )
                                               ( weightSum - currentWeight ) );
                }
              }
              // end of checking if the weight of this combination of
              // *candidateIterator & *listIterator makes it worth keeping,
              // either as bestCandidateSoFar or just as a candidate to build
              // others from.
            }
            // end of if the candidate would not consist of too many objects.
          }
          // end of for loop over existing candidate partitions to build new
          // partitions from with (*listIterator)->first

          candidatesToBuildOn.splice( candidatesToBuildOn.begin(),
                                      extraCandidates );
          // all the new candidate partitions are moved from extraCandidates to
          // candidatesToBuildOn, at the start, though that isn't important.

          // now we need to build all the 2-object candidates with
          // *listIterator & all objects before it in orderedList:
          lowerWeightIterator = orderedList.begin();
          while( listIterator != lowerWeightIterator )
          {
            currentWeight
            = ( listIterator->second + lowerWeightIterator->second );
            if( halfTotalWeight < currentWeight )
              // if this candidate is already over half the weight, it cannot
              // be used to build upon, so we only keep it if it's a better
              // candidate than bestCandidateSoFar.
            {
              if( higherWeight > currentWeight )
                // if this *is* a better candidate...
              {
                // note the weights of this new current best candidate:
                higherWeight = currentWeight;
                lowerWeight = ( weightSum - currentWeight );
                // make the new best candidate & manage the memory:
                bestCandidateSoFar
                = &(vectorManagingCandidateMemory.newEnd());
                // set the candidate's properties:
                bestCandidateSoFar->getVector().assign( 2,
                                                  lowerWeightIterator->first );
                bestCandidateSoFar->getVector()[ 1 ] = listIterator->first;
                bestCandidateSoFar->setWeight( higherWeight );
                bestCandidateSoFar->setInverseWeight( lowerWeight );
              }
              // otherwise, we don't do anything with this combination of
              // *lowerWeightIterator & *listIterator, or any other
              // *lowerWeightIterator, because they'll all have even greater
              // weight, so we end the loop:
              lowerWeightIterator = listIterator;
            }
            else
              // otherwise, we have a new candidate to build upon, even if it's
              // not the best candidate so far...
            {
              // make the new candidate & manage the memory:
              currentCandidate = &(vectorManagingCandidateMemory.newEnd());
              candidatesToBuildOn.push_back( currentCandidate );
              // set the candidate's properties:
              currentCandidate->getVector().assign( 2,
                                                  lowerWeightIterator->first );
              currentCandidate->getVector()[ 1 ] = listIterator->first;
              currentCandidate->setWeight( currentWeight );
              currentCandidate->setInverseWeight( ( weightSum
                                                    - currentWeight ) );
              if( lowerWeight < currentWeight )
                // if this *is* a better candidate...
              {
                // note the weights of this new current best candidate:
                lowerWeight = currentWeight;
                higherWeight = ( weightSum - currentWeight );
                // note the actual new best current candidate:
                bestCandidateSoFar = currentCandidate;
              }

              // now we try the next-lowest weight pairing of *listIterator:
              ++lowerWeightIterator;
            }
            // end of if-else deciding whether to create a new candidate & move
            // on to the next pairing, or to break the loop (possibly after
            // updating bestCandidateSoFar).
          }
          // end of while loop over objects in orderedList before
          // *listIterator.

          // move on to the next object by weight:
          ++listIterator;
        }
        // end of while loop going over the higher-weighted objects, building
        // all the partition candidates which can involve them.

        // at this point, bestCandidateSoFar *is* the best candidate:
        if( halfTotalWeight > bestCandidateSoFar->getWeight() )
        {
          fillVectorFromBestCandidate( lowerWeightPartition );
          makeOtherPartition( higherWeightPartition,
                              lowerWeightPartition );
        }
        else
        {
          fillVectorFromBestCandidate( higherWeightPartition );
          makeOtherPartition( lowerWeightPartition,
                              higherWeightPartition );
        }
      }
    }
    // end of if-else about whether we needed to check all appropriate
    // combinations or not.
  }

  template< typename ClassToBeSorted >
  inline std::vector< unsigned int > const&
  BalancedPartitioner< ClassToBeSorted >::getLowerWeightPartitionIndices(
                                                                        ) const
  // this returns a reference to an internally-held std::vector of
  // unsigned ints corresponding to indices in the vector given to the last
  // call of makeMinimumDifferencePartition(...). the collection of objects
  // from that vector with these indices is the partition with GREATEST weight
  // UNDER (or equal to) half the the total weight of the set of objects last
  // given to makeMinimumDifferencePartition(...).
  // NOTE BENE: this std::vector gets OVER-WRITTEN by subsequent calls to
  // makeMinimumDifferencePartition(...)! if you want to keep the returned
  // set, you have to make your own *copy* of the std::vector referenced by
  // the return value of this function.
  {
    return lowerWeightPartition;
  }

  template< typename ClassToBeSorted >
  inline double
  BalancedPartitioner< ClassToBeSorted >::getLowerWeight() const
  // this returns the sum of absolute weights of the objects with indices in
  // the return reference of getLowerWeightPartitionIndices().
  {
    return lowerWeight;
  }

  template< typename ClassToBeSorted >
  inline std::vector< unsigned int > const&
  BalancedPartitioner< ClassToBeSorted >::getHigherWeightPartitionIndices(
                                                                        ) const
  // this returns a reference to an internally-held std::vector of
  // unsigned ints corresponding to indices in the vector given to the last
  // call of makeMinimumDifferencePartition(...). the collection of objects
  // from that vector with these indices is the partition with LEAST weight
  // OVER (or equal to) half the the total weight of the set of objects last
  // given to makeMinimumDifferencePartition(...).
  // NOTE BENE: this std::vector gets OVER-WRITTEN by subsequent calls to
  // makeMinimumDifferencePartition(...)! if you want to keep the returned
  // set, you have to make your own *copy* of the std::vector referenced by
  // the return value of this function.
  {
    return higherWeightPartition;
  }

  template< typename ClassToBeSorted >
  inline double
  BalancedPartitioner< ClassToBeSorted >::getHigherWeight() const
  // this returns the sum of absolute weights of the objects with indices in
  // the return reference of getHigherWeightPartitionIndices().
  {
    return higherWeight;
  }

  template< typename ClassToBeSorted >
  inline void
  BalancedPartitioner< ClassToBeSorted >::clearSavedCandidates()
  {
    lowerWeightPartition.clear();
    higherWeightPartition.clear();
    orderedList.clear();
    vectorManagingCandidateMemory.clearEntries();
    candidatesToBuildOn.clear();
  }

  template< typename ClassToBeSorted >
  inline bool
  BalancedPartitioner< ClassToBeSorted >::orderByDouble(
                                              IndexWithWeight const& firstPair,
                                            IndexWithWeight const& secondPair )
  // this returns true if firstPair->second <= secondPair->second, false
  // otherwise.
  {
    if( firstPair.second <= secondPair.second )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  template< typename ClassToBeSorted >
  inline void
  BalancedPartitioner< ClassToBeSorted >::fillVectorFromBestCandidate(
                               std::vector< unsigned int >& partitionToBeMade )
  {
    partitionToBeMade = bestCandidateSoFar->getVector();
  }

  template< typename ClassToBeSorted >
  inline void
  BalancedPartitioner< ClassToBeSorted >::makeOtherPartition(
                                std::vector< unsigned int >& partitionToBeMade,
                            std::vector< unsigned int >& partitionAlreadyMade )
  // this REMOVES ALL OF THE PAIRS WITH THE SAME INDICES AS IN
  // partitionAlreadyMade FROM orderedList & then fills partitionToBeMade with
  // the indices that remain.
  {
    IndexWithWeightListIterator removingIterator( orderedList.begin() );
    for( CandidateVectorConstIterator
         alreadyUsedIterator( partitionAlreadyMade.begin() );
         partitionAlreadyMade.end() > alreadyUsedIterator;
         ++alreadyUsedIterator )
    {
    // now we look for *alreadyUsedIterator in orderedList:
      removingIterator = orderedList.begin();
      while( removingIterator->first != *alreadyUsedIterator )
      {
        ++removingIterator;
      }
      // now we've found *alreadyUsedIterator in orderedList.
      orderedList.erase( removingIterator );
      // we remove it from the list.
    }
    // now all the elements of partitionAlreadyMade should have been removed
    // from orderedList.
    for( IndexWithWeightListIterator fillingIterator( orderedList.begin() );
         orderedList.end() != fillingIterator;
         ++fillingIterator )
    {
      partitionToBeMade.push_back( fillingIterator->first );
    }
  }

}  // end of BOL namespace

#endif /* BALANCEDPARTITIONER_HPP_ */

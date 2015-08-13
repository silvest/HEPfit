/*
 * BalancedPartitionCandidate.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef BALANCEDPARTITIONCANDIDATE_HPP_
#define BALANCEDPARTITIONCANDIDATE_HPP_

#include <vector>
#include <list>
#include "UsefulStuff.hpp"
#include "VectorlikeArray.hpp"

namespace BOL
{
  // this class is used by the BalancedPartitioner class to hold candidates for
  // the partition.
  class BalancedPartitionCandidate
  {
  public:
    BalancedPartitionCandidate();
    ~BalancedPartitionCandidate();

    void
    buildFrom( BalancedPartitionCandidate const& basePartition,
               unsigned int const extraIndex,
               double const candidateWeight,
               double const inverseWeight );
    // this sets up this BalancedPartitionCandidate to be
    // basePartition + extraIndex, using the given weights (rather than
    // working them out itself).

    std::vector< unsigned int >&
    getVector();
    double
    getWeight() const;
    void
    setWeight( double const candidateWeight );
    double
    getInverseWeight() const;
    void
    setInverseWeight( double const inverseWeight );


  protected:
    std::vector< unsigned int > indexVector;
    double candidateWeight;
    double inverseWeight;
  };





  inline void
  BalancedPartitionCandidate::buildFrom(
                               BalancedPartitionCandidate const& basePartition,
                                         unsigned int const extraIndex,
                                         double const candidateWeight,
                                         double const inverseWeight )
  // this sets up this BalancedPartitionCandidate to be
  // basePartition + extraIndex, using the given weights (rather than working
  // them out itself).
  {
    indexVector = basePartition.indexVector;
    indexVector.push_back( extraIndex );
    this->candidateWeight = candidateWeight;
    this->inverseWeight = inverseWeight;
  }

  inline std::vector< unsigned int >&
  BalancedPartitionCandidate::getVector()
  {
    return indexVector;
  }

  inline double
  BalancedPartitionCandidate::getWeight() const
  {
    return candidateWeight;
  }

  inline void
  BalancedPartitionCandidate::setWeight( double const candidateWeight )
  {
    this->candidateWeight = candidateWeight;
  }

  inline double
  BalancedPartitionCandidate::getInverseWeight() const
  {
    return inverseWeight;
  }

  inline void
  BalancedPartitionCandidate::setInverseWeight( double const inverseWeight )
  {
    this->inverseWeight = inverseWeight;
  }

}

#endif /* BALANCEDPARTITIONCANDIDATE_HPP_ */

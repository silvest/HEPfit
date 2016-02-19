/*  
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MATCHING_H
#define MATCHING_H

#include <boost/ref.hpp>

template <typename T, typename U>
class Matching {
public:
  Matching(const U& ur): obj(ur),objr(boost::ref(obj)) {}
  T& getObj() { return objr; }
  void setObj(T& obji) {objr = boost::ref(obji);}
private:
  T obj;
  boost::reference_wrapper<T> objr;  	
};

#endif /* MATCHING_H */


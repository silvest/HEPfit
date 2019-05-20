/*  
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MATCHING_H
#define MATCHING_H

template <typename T, typename U>
class Matching {
public:
  Matching(const U& ur): obj(ur),objr(std::ref(obj)) {}
  T& getObj() { return objr; }
  void setObj(T& obji) {objr = std::ref(obji);}
private:
  T obj;
  std::reference_wrapper<T> objr;  	
};

#endif /* MATCHING_H */


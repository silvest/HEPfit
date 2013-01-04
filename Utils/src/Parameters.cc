/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <iostream>
#include "Parameters.h"

const std::string Parameters::TypeList[NumberOfTypes] = {"int", "double",
"complex", "string", "gslpp::matrix<double>", "gslpp::matrix<gslpp::complex>",
"std::vector<double>"};

Parameters::Parameters(Parameters& P)
{
  Ints = P.getInts();
  Doubles = P.getDoubles();
  Complexes = P.getComplexes();
  Strings = P.getStrings();
  ComplexMatrices = P.getComplexMatrices();
  DoubleMatrices = P.getDoubleMatrices();
  DoubleVectors = P.getDoubleVectors();
}

std::map<std::string, int> Parameters::getInts() const
{
  return Ints;
}

std::map<std::string, double> Parameters::getDoubles() const
{
  return Doubles;
}

std::map<std::string, gslpp::complex> Parameters::getComplexes() const
{
  return Complexes;
}

std::map<std::string, std::string> Parameters::getStrings() const
{
  return Strings;
}

std::map<std::string, gslpp::matrix<double> > Parameters::getDoubleMatrices() const
{
  return DoubleMatrices;
}

std::map<std::string, gslpp::matrix<gslpp::complex> > Parameters::getComplexMatrices() const
{
  return ComplexMatrices;
}

std::map<std::string, std::vector<double> > Parameters::getDoubleVectors() const
{
  return DoubleVectors;
}

void Parameters::Set(std::string s, int i)
{
  InOtherMaps(s, INT);
  Ints[s] = i;
}

void Parameters::Set(std::string s, double d)
{
  InOtherMaps(s, DOUBLE);
  Doubles[s] = d;
}

void Parameters::Set(std::string s, gslpp::complex z)
{
  InOtherMaps(s, COMPLEX);
  Complexes[s] = z;
}

void Parameters::Set(std::string s, std::string t)
{
  InOtherMaps(s, STRING);
  Strings[s] = t;
}

void Parameters::Set(std::string s, gslpp::matrix<double> md)
{
  InOtherMaps(s, DOUBLE_MATRIX);
  DoubleMatrices.insert(DoubleMatrices.find(s),std::make_pair(s,md));
}

void Parameters::Set(std::string s, gslpp::matrix<gslpp::complex> mc)
{
  InOtherMaps(s, COMPLEX_MATRIX);
  ComplexMatrices.insert(ComplexMatrices.find(s),std::make_pair(s,mc));
}

void Parameters::Set(std::string s, std::vector<double> v)
{
  InOtherMaps(s, DOUBLE_VECTOR);
  DoubleVectors.insert(DoubleVectors.find(s),std::make_pair(s,v));
}

void Parameters::Get(std::string s, int& i) const
{
  InMap(s, INT);
  i = Ints.find(s)->second;
}

void Parameters::Get(std::string s, double& d) const
{
  InMap(s, DOUBLE);
  d = Doubles.find(s)->second;
}

void Parameters::Get(std::string s, gslpp::complex& z) const
{
  InMap(s, COMPLEX);
  z = Complexes.find(s)->second;
}

void Parameters::Get(std::string s, std::string& t) const
{
  InMap(s, STRING);
  t = Strings.find(s)->second;
}

void Parameters::Get(std::string s, gslpp::matrix<double> & md) const
{
  InMap(s, DOUBLE_MATRIX);
  md = DoubleMatrices.find(s)->second;
}

void Parameters::Get(std::string s, gslpp::matrix<gslpp::complex> & mc) const
{
  InMap(s, COMPLEX_MATRIX);
  mc = ComplexMatrices.find(s)->second;
}

void Parameters::Get(std::string s, std::vector<double> & v) const
{
  InMap(s, DOUBLE_VECTOR);
  v = DoubleVectors.find(s)->second;
}

void Parameters::InOtherMaps(std::string s, MapType m) const
{
  int n = -1;

  if(m != INT && Ints.find(s) != Ints.end()) n = INT;
  else if(m != DOUBLE && Doubles.find(s) != Doubles.end()) n = DOUBLE;
  else if(m != COMPLEX && Complexes.find(s) != Complexes.end()) n = COMPLEX;
  else if(m != STRING && Strings.find(s) != Strings.end()) n = STRING;
  else if(m != DOUBLE_MATRIX && DoubleMatrices.find(s) != DoubleMatrices.end())
      n = DOUBLE_MATRIX;
  else if(m != COMPLEX_MATRIX && ComplexMatrices.find(s) != ComplexMatrices.end())
      n = COMPLEX_MATRIX;
  else if(m != DOUBLE_VECTOR && DoubleVectors.find(s) != DoubleVectors.end())
      n = DOUBLE_VECTOR;
  if(n != -1)
    {
      std::cout << "ERROR: requested "<< TypeList[m] << " key \""
      << s << "\" already exists as " << TypeList[n] << std::endl;
      exit(EXIT_FAILURE);
    }
}

void Parameters::InMap(std::string s, MapType m) const
{
  switch(m)
  {
    case INT:
      if(Ints.find(s) != Ints.end()) return;
    case DOUBLE:
      if(Doubles.find(s) != Doubles.end()) return;
    case COMPLEX:
      if(Complexes.find(s) != Complexes.end()) return;
    case STRING:
      if(Strings.find(s) != Strings.end()) return;
    case DOUBLE_MATRIX:
      if(DoubleMatrices.find(s) != DoubleMatrices.end()) return;
    case COMPLEX_MATRIX:
      if(ComplexMatrices.find(s) != ComplexMatrices.end()) return;
    case DOUBLE_VECTOR:
      if(DoubleVectors.find(s) != DoubleVectors.end()) return;
  }
  std::cout << "ERROR: wrong type or key \"" << s << "\" not found" << std::endl;
  exit(EXIT_FAILURE);
}

int Parameters::Find(std::string s) const {
    if(Ints.find(s) != Ints.end()) return(INT);
    if(Doubles.find(s) != Doubles.end()) return(DOUBLE);
    if(Complexes.find(s) != Complexes.end()) return(COMPLEX);
    if(Strings.find(s) != Strings.end()) return(STRING);
    if(DoubleMatrices.find(s) != DoubleMatrices.end()) return(DOUBLE_MATRIX);
    if(ComplexMatrices.find(s) != ComplexMatrices.end()) return(COMPLEX_MATRIX);
    if(DoubleVectors.find(s) != DoubleVectors.end()) return(DOUBLE_VECTOR);
    return(-1);
}

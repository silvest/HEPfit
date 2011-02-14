#include <cstdlib>
#include <iostream>
#include "Parameters.h"

const std::string Parameters::TypeList[NumberOfTypes] = {"int", "double",
"complex", "string", "gslpp::matrix<double>", "gslpp::matrix<gslpp::complex>"};

Parameters::Parameters(Parameters& P)
{
  Ints = P.getInts();
  Doubles = P.getDoubles();
  Complexes = P.getComplexes();
  Strings = P.getStrings();
  ComplexMatrices = P.getComplexMatrices();
  DoubleMatrices = P.getDoubleMatrices();
}

std::map<std::string, int> Parameters::getInts()
{
  return Ints;
}

std::map<std::string, double> Parameters::getDoubles()
{
  return Doubles;
}

std::map<std::string, gslpp::complex> Parameters::getComplexes()
{
  return Complexes;
}

std::map<std::string, std::string> Parameters::getStrings()
{
  return Strings;
}

std::map<std::string, gslpp::matrix<double> > Parameters::getDoubleMatrices()
{
  return DoubleMatrices;
}

std::map<std::string, gslpp::matrix<gslpp::complex> > Parameters::getComplexMatrices()
{
  return ComplexMatrices;
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

void Parameters::Get(std::string s, int& i)
{
  InMap(s, INT);
  i = Ints[s];
}

void Parameters::Get(std::string s, double& d)
{
  InMap(s, DOUBLE);
  d = Doubles[s];
}

void Parameters::Get(std::string s, gslpp::complex& z)
{
  InMap(s, COMPLEX);
  z = Complexes[s];
}

void Parameters::Get(std::string s, std::string& t)
{
  InMap(s, STRING);
  t = Strings[s];
}

void Parameters::Get(std::string s, gslpp::matrix<double> & md)
{
  InMap(s, DOUBLE_MATRIX);
  md = DoubleMatrices.find(s)->second;
}

void Parameters::Get(std::string s, gslpp::matrix<gslpp::complex> & mc)
{
  InMap(s, COMPLEX_MATRIX);
  mc = ComplexMatrices.find(s)->second;
}

void Parameters::InOtherMaps(std::string s, MapType m)
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
  if(n != -1)
    {
      std::cout << "ERROR: requested "<< TypeList[m] << " key \""
      << s << "\" already exists as " << TypeList[n] << std::endl;
      exit(EXIT_FAILURE);
    }
}

void Parameters::InMap(std::string s, MapType m)
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
  }
  std::cout << "ERROR: wrong type or key \"" << s << "\" not found" << std::endl;
  exit(EXIT_FAILURE);
}

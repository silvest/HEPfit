/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Expanded.h"

template <class T>
Expanded<T>::Expanded() {
    minord1 = 0;
    minord2 = 0;
    n1 = 0;
}

template Expanded<double>::Expanded();
template Expanded<complex>::Expanded();
template Expanded<vector<double> >::Expanded();
template Expanded<vector<complex> >::Expanded();
template Expanded<matrix<double> >::Expanded();
template Expanded<matrix<complex> >::Expanded();

template <class T>
Expanded<T>::Expanded(std::vector<T> & dinp, int minord2_i)
{
    minord1 = 0;
    n1 = 1;
    minord2 = minord2_i;
    n2.push_back(dinp.size());
    data.push_back(dinp);
}

template Expanded<double>::Expanded(std::vector<double> & dinp, int minord2_i);
template Expanded<complex>::Expanded(std::vector<complex> & dinp, int minord2_i);
template Expanded<vector<double> >::Expanded(std::vector<vector<double> > & dinp, int minord2_i);
template Expanded<vector<complex> >::Expanded(std::vector<vector<complex> > & dinp, int minord2_i);
template Expanded<matrix<double> >::Expanded(std::vector<matrix<double> > & dinp, int minord2_i);
template Expanded<matrix<complex> >::Expanded(std::vector<matrix<complex> > & dinp, int minord2_i);

template <class T>
Expanded<T>::Expanded(std::vector<std::vector<T> > & dinp, int minord2_i, int minord1_i) {
        minord1 = minord1_i;
        minord2 = minord2_i;
        n1 = dinp.size();
        for (int i = 0; i < n1; i++)
            n2.push_back(dinp[i].size());
        data = dinp;
    }

template Expanded<double>::Expanded(std::vector<std::vector<double> >& dinp, int minord2_i, int minord1_i);
template Expanded<complex>::Expanded(std::vector<std::vector<complex> >& dinp, int minord2_i, int minord1_i);
template Expanded<vector<double> >::Expanded(std::vector<std::vector<vector<double> > >& dinp, int minord2_i, int minord1_i);
template Expanded<vector<complex> >::Expanded(std::vector<std::vector<vector<complex> > >& dinp, int minord2_i, int minord1_i);
template Expanded<matrix<double> >::Expanded(std::vector<std::vector<matrix<double> > >& dinp, int minord2_i, int minord1_i);
template Expanded<matrix<complex> >::Expanded(std::vector<std::vector<matrix<complex> > >& dinp, int minord2_i, int minord1_i);

template <> Expanded<complex> Expanded<complex>::conjugate() const
{
    std::vector<std::vector<complex> > res;
    for (int i = 0; i < n1; i++)
    {
        res.push_back(std::vector<complex>(n2.at(i)));
        for (int j = 0; j < n2.at(i); j++)
            res[i][j] = data[i][j].conjugate();
    }
    return Expanded<complex>(res, minord2, minord1);
}

template <> Expanded<matrix<double> > Expanded<matrix<double> >::transpose() const
{
    std::vector<std::vector<matrix<double> > > res;
    for (int i1 = 0; i1 < n1; i1++)
    {
        res.push_back(std::vector<matrix<double> >(n2.at(i1), matrix<double>(data[0][0].size_i(), data[0][0].size_j(), 0.)));
        for (int i2 = 0; i2 < n2.at(i1); i2++)
            res[i1][i2] = data[i1][i2].transpose();
    }
    return Expanded<matrix<double>>(res, minord2, minord1);
}

template <> Expanded<matrix<complex> > Expanded<matrix<complex> >::transpose() const
{
    std::vector<std::vector<matrix<complex> > > res;
    for (int i1 = 0; i1 < n1; i1++)
    {
        res.push_back(std::vector<matrix<complex> >(n2.at(i1), matrix<complex>(data[0][0].size_i(), data[0][0].size_j(), 0.)));
        for (int i2 = 0; i2 < n2.at(i1); i2++)
            res[i1][i2] = data[i1][i2].transpose();
    }
    return Expanded<matrix < complex >> (res, minord2, minord1);
}

template <> Expanded<double> Expanded<complex>::real() const
{
    std::vector<std::vector<double> > res;
    for (int i = 0; i < n1; i++)
    {
        res.push_back(std::vector<double>(n2.at(i), 0.));
        for (int j = 0; j < n2.at(i); j++)
            res[i][j] = data[i][j].real();
    }
    return Expanded<double>(res, minord2, minord1);
}

template <> Expanded<double> Expanded<complex>::abs2() const
{
    return ((*this) * ((*this).conjugate())).real();
}

template <>
matrix<double> Expanded<matrix<double> >::Zero() const
{
    int id, jd;
    id = data[0][0].size_i();
    jd = data[0][0].size_j();
    return matrix<double>(id, jd, 0.);
}

template <>
matrix<complex> Expanded<matrix<complex> >::Zero() const
{
    int id, jd;
    id = data[0][0].size_i();
    jd = data[0][0].size_j();
    return matrix<complex>(id, jd, 0.);
}

template<>
vector<double> Expanded<vector<double> >::Zero() const
{
    int id;
    id = data[0][0].size();
    return vector<double>(id, 0.);
}

template<>
vector<complex> Expanded<vector<complex> >::Zero() const
{
    int id;
    id = data[0][0].size();
    return vector<complex>(id, 0.);
}

template<>
double Expanded<double>::invCoeff(int i, int j) const {
    if(i == 0 && j == 0) return (1. / data[0][0]);
    double res = 0.;
    for(int n = 0; n <= i; n++)
        for(int m = 0; m <= j; m++) {
            if(n == i && m == j) continue;
            res -= invCoeff(n, m) * data[i - n][j - m];
        }
    return (res / data[0][0]);
}

template<>
complex Expanded<complex>::invCoeff(int i, int j) const {
    if(i == 0 && j == 0) return (1. / data[0][0]);
    complex res = 0.;
    for(int n = 0; n <= i; n++)
        for(int m = 0; m <= j; m++) {
            if(n == i && m == j) continue;
            res -= invCoeff(n, m) * data[i - n][j - m];
        }
    return (res / data[0][0]);
}

template <>
Expanded<double> Expanded<double>::inverse() const
{
    std::vector<std::vector<double> > res;
    for (int i1 = 0; i1 < n1; i1++)
    {
        res.push_back(std::vector<double> (n2.at(i1), 0.));
        for (int i2 = 0; i2 < n2.at(i1); i2++)
            res[i1][i2] = invCoeff(i1, i2);
    }
    return Expanded<double>(res, -minord2, -minord1);
}

template <>
Expanded<complex> Expanded<complex>::inverse() const
{
    std::vector<std::vector<complex> > res;
    for (int i1 = 0; i1 < n1; i1++)
    {
        res.push_back(std::vector<complex> (n2.at(i1), 0.));
        for (int i2 = 0; i2 < n2.at(i1); i2++)
            res[i1][i2] = invCoeff(i1, i2);
    }
    return Expanded<complex>(res, -minord2, -minord1);
}

/***************** Begin UnExpanded * Expanded *****************/

// scalar * X -> commute

Expanded<complex> operator*(const complex& ue, const Expanded<double>& ex)
{
    return (ex * ue);
}

Expanded<vector<complex> > operator*(const complex& ue, const Expanded<vector<double> >& ex)
{
    return (ex * ue);
}

Expanded<matrix<complex> > operator*(const complex& ue, const Expanded<matrix<double> >& ex)
{
    return (ex * ue);
}

// X * scalar -> commute

Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<double>& ex)
{
    return (ex * ue);
}

Expanded<vector<complex> > operator*(const vector<double>& ue, const Expanded<complex>& ex)
{
    return (ex * ue);
}

Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<double>& ex)
{
    return (ex * ue);
}

Expanded<matrix<complex> > operator*(const matrix<double>& ue, const Expanded<complex>& ex)
{
    return (ex * ue);
}

// both vector -> commute

Expanded<double> operator*(const vector<double>& ue, const Expanded<vector<double> >& ex)
{
    return (ex * ue);
}

//  vector * matrix

Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<matrix<double> >& ex)
{
    return (ex.transpose() * ue);
}

//  matrix * vector

Expanded<vector<double> > operator*(const matrix<double>& ue, const Expanded<vector<double> >& ex)
{
    return (ex * ue.transpose());
}

//  matrix * matrix

Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<matrix<double> >& ex)
{
    return ((ex.transpose() * ue.transpose()).transpose());
}

/***************** End UnExpanded * Expanded *****************/

/* 
 * File:   EvolBsmm.h
 * Author: claudio
 *
 * Created on 30 settembre 2015, 18.44
 */

#ifndef EVOLBSMM_H
#define	EVOLBSMM_H

#include <RGEvolutor.h>
#include <StandardModel.h>
#include <sstream>

class EvolBsmm : public RGEvolutor {
public:
    EvolBsmm(unsigned int dim, schemes scheme, orders order, orders_ew order_ew, const StandardModel& model);
    virtual ~EvolBsmm();
    gslpp::matrix<double> AnomalousDimension(int gam, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double>& Df1Evol(double mu, double M, orders order, orders_ew order_ew, schemes scheme = NDR);
    double alphatilde_e(double mu);
    double alphatilde_s(double mu);


private:
    int nu, nd;
    double a[4][8], b[4][8][8][8];
    const StandardModel& model;
    void Df1Evol(double mu, double M, double nf, schemes scheme);
    gslpp::matrix <gslpp::complex> V, Vi, AA, BB, CC, DD, EE, FF, RR;
    gslpp::vector<gslpp::complex> e;
    std::vector<double> vavi, vbvi, vcvi, vdvi, vevi, vfvi, vrvi, vaevi, vbbvi,
                        vbdvi, vbevi, vdbvi, vdevi, veavi, vebvi, vedvi, veevi, vbeevi,
                        vebevi, veebvi, vbbevi, vbebvi, vebbvi;
    unsigned int dim;
    double alsM;
    double alsmu;

    double eta;
    double logeta;

    double F(unsigned int i, unsigned int j, int x, double mu, double M, double nf);
    double R(unsigned int i, unsigned int j, int x, double mu, double M, double nf);
    double G(unsigned int i, unsigned int p, unsigned int j, int x, int y, double mu, double M, double nf);
    double H(unsigned int i, unsigned int p, unsigned int q, unsigned int j, int x, int y, int z, double mu, double M, double nf);
    gslpp::matrix<double> BuiltB(char letter, unsigned int n_u, unsigned int n_d);
};

#endif	/* EVOLBSMM_H */

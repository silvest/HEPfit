/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <map>
#include <vector>
#include <initializer_list>
#include "QCD.h"
#include "EvolDF1.h"

#define EPS 1.e-15

std::map<std::string, uint> blocks_nops = {
    {"C", 2},
    {"CP", 6},
    {"CPM", 8},
    {"L", 2},
    {"CPL", 8},
    {"CPQ", 10},
    {"CPML", 10},
    {"CPQB", 11},
    {"CPMQB", 13},
    {"CPMLQB", 15}
};

EvolDF1::EvolDF1(std::string reqblocks, schemes scheme, const StandardModel& model_i, qcd_orders ord_qcd, qed_orders ord_qed)
: RGEvolutorNew(blocks_nops.at(reqblocks), scheme, ord_qcd, ord_qed), model(model_i), blocks(reqblocks),
evec(blocks_nops.at(reqblocks), 0.), evec_i(blocks_nops.at(reqblocks), 0.), js(blocks_nops.at(reqblocks), 0.),
h(blocks_nops.at(reqblocks), 0.), gg(blocks_nops.at(reqblocks), 0.), s_s(blocks_nops.at(reqblocks), 0.),
jssv(blocks_nops.at(reqblocks), 0.), jss(blocks_nops.at(reqblocks), 0.), jv(blocks_nops.at(reqblocks), 0.),
vij(blocks_nops.at(reqblocks), 0.), eval(blocks_nops.at(reqblocks), 0.), evecc(blocks_nops.at(reqblocks), 0.),
evalc(blocks_nops.at(reqblocks), 0.)
{
    //    blocks_ord = {
    //        {"C", NNLO},
    //        {"CP", NNLO},
    //        {"CPM", NNLO},
    //        {"L", NNLO},
    //        {"CPML", NNLO},
    //        {"CPQB", NLO},
    //        {"CPMQB", NLO},
    //        {"CPMLQB", NLO}};

    //   if (blocks_nops[blocks] != nops)
    //        throw std::runtime_error("EvolDF1(): number of operators does not match block specification");

    this->nops = blocks_nops.at(reqblocks);
    uint nf, nnf, nu, nd, a, b, i, j, p, q;
    double b0, b0e, b1, b2, b3, b4, term;

    alsM_cache = 0.;
    MAls_cache = 0.;

    gslpp::matrix<double> W10(nops, nops, 0.), W20(nops, nops, 0.), W30(nops, nops, 0.),
            W01(nops, nops, 0.), W02(nops, nops, 0.), W11(nops, nops, 0.), W21(nops, nops, 0.);
    gslpp::matrix<double> M1(nops, nops, 0.), M2(nops, nops, 0.), M3(nops, nops, 0.), M4(nops, nops, 0.),
            M5(nops, nops, 0.), M6(nops, nops, 0.);

    if (order_qed == QED0 && blocks.find("L") == std::string::npos &&
            blocks.find("Q") == std::string::npos && blocks.find("B") == std::string::npos)
    {
        nfmin = 3;
        nfmax = 6;
    } else
        nfmin = nfmax = 5;


    for (nf = nfmin; nf <= nfmax; nf++)
    {
        nnf = nf - nfmin;
        nu = nf / 2;
        nd = nf % 2 == 0 ? nf / 2 : nf / 2 + 1;

        b0 = model.Beta_s(00, nf);
        b1 = model.Beta_s(10, nf) / 2. / b0 / b0;
        b2 = model.Beta_s(20, nf) / 4. / b0 / b0 / b0 - b1 * b1;

        W10 = AnomalousDimension(10, nu, nd).transpose() / 2. / b0;
        W20 = AnomalousDimension(20, nu, nd).transpose() / 4. / b0 / b0;
        W30 = AnomalousDimension(30, nu, nd).transpose() / 8. / b0 / b0 / b0;

        //        std::cout << AnomalousDimension(01, nu, nd).transpose() << std::endl;
        //        std::cout << W10 << std::endl;

        // Misiak-Munz basis, defined as in T. Huber et al., hep-ph/0512066
        W10.eigensystem(evecc, evalc);
        for (size_t jj = 0; jj < evalc.size(); jj++)
            if (fabs(evalc(jj).imag()) > EPS) throw ("check the imaginary part of eigenvalues of W10");
        eval = evalc.real();
        evec = evecc.real();
        evec_i = (evecc.inverse()).real();

        // QCD magic numbers
        // M2: B(-2), M1: B(-1)_10
        M2 = evec_i * (W30 - b1 * W20 - b2 * W10) * evec;
        M1 = evec_i * (W20 - b1 * W10) * evec;

        for (a = 0; a < nops; a++)
        {
            ai[nnf].insert(std::pair<uint, double > (a, eval(a)));
            for (b = 0; b < nops; b++)
                for (i = 0; i < nops; i++)
                {
                    if (fabs(term = evec(a, i) * evec_i(i, b)) > EPS)
                        vM0vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i}, term)); // QCD LO evolutor
                    for (j = 0; j < nops; j++)
                    {
                        if (fabs(term = evec(a, i) * M1(i, j) * evec_i(j, b)) > EPS)
                            vM1vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j}, term)); // QCD NLO evolutor
                        if (fabs(term = evec(a, i) * M2(i, j) * evec_i(j, b)) > EPS)
                            vM2vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j}, term)); // QCD NNLO evolutor
                        for (p = 0; p < nops; p++)
                            if (fabs(term = evec(a, i) * M1(i, p) * M1(p, j) * evec_i(j, b)) > EPS)
                                vM11vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term)); // QCD NNLO evolutor                        
                    }
                }
        }

        if (order_qed != QED0)
        {
            b0e = model.Beta_e(00, nf);
            b3 = model.Beta_s(01, nf) / 2. / b0 / b0e;
            b4 = model.Beta_s(11, nf) / 4. / b0 / b0 / b0e - 2. * b1*b3;
            W01 = AnomalousDimension(01, nu, nd).transpose() / 2. / b0e;
            W02 = AnomalousDimension(02, nu, nd).transpose() / 4. / b0e / b0e;
            W11 = AnomalousDimension(11, nu, nd).transpose() / 4. / b0 / b0e;
            W21 = AnomalousDimension(21, nu, nd).transpose() / 8. / b0 / b0 / b0e;

            // QED magic numbers
            // M3: B(1)_01, B(2)_02, B(1)_02, R5_12, M4: B(0)_11, B(0)_12, M5: B(-1)_21, M6: B(1)_12 - proper powers of omega and lambda added in DF1Evol
            M3 = evec_i * W01 * evec;
            M4 = evec_i * (W11 - b1 * W01 - b3 * W10) * evec;
            M5 = evec_i * (W21 - b1 * W11 - b2 * W01 - b3 * W20 - b4 * W10) * evec;
            M6 = evec_i * (W02 + W11 - (b1 + b3) * W01 - b3 * W10) * evec;
            for (a = 0; a < nops; a++)
                for (b = 0; b < nops; b++)
                    for (i = 0; i < nops; i++)
                        for (j = 0; j < nops; j++)
                        {
                            if (fabs(term = evec(a, i) * M3(i, j) * evec_i(j, b)) > EPS)
                                vM3vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j}, term));
                            if (fabs(term = evec(a, i) * M4(i, j) * evec_i(j, b)) > EPS)
                                vM4vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j}, term));
                            if (fabs(term = evec(a, i) * M5(i, j) * evec_i(j, b)) > EPS)
                                vM5vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j}, term));
                            if (fabs(term = evec(a, i) * M6(i, j) * evec_i(j, b)) > EPS)
                                vM6vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j}, term));
                            for (p = 0; p < nops; p++)
                            {
                                if (fabs(term = evec(a, i) * M3(i, p) * M3(p, j) * evec_i(j, b)) > EPS)
                                    vM33vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M3(i, p) * M1(p, j) * evec_i(j, b)) > EPS)
                                    vM31vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M1(i, p) * M3(p, j) * evec_i(j, b)) > EPS)
                                    vM13vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M3(i, p) * M4(p, j) * evec_i(j, b)) > EPS)
                                    vM34vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M4(i, p) * M3(p, j) * evec_i(j, b)) > EPS)
                                    vM43vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M2(i, p) * M3(p, j) * evec_i(j, b)) > EPS)
                                    vM23vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M3(i, p) * M2(p, j) * evec_i(j, b)) > EPS)
                                    vM32vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M1(i, p) * M4(p, j) * evec_i(j, b)) > EPS)
                                    vM14vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                if (fabs(term = evec(a, i) * M4(i, p) * M1(p, j) * evec_i(j, b)) > EPS)
                                    vM41vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p}, term));
                                for (q = 0; q < nops; q++)
                                {
                                    if (fabs(term = evec(a, i) * M1(i, p) * M1(p, q) * M3(q, j) * evec_i(j, b)) > EPS)
                                        vM113vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p, q}, term));
                                    if (fabs(term = evec(a, i) * M1(i, p) * M3(p, q) * M1(q, j) * evec_i(j, b)) > EPS)
                                        vM131vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p, q}, term));
                                    if (fabs(term = evec(a, i) * M3(i, p) * M1(p, q) * M1(q, j) * evec_i(j, b)) > EPS)
                                        vM311vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p, q}, term));
                                    if (fabs(term = evec(a, i) * M3(i, p) * M3(p, q) * M1(q, j) * evec_i(j, b)) > EPS)
                                        vM331vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p, q}, term));
                                    if (fabs(term = evec(a, i) * M3(i, p) * M1(p, q) * M3(q, j) * evec_i(j, b)) > EPS)
                                        vM313vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p, q}, term));
                                    if (fabs(term = evec(a, i) * M1(i, p) * M3(p, q) * M3(q, j) * evec_i(j, b)) > EPS)
                                        vM133vi[nnf].insert(std::pair<std::vector<uint>, double > ({a, b, i, j, p, q}, term));
                                }
                            }
                        }
        }
    }
}

EvolDF1::~EvolDF1()
{
}

/* Delta F = 1 anomalous dimension in Misiak basis, in the effective basis (C7eff, C8eff)
       ref. for CC,CP,PP QED NLO, QP,QQ,BP,BB QCD NLO, CL,PL NNLO, PQ,QL,LL,BL NLO: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
       ref. for CC,CP,PP QCD NNLO with nf: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
       ref. for CM,PM QCD NNLO with nf: Czakon, Haisch, Misiak, JHEP 0703, 008, hep-ph/0612329
       ref. for MM QCD NNLO with nf: Gorbahn, Haisch, Misiak, Phys. Rev. Lett. 95, 102004, hep-ph/0504194
       ref. for CM,PM QED LO, QM QCD LO: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
       ref. for MM,QP,QQ,BP,BB QED NLO, BQ NLO: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
                       QM QED,BM??? (in 031209 C7,8 are normalised with alpha_s and not in the effective basis)
 */

double EvolDF1::f_f(uint nnf, uint i, uint j, int k, double eta)
{
    for (int il = 0; il < F_iCacheSize; il++)
        if (f_f_c[0][il] == (int) nnf && f_f_c[1][il] == (int) i && f_f_c[2][il] == (int) j &&
                f_f_c[3][il] == k && f_f_d[0][il] == eta)
            return f_f_d[1][il];

    double aii = ai[nnf].at(i), aij = ai[nnf].at(j);
    double den = aij + k - aii, ret;

    if (fabs(den) < EPS)
        ret = pow(eta, aii) * log(eta);
    else
        ret = (pow(eta, aij + k) - pow(eta, aii)) / den;

    model.CacheShift(f_f_c, 4);
    model.CacheShift(f_f_d, 2);
    f_f_c[0][0] = (int) nnf;
    f_f_c[1][0] = (int) i;
    f_f_c[2][0] = (int) j;
    f_f_c[3][0] = k;
    f_f_d[0][0] = eta;
    f_f_d[1][0] = ret;

    return ret;
}

double EvolDF1::f_r(uint nnf, uint i, uint j, int k, double eta)
{
    double ll = log(eta), den = ai[nnf].at(j) + k - ai[nnf].at(i);

    if (fabs(den) < EPS)
        return (0.5 * pow(eta, ai[nnf].at(i)) * ll * ll);
    else
        return ((pow(eta, ai[nnf].at(j) + k) * ll - f_f(nnf, i, j, k, eta)) / den);
}

double EvolDF1::f_g(uint nnf, uint i, uint p, uint j, int k, int l, double eta)
{
    double den = ai[nnf].at(j) + l - ai[nnf].at(p);

    if (fabs(den) < EPS)
        return (f_r(nnf, i, p, k, eta));
    else
        return ((f_f(nnf, i, j, k + l, eta) - f_f(nnf, i, p, k, eta)) / den);
}

double EvolDF1::f_h(uint nnf, uint i, uint p, uint q, uint j, int k, int l, int m, double eta)
{
    double ll = log(eta), den1 = ai[nnf].at(j) + m - ai[nnf].at(q),
            den2 = ai[nnf].at(q) + l - ai[nnf].at(p),
            den3 = ai[nnf].at(p) + k - ai[nnf].at(i);

    if (fabs(den1) < EPS && fabs(den2) < EPS && fabs(den3) < EPS)
        return (pow(eta, ai[nnf].at(i)) * ll * ll * ll / 6.);
    else if (fabs(den1) < EPS && fabs(den2) < EPS)
        return ((0.5 * pow(eta, ai[nnf].at(p) + k) * ll * ll - f_r(nnf, i, p, k, eta)) / den3);
    else if (fabs(den1) < EPS)
        return ((f_r(nnf, i, q, k + l, eta) - f_g(nnf, i, p, q, k, l, eta)) / den2);
    else
        return ((f_g(nnf, i, p, j, k, l + m, eta) - f_g(nnf, i, p, q, k, l, eta)) / den1);
}

void EvolDF1::CheckNf(indices nm, uint nf) const
{
    if (nm % 10 == 0)
    {
        if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6))
            throw std::runtime_error("EvolDF1::CheckNf(): Wrong number of flavours in anoumalous dimensions");
    } else if (nf != 5)
        throw std::runtime_error("EvolDF1::CheckNf(): Wrong number of flavours in anoumalous dimensions");
}

gslpp::matrix<double> EvolDF1::GammaCC(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    CheckNf(nm, nf);

    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QCD
            // ref.: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
        case 10:
            gammaDF1(0, 0) = -4.;
            gammaDF1(0, 1) = 8. / 3.;
            gammaDF1(1, 0) = 12.;
            break;
        case 20:
            gammaDF1(0, 0) = -145. / 3. + nf * 16. / 9.; // -355./9.
            gammaDF1(0, 1) = -26. + nf * 40. / 27.; // -502./27.
            gammaDF1(1, 0) = -45. + nf * 20. / 3.; // -35./3.
            gammaDF1(1, 1) = -28. / 3.;
            break;
        case 30:
            gammaDF1(0, 0) = -1927. / 2. + nf * 257. / 9. + nf * nf * 40. / 9. + z3 * (224. + nf * 160. / 3.); // -12773./18. + z3*1472./3.
            gammaDF1(0, 1) = 475. / 9. + nf * 362. / 27. - nf * nf * 40. / 27. - z3 * (896. / 3. + nf * 320. / 9.); // 745./9. - z3*4288./9.
            gammaDF1(1, 0) = 307. / 2. + nf * 361. / 3. - nf * nf * 20. / 3. - z3 * (1344. + nf * 160.); // 1177./2. - z3*2144.
            gammaDF1(1, 1) = 1298. / 3. - nf * 76. / 3. - z3 * 224.; // 306. + z3*224.
            break;
            // QED
            // only available for nf = 5
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = -8. / 3.;
            gammaDF1(1, 1) = -8. / 3.;
            break;
        case 11:
            gammaDF1(0, 0) = 169. / 9.;
            gammaDF1(0, 1) = 100. / 27.;
            gammaDF1(1, 0) = 50. / 3.;
            gammaDF1(1, 1) = -8. / 3.;
            break;
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCC(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCP(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    CheckNf(nm, nf);

    gslpp::matrix<double> gammaDF1(2, 4, 0.);
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QCD
            // ref.: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
        case 10:
            gammaDF1(0, 1) = -2. / 9.;
            gammaDF1(1, 1) = 4. / 3.;
            break;
        case 20:
            gammaDF1(0, 0) = -1412. / 243.;
            gammaDF1(0, 1) = -1369. / 243.;
            gammaDF1(0, 2) = 134. / 243.;
            gammaDF1(0, 3) = -35. / 162.;
            gammaDF1(1, 0) = -416. / 81.;
            gammaDF1(1, 1) = 1280. / 81.;
            gammaDF1(1, 2) = 56. / 81.;
            gammaDF1(1, 3) = 35. / 27.;
            break;
        case 30:
            gammaDF1(0, 0) = 269107. / 13122. - nf * 2288. / 729. - z3 * 1360. / 81.; // 63187./13122. - z3*1360./81.
            gammaDF1(0, 1) = -2425817. / 13122. + nf * 30815. / 4374. - z3 * 776. / 81.; // -981796./6561. - z3*776./81.
            gammaDF1(0, 2) = -343783. / 52488. + nf * 392. / 729. + z3 * 124. / 81.; // -202663./52488. + z3*124./81.
            gammaDF1(0, 3) = -37573. / 69984. + nf * 35. / 972. + z3 * 100. / 27.; // -24973./69984. + z3*100./27.
            gammaDF1(1, 0) = 69797. / 2187. + nf * 904. / 243. + z3 * 2720. / 27.; // 110477./2187. + z3*2720./.27.
            gammaDF1(1, 1) = 1457549. / 8748. - nf * 22067. / 729. - z3 * 2768. / 27.; // 133529./8748. - z3*2768./27.
            gammaDF1(1, 2) = -37889. / 8748. - nf * 28. / 243. - z3 * 248. / 27.; // -42929./8748. - z3*248./27.
            gammaDF1(1, 3) = 366919. / 11664. - nf * 35. / 162. - z3 * 110. / 9.; // 354319./11664. - z3*110./9.
            break;
            // QED
            // only available for nf = 5
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            break;
        case 11:
            gammaDF1(0, 1) = 254. / 729.;
            gammaDF1(1, 1) = 1076. / 243.;
            break;
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCP(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCM(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    CheckNf(nm, nf);

    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    double Qu = 2. / 3.;
    double Qd = -1. / 3.;
    double Qbar = n_u * Qu + n_d*Qd;
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QCD
            // ref.: Czakon, Haisch, Misiak, JHEP 0703, 008, hep-ph/0612329
        case 10:
            gammaDF1(0, 0) = 8. / 243. - Qu * 4. / 3.;
            gammaDF1(0, 1) = 173. / 162.;
            gammaDF1(1, 0) = -16. / 81. + Qu * 8.;
            gammaDF1(1, 1) = 70. / 27.;
            break;
        case 20:
            gammaDF1(0, 0) = 12614. / 2187. - nf * 64. / 2187. - Qu * 374. / 27. + nf * Qu * 2. / 27.;
            gammaDF1(0, 1) = 65867. / 5832. + nf * 431. / 5832.;
            gammaDF1(1, 0) = -2332. / 729. + nf * 128. / 729. + Qu * 136. / 9. - nf * Qu * 4. / 9.;
            gammaDF1(1, 1) = 10577. / 486. - nf * 917. / 972.;
            break;
        case 30:
            gammaDF1(0, 0) = 77506102. / 531441. - nf * 875374. / 177147. + nf * nf * 560. / 19683. - Qu * 9731. / 162. +
                    nf * Qu * 11045. / 729. + nf * nf * Qu * 316. / 729. + Qbar * 3695. / 486. + z3 * (-112216. / 6561. + nf * 728. / 729. +
                    Qu * 25508. / 81. - nf * Qu * 64. / 81. - Qbar * 100. / 27.);
            gammaDF1(0, 1) = -421272953. / 1417176. - nf * 8210077. / 472392. - nf * nf * 1955. / 6561. + z3 * (-953042. / 2187. -
                    nf * 10381. / 486.);
            gammaDF1(1, 0) = -15463055. / 177147. + nf * 242204. / 59049. - nf * nf * 1120. / 6561. + Qu * 55748. / 27. -
                    nf * Qu * 33970. / 243. - nf * nf * Qu * 632. / 243. - Qbar * 3695. / 81. + z3 * (365696. / 2187. - nf * 1168. / 243. -
                    Qu * 51232. / 27. - nf * Qu * 1024. / 27. + Qbar * 200. / 9.);
            gammaDF1(1, 1) = 98548513. / 472392. - nf * 5615165. / 78732. - nf * nf * 2489. / 2187. + z3 * (-607103. / 729. -
                    nf * 1679. / 81.);
            break;
            // QED
            // only available for nf = 5
            // ref.: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
        case 01:
            gammaDF1(0, 0) = -832. / 729.;
            gammaDF1(0, 1) = 22. / 243.;
            gammaDF1(1, 0) = -208. / 243.;
            gammaDF1(1, 1) = -116. / 81.;
            break;
        case 11:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCM(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCL(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaCL(): Wrong number of flavours in anoumalous dimensions");


    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = -32. / 27.;
            gammaDF1(1, 0) = -8. / 9.;
            break;
        case 11:
            gammaDF1(0, 0) = -2272. / 729.;
            gammaDF1(1, 0) = 1952. / 243.;
            break;
        case 21:
            gammaDF1(0, 0) = -1359190. / 19683. + z3 * 6976. / 243.;
            gammaDF1(1, 0) = -229696. / 6561. - z3 * 3584. / 81.;
            break;
        case 02:
            gammaDF1(0, 0) = -11680. / 2187.;
            gammaDF1(1, 0) = -2920. / 729.;
            gammaDF1(0, 1) = -416. / 81.;
            gammaDF1(1, 1) = -104. / 27.;
            break;
        case 10:
        case 20:
        case 30:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCL(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCQ(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaCQ(): Wrong number of flavours in anoumalous dimensions");


    gslpp::matrix<double> gammaDF1(2, 4, 0.);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = 32. / 27.;
            gammaDF1(1, 0) = 8. / 9.;
            break;
        case 11:
            gammaDF1(0, 0) = 2272. / 729.;
            gammaDF1(0, 1) = 122. / 81.;
            gammaDF1(0, 3) = 49. / 81.;
            gammaDF1(1, 0) = -1952. / 243.;
            gammaDF1(1, 1) = -748. / 27.;
            gammaDF1(1, 3) = 82. / 27.;
            break;
        case 10:
        case 20:
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCQ(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPP(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    CheckNf(nm, nf);

    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QCD
            // ref.: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
        case 10:
            gammaDF1(0, 1) = -52. / 3.;
            gammaDF1(0, 3) = 2.;
            gammaDF1(1, 0) = -40. / 9.;
            gammaDF1(1, 1) = -160. / 9. + (double) nf * 4. / 3.; // -100./9.
            gammaDF1(1, 2) = 4. / 9.;
            gammaDF1(1, 3) = 5. / 6.;
            gammaDF1(2, 1) = -256. / 3.;
            gammaDF1(2, 3) = 20.;
            gammaDF1(3, 0) = -256. / 9.;
            gammaDF1(3, 1) = -544. / 9. + (double) nf * 40. / 3.; // 56./9.
            gammaDF1(3, 2) = 40. / 9.;
            gammaDF1(3, 3) = -2. / 3.;
            break;
        case 20:
            gammaDF1(0, 0) = -4468. / 81.;
            gammaDF1(0, 1) = -29129. / 81. - nf * 52. / 9.; // -31469./81.
            gammaDF1(0, 2) = 400. / 81.;
            gammaDF1(0, 3) = 3493. / 108. - nf * 2. / 9.; // 3373./108.
            gammaDF1(1, 0) = -13678. / 243. + nf * 368. / 81.; // -8158./243. 
            gammaDF1(1, 1) = -79409. / 243. + nf * 1334. / 81.; // -59399./243.
            gammaDF1(1, 2) = 509. / 486. - nf * 8. / 81.; // 269./486.
            gammaDF1(1, 3) = 13499. / 648. - nf * 5. / 27.; // 12899./648.
            gammaDF1(2, 0) = -244480. / 81. - nf * 160. / 9.; // -251680./81.
            gammaDF1(2, 1) = -29648. / 81. - nf * 2200. / 9.; // -128648./81.
            gammaDF1(2, 2) = 23116. / 81. + nf * 16. / 9.; // 23836./81.
            gammaDF1(2, 3) = 3886. / 27. + nf * 148. / 9.; // 6106./27.
            gammaDF1(3, 0) = 77600. / 243. - nf * 1264. / 81.; // 58640./243.
            gammaDF1(3, 1) = -28808. / 243. + nf * 164. / 81.; // -26348./243.
            gammaDF1(3, 2) = -20324. / 243. + nf * 400. / 81.; // -14324./243.
            gammaDF1(3, 3) = -21211. / 162. + nf * 622. / 27.; // -2551./162.
            break;
        case 30:
            gammaDF1(0, 0) = -4203068. / 2187. + nf * 14012. / 243. - z3 * 608. / 27.; // -3572528./2187. - z3*608./27.
            gammaDF1(0, 1) = -18422762. / 2187. + nf * 888605. / 2916. + nf * nf * 272. / 27. + z3 * (39824. / 27. + nf * 160.); // -58158773./8748. + z3*61424./27.
            gammaDF1(0, 2) = 674281. / 4374. - nf * 1352. / 243. - z3 * 496. / 27.; // 552601./4374. - z3*496./27.
            gammaDF1(0, 3) = 9284531. / 11664. - nf * 2798. / 81. - nf * nf * 26. / 27. - z3 * (1921. / 9. + nf * 20.); // 6989171./11664. - z3*2821./9.
            gammaDF1(1, 0) = -5875184. / 6561. + nf * 217892. / 2187. + nf * nf * 472. / 81. + z3 * (27520. / 81. + nf * 1360. / 9.); // -1651004./6561. + z3*88720./81.
            gammaDF1(1, 1) = -70274587. / 13122. + nf * 8860733. / 17496. - nf * nf * 4010. / 729. + z3 * (16592. / 81. + nf * 2512. / 27.); // -155405353./52488. + z3*54272./81.
            gammaDF1(1, 2) = 2951809. / 52488. - nf * 31175. / 8748. - nf * nf * 52. / 81. - z3 * (3154. / 81. + nf * 136. / 9.); // 1174159./52488. - z3*9274./81.
            gammaDF1(1, 3) = 3227801. / 8748. - nf * 105293. / 11664. - nf * nf * 65. / 54. + z3 * (200. / 27. - nf * 220. / 9.); // 10278809./34992. - z3*3100./27.
            gammaDF1(2, 0) = -194951552. / 2187. + nf * 358672. / 81. - nf * nf * 2144. / 81. + z3 * 87040. / 27.; // -147978032./2187. + z3*87040./27.
            gammaDF1(2, 1) = -130500332. / 2187. - nf * 2949616. / 729. + nf * nf * 3088. / 27. + z3 * (238016. / 27. + nf * 640.); // -168491372./2187. + z3*324416./27.
            gammaDF1(2, 2) = 14732222. / 2187. - nf * 27428. / 81. + nf * nf * 272. / 81. - z3 * 13984. / 27.; // 11213042./2187. - z3*13984./27.
            gammaDF1(2, 3) = 16521659. / 2916. + nf * 8081. / 54. - nf * nf * 316. / 27. - z3 * (22420. / 9. + nf * 200.); // 17850329./2916. - z3*31420./9.
            gammaDF1(3, 0) = 162733912. / 6561. - nf * 2535466. / 2187. + nf * nf * 17920. / 243. + z3 * (174208. / 81. + nf * 12160. / 9.); // 136797922./6561. + z3*721408./81.
            gammaDF1(3, 1) = 13286236. / 6561. - nf * 1826023. / 4374. - nf * nf * 159548. / 729. - z3 * (24832. / 81. + nf * 9440. / 27.); // -72614473./13122. - z3*166432./81.
            gammaDF1(3, 2) = -22191107. / 13122. + nf * 395783. / 4374. - nf * nf * 1720. / 243. - z3 * (33832. / 81. + nf * 1360. / 9.); // -9288181./6561. - z3*95032./81.
            gammaDF1(3, 3) = -32043361. / 8748. + nf * 3353393. / 5832. - nf * nf * 533. / 81. + z3 * (9248. / 27. - nf * 1120. / 9.); // -16664027./17496. - z3*7552./27.
            break;
            // QED
            // only available for nf = 5
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            break;
        case 11:
            gammaDF1(0, 1) = 11116. / 243.;
            gammaDF1(0, 3) = -14. / 3.;
            gammaDF1(1, 0) = 280. / 27.;
            gammaDF1(1, 1) = 18763. / 729.;
            gammaDF1(1, 2) = -28. / 27.;
            gammaDF1(1, 3) = -35. / 18.;
            gammaDF1(2, 1) = 111136. / 243.;
            gammaDF1(2, 3) = -140. / 3.;
            gammaDF1(3, 0) = 2944. / 27.;
            gammaDF1(3, 1) = 193312. / 729.;
            gammaDF1(3, 2) = -280. / 27.;
            gammaDF1(3, 3) = -175. / 9.;
            break;
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPP(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPM(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    CheckNf(nm, nf);

    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    double Qu = 2. / 3.;
    double Qd = -1. / 3.;
    double Qbar = n_u * Qu + n_d*Qd;
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QCD
            // ref.: Czakon, Haisch, Misiak, JHEP 0703, 008, hep-ph/0612329
        case 10:
            gammaDF1(0, 0) = -176. / 81.;
            gammaDF1(0, 1) = 14. / 27.;
            gammaDF1(1, 0) = 88. / 243. - nf * 16. / 81.;
            gammaDF1(1, 1) = 74. / 81. - nf * 49. / 54.;
            gammaDF1(2, 0) = -6272. / 81.;
            gammaDF1(2, 1) = 1736. / 27. + nf * 36.;
            gammaDF1(3, 0) = 3136. / 243. - nf * 160. / 81. + Qbar * 48.;
            gammaDF1(3, 1) = 2372. / 81. + nf * 160. / 27.;
            break;
        case 20:
            gammaDF1(0, 0) = 97876. / 729. - nf * 4352. / 729. - Qbar * 112. / 3.;
            gammaDF1(0, 1) = 42524. / 243. - nf * 2398. / 243.;
            gammaDF1(1, 0) = -70376. / 2187. - nf * 15788. / 2187. + nf * nf * 32. / 729. - Qbar * 140. / 9.;
            gammaDF1(1, 1) = -159718. / 729. - nf * 39719. / 5832. - nf * nf * 253. / 486.;
            gammaDF1(2, 0) = 1764752. / 729. - nf * 65408. / 729. - Qbar * 3136. / 3.;
            gammaDF1(2, 1) = 2281576. / 243. + nf * 140954. / 243. - nf * nf * 14.;
            gammaDF1(3, 0) = 4193840. / 2187. - nf * 324128. / 2187. + nf * nf * 896. / 729. - Qbar * 1136. / 9. - nf * Qbar * 56. / 3.;
            gammaDF1(3, 1) = -3031517. / 729. - nf * 15431. / 1458. - nf * nf * 6031. / 486.;
            break;
        case 30:
            gammaDF1(0, 0) = 102439553. / 177147. - nf * 12273398 / 59049. + nf * nf * 5824. / 6561. + Qbar * 26639. / 81. - nf * Qbar * 8. / 27. +
                    z3 * (3508864. / 2187. - nf * 1904. / 243. - Qbar * 1984. / 9. - nf * Qbar * 64. / 9.);
            gammaDF1(0, 1) = 3205172129. / 472392. - nf * 108963529. / 314928. + nf * nf * 58903. / 4374. + z3 * (-1597588. / 729. +
                    nf * 13028. / 81. - nf * nf * 20. / 9.);
            gammaDF1(1, 0) = -2493414077. / 1062882. - nf * 9901031. / 354294. + nf * nf * 243872. / 59049. - nf * nf * nf * 1184. / 6561. -
                    Qbar * 49993. / 972. + nf * Qbar * 305. / 27. + z3 * (-1922264. / 6561. + nf * 308648. / 2187. - nf * nf * 1280. / 243. +
                    Qbar * 1010. / 9. - nf * Qbar * 200. / 27.);
            gammaDF1(1, 1) = -6678822461. / 2834352. + nf * 127999025. / 1889568. + nf * nf * 1699073. / 157464. + nf * nf * nf * 505. / 4374. +
                    z3 * (2312684. / 2187. + nf * 128347. / 729. + nf * nf * 920. / 81.);
            gammaDF1(2, 0) = 8808397748. / 177147. - nf * 174839456. / 59049. + nf * nf * 1600. / 729. - Qbar * 669694. / 81. + nf * Qbar * 10672. / 27. +
                    z3 * (123543040. / 2187. - nf * 207712. / 243. + nf * nf * 128. / 27. - Qbar * 24880. / 9. - nf * Qbar * 640. / 9.);
            gammaDF1(2, 1) = 29013624461. / 118098. - nf * 64260772. / 19683. - nf * nf * 230962. / 243. - nf * nf * nf * 148. / 27. +
                    z3 * (-69359224. / 729. - nf * 885356. / 81. - nf * nf * 5080. / 9.);
            gammaDF1(3, 0) = 7684242746. / 531441. - nf * 351775414. / 177147. - nf * nf * 479776. / 59049. - nf * nf * nf * 11456. / 6561. +
                    Qbar * 3950201. / 243. - nf * Qbar * 130538. / 81. - nf * nf * Qbar * 592. / 81. + z3 * (7699264. / 6561. + nf * 2854976. / 2187. -
                    nf * nf * 12320. / 243. - Qbar * 108584. / 9. - nf * Qbar * 1136. / 27.);
            gammaDF1(3, 1) = -72810260309. / 708588. + nf * 2545824851. / 472392. - nf * nf * 33778271. / 78732. - nf * nf * nf * 3988. / 2187. +
                    z3 * (-61384768. / 2187. - nf * 685472. / 729. + nf * nf * 350. / 81.);
            break;
            // QED    
            // only available for nf = 5
            // ref.: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
        case 01:
            gammaDF1(0, 0) = -20. / 243.;
            gammaDF1(0, 1) = 20. / 81.;
            gammaDF1(1, 0) = -176. / 729.;
            gammaDF1(1, 1) = 14. / 243.;
            gammaDF1(2, 0) = -22712. / 243.;
            gammaDF1(2, 1) = 1328. / 81.;
            gammaDF1(3, 0) = -6272. / 729.;
            gammaDF1(3, 1) = -1180. / 243.;
            break;
        case 11:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPM(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPL(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaPL(): Wrong number of flavours in anoumalous dimensions");


    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066

        case 01:
            gammaDF1(0, 0) = -16. / 9.;
            gammaDF1(1, 0) = 32. / 27.;
            gammaDF1(2, 0) = -112. / 9.;
            gammaDF1(3, 0) = 512. / 27.;
            break;
        case 11:
            gammaDF1(0, 0) = -6752. / 243.;
            gammaDF1(1, 0) = -2192. / 729.;
            gammaDF1(2, 0) = -84032. / 243.;
            gammaDF1(3, 0) = -37856. / 729.;
            break;
        case 21:
            gammaDF1(0, 0) = -1290092. / 6561. + z3 * 3200. / 81.;
            gammaDF1(1, 0) = -819971. / 19683. - z3 * 19936. / 243.;
            gammaDF1(2, 0) = -16821944. / 6561. + z3 * 30464. / 81.;
            gammaDF1(3, 0) = -17787368. / 19683. - z3 * 286720. / 243.;
            break;
        case 02:
            gammaDF1(0, 0) = -39752. / 729.;
            gammaDF1(1, 0) = 1024. / 2187.;
            gammaDF1(2, 0) = -381344. / 729.;
            gammaDF1(3, 0) = 24832. / 2187.;
            gammaDF1(0, 1) = -136. / 27.;
            gammaDF1(1, 1) = -448. / 81.;
            gammaDF1(2, 1) = -15616. / 27.;
            gammaDF1(3, 1) = -7936. / 81.;
            break;
        case 10:
        case 20:
        case 30:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPL(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPQ(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaPQ(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(4, 4, 0.);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = 76. / 9.;
            gammaDF1(0, 2) = -2. / 3.;
            gammaDF1(1, 0) = -32. / 27.;
            gammaDF1(1, 1) = 20. / 3.;
            gammaDF1(1, 3) = -2. / 3.;
            gammaDF1(2, 0) = 496. / 9.;
            gammaDF1(2, 2) = -20. / 3.;
            gammaDF1(3, 0) = -512. / 27.;
            gammaDF1(3, 1) = 128. / 3.;
            gammaDF1(3, 3) = -20. / 3.;
            break;
        case 11:
            gammaDF1(0, 0) = -23488. / 243.;
            gammaDF1(0, 1) = 6280. / 27.;
            gammaDF1(0, 2) = 112. / 9.;
            gammaDF1(0, 3) = -538. / 27.;
            gammaDF1(1, 0) = 31568. / 729.;
            gammaDF1(1, 1) = 9481. / 81.;
            gammaDF1(1, 2) = -92. / 27.;
            gammaDF1(1, 3) = -1012. / 81.;
            gammaDF1(2, 0) = -233920. / 243.;
            gammaDF1(2, 1) = 68848. / 27.;
            gammaDF1(2, 2) = 1120. / 9.;
            gammaDF1(2, 3) = -5056. / 27.;
            gammaDF1(3, 0) = 352352. / 729.;
            gammaDF1(3, 1) = 116680. / 81.;
            gammaDF1(3, 2) = -752. / 27.;
            gammaDF1(3, 3) = -10147. / 81.;
            break;
        case 10:
        case 20:
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPQ(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaMM(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    CheckNf(nm, nf);

    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    double Qu = 2. / 3.;
    double Qd = -1. / 3.;
    double Qbar = n_u * Qu + n_d*Qd;
    double z3 = gslpp_special_functions::zeta(3);

    switch (nm)
    {
            // QCD
            //ref.: Gorbahn, Haisch, Misiak, Phys. Rev. Lett. 95, 102004, hep-ph/0504194
        case 10:
            gammaDF1(0, 0) = 32. / 3.;
            gammaDF1(1, 0) = Qd * 32. / 3.;
            gammaDF1(1, 1) = 28. / 3.;
            break;
        case 20:
            gammaDF1(0, 0) = 1936. / 9. - nf * 224. / 27.;
            gammaDF1(1, 0) = Qd * 368. / 3. - nf * Qd * 224. / 27.;
            gammaDF1(1, 1) = 1456. / 9. - nf * 61. / 27.;
            break;
        case 30:
            gammaDF1(0, 0) = 307448. / 81. - nf * 23776. / 81. - nf * nf * 352. / 81. + z3 * (-1856. / 27. - nf * 1280. / 9.);
            gammaDF1(1, 0) = -Qbar * 1600. / 27. + Qd * 159872. / 81. - nf * Qd * 17108. / 81. - nf * nf * Qd * 352. / 81. + z3 * (Qbar * 640. / 9. -
                    Qd * 1856. / 27. - nf * Qd * 1280. / 9.);
            gammaDF1(1, 1) = 268807. / 81. - nf * 4343. / 27. - nf * nf * 461. / 81. + z3 * (-28624. / 27. - nf * 1312. / 9.);
            break;
            // QED
            // only available for nf = 5
            // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
        case 01:
            gammaDF1(0, 0) = 16. / 9.;
            gammaDF1(0, 1) = -8. / 3.;
            gammaDF1(1, 1) = 8. / 9.;
            break;
        case 11:
            gammaDF1(0, 0) = -256. / 27.;
            gammaDF1(0, 1) = -52. / 9.;
            gammaDF1(1, 0) = 128. / 81.;
            gammaDF1(1, 1) = -40. / 27.;
            break;
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaMM(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaLL(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaLL(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(2, 2, 0.);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = 8.;
            gammaDF1(0, 1) = -4.;
            gammaDF1(1, 0) = -4.;
            break;
        case 11:
            gammaDF1(0, 1) = 16.;
            gammaDF1(1, 0) = 16.;
            break;
        case 10:
        case 20:
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaLL(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQP(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaQP(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(4, 4, 0.);

    switch (nm)
    {
            // QCD
            //ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 10:
            gammaDF1(0, 1) = -8. / 9.;
            gammaDF1(1, 1) = 16. / 27.;
            gammaDF1(2, 1) = -128. / 9.;
            gammaDF1(3, 1) = 184. / 27.;
            break;
        case 20:
            gammaDF1(0, 0) = 832. / 243.;
            gammaDF1(0, 1) = -4000. / 243.;
            gammaDF1(0, 2) = -112. / 243.;
            gammaDF1(0, 3) = -70. / 81.;
            gammaDF1(1, 0) = 3376. / 729.;
            gammaDF1(1, 1) = 6344. / 729.;
            gammaDF1(1, 2) = -280. / 729.;
            gammaDF1(1, 3) = 55. / 486.;
            gammaDF1(2, 0) = 2272. / 243.;
            gammaDF1(2, 1) = -72088. / 243.;
            gammaDF1(2, 2) = -688. / 243.;
            gammaDF1(2, 3) = -1240. / 81.;
            gammaDF1(3, 0) = 45424. / 729.;
            gammaDF1(3, 1) = 84236. / 729.;
            gammaDF1(3, 2) = -3880. / 729.;
            gammaDF1(3, 3) = 1220. / 243.;
            break;
            // QED
            // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
        case 01:
            gammaDF1(0, 0) = 40. / 27.;
            gammaDF1(0, 2) = -4. / 27.;
            gammaDF1(1, 1) = 40. / 27.;
            gammaDF1(1, 3) = -4. / 27.;
            gammaDF1(2, 0) = 256. / 27.;
            gammaDF1(2, 2) = -40. / 27.;
            gammaDF1(3, 1) = 256. / 27.;
            gammaDF1(3, 3) = -40. / 27.;
            break;
        case 11:
            gammaDF1(0, 0) = -2240. / 81.;
            gammaDF1(0, 1) = 39392. / 729.;
            gammaDF1(0, 2) = 224. / 81.;
            gammaDF1(0, 3) = -92. / 27.;
            gammaDF1(1, 0) = 2176. / 243.;
            gammaDF1(1, 1) = 84890. / 2187.;
            gammaDF1(1, 2) = -184. / 243.;
            gammaDF1(1, 3) = -224. / 81.;
            gammaDF1(2, 0) = -23552. / 81.;
            gammaDF1(2, 1) = 399776. / 729.;
            gammaDF1(2, 2) = 2240. / 81.;
            gammaDF1(2, 3) = -752. / 27.;
            gammaDF1(3, 0) = 23296. / 243.;
            gammaDF1(3, 1) = 933776. / 2187.;
            gammaDF1(3, 2) = -1504. / 243.;
            gammaDF1(3, 3) = -2030. / 81.;
            break;
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQP(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQM(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaQM(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(4, 2, 0.);

    switch (nm)
    {
            // QCD
            // ref.: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
        case 10:
            gammaDF1(0, 0) = 176. / 243.;
            gammaDF1(0, 1) = -14. / 81.;
            gammaDF1(1, 0) = -136. / 729.;
            gammaDF1(1, 1) = -295. / 486.;
            gammaDF1(2, 0) = 6272. / 243.;
            gammaDF1(2, 1) = -764. / 81.;
            gammaDF1(3, 0) = 39152. / 729.;
            gammaDF1(3, 1) = -1892. / 243.;
            break;
        case 20:
        case 30:
        case 01:
        case 11:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQM(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQL(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaQL(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(4, 2, 0.);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = -272. / 27.;
            gammaDF1(1, 0) = -32. / 81.;
            gammaDF1(2, 0) = -2768. / 27.;
            gammaDF1(3, 0) = -512. / 81.;
            break;
        case 11:
            gammaDF1(0, 0) = -24352. / 729.;
            gammaDF1(1, 0) = 54608. / 2187.;
            gammaDF1(2, 0) = -227008. / 729.;
            gammaDF1(3, 0) = 551648. / 2187.;
            break;
        case 10:
        case 20:
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQL(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQQ(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaQQ(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(4, 4, 0.);

    switch (nm)
    {
            // QCD
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 10:
            gammaDF1(0, 1) = -20.;
            gammaDF1(0, 3) = 2.;
            gammaDF1(1, 0) = -40. / 9.;
            gammaDF1(1, 1) = -52. / 3.;
            gammaDF1(1, 2) = 4. / 9.;
            gammaDF1(1, 3) = 5. / 6.;
            gammaDF1(2, 1) = -128.;
            gammaDF1(2, 3) = 20.;
            gammaDF1(3, 0) = -256. / 9.;
            gammaDF1(3, 1) = -160. / 3.;
            gammaDF1(3, 2) = 40. / 9.;
            gammaDF1(3, 3) = -2. / 3.;
            break;
        case 20:
            gammaDF1(0, 0) = -404. / 9.;
            gammaDF1(0, 1) = -3077. / 9.;
            gammaDF1(0, 2) = 32. / 9.;
            gammaDF1(0, 3) = 1031. / 36.;
            gammaDF1(1, 0) = -2698. / 81.;
            gammaDF1(1, 1) = -8035. / 27.;
            gammaDF1(1, 2) = -49. / 162.;
            gammaDF1(1, 3) = 4493. / 216.;
            gammaDF1(2, 0) = -19072. / 9.;
            gammaDF1(2, 1) = -14096. / 9.;
            gammaDF1(2, 2) = 1708. / 9.;
            gammaDF1(2, 3) = 1622. / 9.;
            gammaDF1(3, 0) = 32288. / 81.;
            gammaDF1(3, 1) = -15976. / 27.;
            gammaDF1(3, 2) = -6692. / 81.;
            gammaDF1(3, 3) = -2437. / 54.;
            break;
            // QED
            // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
        case 01:
            gammaDF1(0, 0) = 332. / 27.;
            gammaDF1(0, 2) = -2. / 9.;
            gammaDF1(1, 0) = 32. / 81.;
            gammaDF1(1, 1) = 20. / 9.;
            gammaDF1(1, 3) = -2. / 9.;
            gammaDF1(2, 0) = 3152. / 27.;
            gammaDF1(2, 2) = -20. / 9.;
            gammaDF1(3, 0) = 512. / 81.;
            gammaDF1(3, 1) = 128. / 9.;
            gammaDF1(3, 3) = -20. / 9.;
            break;
        case 11:
            gammaDF1(0, 0) = -5888. / 729.;
            gammaDF1(0, 1) = 13916. / 81.;
            gammaDF1(0, 2) = 112. / 27.;
            gammaDF1(0, 3) = -812. / 81.;
            gammaDF1(1, 0) = -2552. / 2187.;
            gammaDF1(1, 1) = 15638. / 243.;
            gammaDF1(1, 2) = -176. / 81.;
            gammaDF1(1, 3) = -2881. / 486.;
            gammaDF1(2, 0) = -90944. / 729.;
            gammaDF1(2, 1) = 90128. / 81.;
            gammaDF1(2, 2) = 1120. / 27.;
            gammaDF1(2, 3) = -1748. / 81.;
            gammaDF1(3, 0) = 1312. / 2187.;
            gammaDF1(3, 1) = 102488. / 243.;
            gammaDF1(3, 2) = -1592. / 81.;
            gammaDF1(3, 3) = -6008. / 243.;
            break;
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQQ(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBP(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaBP(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(1, 4, 0.);

    switch (nm)
    {
            // QCD
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 10:
            gammaDF1(0, 1) = 4. / 3.;
            break;
        case 20:
            gammaDF1(0, 0) = -1576. / 81.;
            gammaDF1(0, 1) = 446. / 27.;
            gammaDF1(0, 2) = 172. / 81.;
            gammaDF1(0, 3) = 40. / 27.;
            break;
            // QED
            // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
        case 01:
            break;
        case 11:
            gammaDF1(0, 1) = -232. / 81.;
            break;
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBP(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBL(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaBL(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(1, 2, 0.);

    switch (nm)
    {
            // QED
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 01:
            gammaDF1(0, 0) = 16. / 9.;
            break;
        case 11:
            gammaDF1(0, 0) = -8. / 9.;
            break;
        case 10:
        case 20:
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBL(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBQ(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaBQ(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(1, 4, 0.);

    switch (nm)
    {
            // QED
            // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
        case 01:
            gammaDF1(0, 0) = -16. / 9.;
            break;
        case 11:
            gammaDF1(0, 1) = 580. / 27.;
            gammaDF1(0, 3) = -94. / 27.;
            break;
        case 10:
        case 20:
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBQ(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBB(indices nm, uint n_u, uint n_d) const
{
    uint nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (nf != 5)
        throw std::runtime_error("EvolDF1::GammaBB(): Wrong number of flavours in anoumalous dimensions");

    gslpp::matrix<double> gammaDF1(1, 1, 0.);

    switch (nm)
    {
            // QCD
            // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
        case 10:
            gammaDF1(0, 0) = 4.;
            break;
        case 20:
            gammaDF1(0, 0) = 325. / 9.;
            break;
            // QED
            // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
        case 01:
            gammaDF1(0, 0) = 4. / 3.;
            break;
        case 11:
            gammaDF1(0, 0) = -388. / 9.;
            break;
        case 30:
        case 21:
        case 02:
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBB(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::AnomalousDimension(indices nm, uint n_u, uint n_d) const
{
    gslpp::matrix<double> gammaDF1(nops, nops, 0.);

    // assign blocks according to user request: "C", "CP", "CPM", "L", "CPML", "CPQB", "CPMQB", "CPMLQB"

    if (blocks.compare("C") == 0)
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
    else if (blocks.compare("CP") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
    } else if (blocks.compare("CPM") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCM(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPM(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaMM(nm, n_u, n_d));
    } else if (blocks.compare("CPQ") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCQ(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPQ(nm, n_u, n_d));
        gammaDF1.assign(6, 2, GammaQP(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaQQ(nm, n_u, n_d));
    } else if (blocks.compare("L") == 0)
    {
        gammaDF1.assign(0, 0, GammaLL(nm, n_u, n_d));
    } else if (blocks.compare("CPL") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCL(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPL(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaLL(nm, n_u, n_d));
    } else if (blocks.compare("CPML") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCM(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPM(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaMM(nm, n_u, n_d));
        gammaDF1.assign(0, 8, GammaCL(nm, n_u, n_d));
        gammaDF1.assign(2, 8, GammaPL(nm, n_u, n_d));
        gammaDF1.assign(8, 8, GammaLL(nm, n_u, n_d));
    } else if (blocks.compare("CPQB") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCQ(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPQ(nm, n_u, n_d));
        gammaDF1.assign(6, 2, GammaQP(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaQQ(nm, n_u, n_d));
        gammaDF1.assign(10, 2, GammaBP(nm, n_u, n_d));
        gammaDF1.assign(10, 6, GammaBQ(nm, n_u, n_d));
        gammaDF1.assign(10, 10, GammaBB(nm, n_u, n_d));
    } else if (blocks.compare("CPMQB") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCM(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPM(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaMM(nm, n_u, n_d));
        gammaDF1.assign(0, 8, GammaCQ(nm, n_u, n_d));
        gammaDF1.assign(2, 8, GammaPQ(nm, n_u, n_d));
        gammaDF1.assign(8, 2, GammaQP(nm, n_u, n_d));
        gammaDF1.assign(8, 6, GammaQM(nm, n_u, n_d));
        gammaDF1.assign(8, 8, GammaQQ(nm, n_u, n_d));
        gammaDF1.assign(12, 2, GammaBP(nm, n_u, n_d));
        gammaDF1.assign(12, 8, GammaBQ(nm, n_u, n_d));
        gammaDF1.assign(12, 12, GammaBB(nm, n_u, n_d)); // *** does BM exists?
    } else if (blocks.compare("CPMLQB") == 0)
    {
        gammaDF1.assign(0, 0, GammaCC(nm, n_u, n_d));
        gammaDF1.assign(0, 2, GammaCP(nm, n_u, n_d));
        gammaDF1.assign(2, 2, GammaPP(nm, n_u, n_d));
        gammaDF1.assign(0, 6, GammaCM(nm, n_u, n_d));
        gammaDF1.assign(2, 6, GammaPM(nm, n_u, n_d));
        gammaDF1.assign(6, 6, GammaMM(nm, n_u, n_d));
        gammaDF1.assign(0, 8, GammaCL(nm, n_u, n_d));
        gammaDF1.assign(2, 8, GammaPL(nm, n_u, n_d));
        gammaDF1.assign(8, 8, GammaLL(nm, n_u, n_d));

        gammaDF1.assign(0, 10, GammaCQ(nm, n_u, n_d));
        gammaDF1.assign(2, 10, GammaPQ(nm, n_u, n_d));
        gammaDF1.assign(10, 2, GammaQP(nm, n_u, n_d));
        gammaDF1.assign(10, 6, GammaQM(nm, n_u, n_d));
        gammaDF1.assign(10, 8, GammaQL(nm, n_u, n_d));
        gammaDF1.assign(10, 10, GammaQQ(nm, n_u, n_d));
        gammaDF1.assign(14, 2, GammaBP(nm, n_u, n_d));
        gammaDF1.assign(14, 8, GammaBL(nm, n_u, n_d));
        gammaDF1.assign(14, 10, GammaBQ(nm, n_u, n_d));
        gammaDF1.assign(14, 14, GammaBB(nm, n_u, n_d)); // *** does BM exists?
    } else
        throw std::runtime_error("EvolDF1::AnomalousDimension(): block not implemented");

    return (gammaDF1);
}

const Expanded<gslpp::matrix<double> >& EvolDF1::DF1Evol(double mu, double M, schemes scheme)
{
    if (nfmin == 5 && nfmax == 5 && (model.Nf(mu) != 5. || model.Nf(M) != 5.))
        throw std::runtime_error("EvolDF1::Df1Evol(): only nf = 5 available.");

    switch (scheme)
    {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1::Df1Evol(): scheme " + out.str() + " not implemented ");
    }

    double alsM = model.getAlsM();
    double MAls = model.getMAls();
    if (alsM == alsM_cache && MAls == MAls_cache)
    {
        if (mu == this->mu && M == this->M && scheme == this->scheme)
            return (getEvol());
    }
    alsM_cache = alsM;
    MAls_cache = MAls;

    if (M < mu)
    {
        std::stringstream out;
        out << "M = " << M << " < mu = " << mu;
        throw std::runtime_error("EvolDF1::Df1Evol(): " + out.str() + ".");
    }

    setScales(mu, M); // also assign evol to identity

    double m_down = mu;
    double m_up = model.AboveTh(m_down);
    double nf = model.Nf(m_down);

    while (m_up < M)
    { // where are the nf thresholds? ???? <<<<<<<<<<
        DF1Ev(m_down, m_up, (int) nf, scheme);
        m_down = m_up;
        m_up = model.AboveTh(m_down);
        nf += 1.;
    }
    DF1Ev(m_down, M, (int) nf, scheme);

    return (getEvol());
}

//gslpp::matrix<double>& EvolDF1::DF1Evol(double mu, double M, orders_qed ord, schemes scheme)
//{
//    if(ord > order_qed)
//        throw std::runtime_error("EvolDF1::Df1Evol(): order not present in this Hamiltonian.");
//    double MAls = model.getMAls();
//    if(model.Nf(mu) != 5. || model.Nf(M) != 5. || model.Nf(MAls) != 5.)
//        throw std::runtime_error("EvolDF1::Df1Evol(): only nf = 5 available.");
//        
//    switch (scheme) {
//        case NDR:
//            break;
//        case LRI:
//        case HV:
//        default:
//            std::stringstream out;
//            out << scheme;
//            throw std::runtime_error("EvolDF1::Df1Evol(): scheme " + out.str() + " not implemented "); 
//    }
//    
//// WRONG caching! Rethink *************
//    
//    double alsM = model.getAlsM();
//    if(alsM == alsM_cache && MAls == MAls_cache) {
//        if (mu == this->mu && M == this->M && scheme == this->scheme)
//            return (*Evol(ord));        
//    }
//    alsM_cache = alsM;
//    MAls_cache = MAls;
//        
//    if (M < mu) {
//        std::stringstream out;
//        out << "M = " << M << " < mu = " << mu;
//        throw std::runtime_error("EvolDF1::Df1Evol(): " + out.str() + ".");
//    }
//
//    setScales(mu, M); // also assign evol to identity
//
//    DF1Ev(mu, M, 5, scheme); // only 5 flavour *******************
//
//    return (*Evol(ord));
//}

void EvolDF1::DF1Ev(double mu, double M, int nf, schemes scheme)
{
    gslpp::matrix<double> mtmp(nops, 0.);
    std::vector<std::vector<gslpp::matrix<double> > > vtmp2;
    std::vector<gslpp::matrix<double> > vtmp;
    for (int j = 0; j <= order_qed; j++)
        vtmp.push_back(mtmp);
    for (int i = 0; i <= order_qcd; i++)
        vtmp2.push_back(vtmp);
    Expanded<gslpp::matrix<double> > res(vtmp2);
    //     gslpp::matrix<double> res01(nops, 0.), res02(nops, 0.), res11(nops, 0.), res12(nops, 0.),
    //            res21(nops, 0.), resLO(nops, 0.), resNLO(nops, 0.), resNNLO(nops, 0.);

    //    uint a, b, i, j, p, q
    uint nnf = nf - nfmin;
    double b0, b0e, b5, alsM, eta, omega, lambda; //, term;
    std::map< std::vector<uint>, double >::iterator itr;
    //    std::vector<uint> v;


    //    alsM = model.Als(M) / 4. / M_PI;
    //    double alsmu = model.Als(mu) / 4. / M_PI;
    b0 = model.Beta_s(00, nf);
    alsM = model.Als(M, FULLNNNLO, true, order_qed == QED0 ? false : true);
    eta = alsM / model.Als(mu, FULLNNNLO, true, order_qed == QED0 ? false : true);
    //    eta = alsM / model.Als(mu);
    omega = 2. * b0 * alsM / 4. / M_PI;

    for (itr = vM0vi[nnf].begin(); itr != vM0vi[nnf].end(); ++itr)
    {
        const std::vector<uint> &v = itr->first;
        const uint &a = v[0];
        const uint &b = v[1];
        const uint &i = v[2];

        res.setMatrixElement(QCD0, QED0, a, b, res.getOrd(QCD0, QED0)(a, b) + itr->second * pow(eta, ai[nnf].at(i)));
    }

    for (itr = vM1vi[nnf].begin(); itr != vM1vi[nnf].end(); ++itr)
    {
        const std::vector<uint> &v = itr->first;
        const uint &a = v[0];
        const uint &b = v[1];
        const uint &i = v[2];
        const uint &j = v[3];

        res.setMatrixElement(QCD1, QED0, a, b, res.getOrd(QCD1, QED0)(a, b) + omega * itr->second * f_f(nnf, i, j, -1, eta));
    }

    for (itr = vM2vi[nnf].begin(); itr != vM2vi[nnf].end(); ++itr)
    {
        const std::vector<uint> &v = itr->first;
        const uint &a = v[0];
        const uint &b = v[1];
        const uint &i = v[2];
        const uint &j = v[3];

        res.setMatrixElement(QCD2, QED0, a, b, res.getOrd(QCD2, QED0)(a, b) + omega * omega * itr->second * f_f(nnf, i, j, -2, eta));
    }

    for (itr = vM11vi[nnf].begin(); itr != vM11vi[nnf].end(); ++itr)
    {
        const std::vector<uint> &v = itr->first;
        const uint &a = v[0];
        const uint &b = v[1];
        const uint &i = v[2];
        const uint &j = v[3];
        const uint &p = v[4];

        res.setMatrixElement(QCD2, QED0, a, b, res.getOrd(QCD2, QED0)(a, b) + omega * omega * itr->second * f_g(nnf, i, p, j, -1, -1, eta));
    }

    if (order_qed != QED0)
    {
        b0e = model.Beta_e(00, nf);
        b5 = model.Beta_e(01, nf) / 2. / b0 / b0e - model.Beta_s(10, nf) / 2. / b0 / b0;
        lambda = b0e * model.Ale(M, FULLNLO) / b0 / model.Als(M, FULLNNNLO, true, true); // WARNING: CHANGE ME!!!

        for (itr = vM3vi[nnf].begin(); itr != vM3vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const double &term = itr->second;

            res.setMatrixElement(QCD0, QED1, a, b, res.getOrd(QCD0, QED1)(a, b) + lambda * term * f_f(nnf, i, j, 1, eta));
            res.setMatrixElement(QCD0, QED2, a, b, res.getOrd(QCD0, QED2)(a, b) + lambda * lambda * term * (f_f(nnf, i, j, 2, eta) - f_f(nnf, i, j, 1, eta)));
            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * b5 * term * f_r(nnf, i, j, 1, eta));
        }

        for (itr = vM4vi[nnf].begin(); itr != vM4vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const double &term = itr->second;

            res.setMatrixElement(QCD1, QED1, a, b, res.getOrd(QCD1, QED1)(a, b) + omega * lambda * term * f_f(nnf, i, j, 0, eta));
            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) - omega * lambda * lambda * term * f_f(nnf, i, j, 0, eta));
        }

        for (itr = vM5vi[nnf].begin(); itr != vM5vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_f(nnf, i, j, -1, eta));
        }

        for (itr = vM6vi[nnf].begin(); itr != vM6vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];

            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * itr->second * f_f(nnf, i, j, 1, eta));
        }

        for (itr = vM33vi[nnf].begin(); itr != vM33vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD0, QED2, a, b, res.getOrd(QCD0, QED2)(a, b) + lambda * lambda * itr->second * f_g(nnf, i, p, j, 1, 1, eta));
        }

        for (itr = vM13vi[nnf].begin(); itr != vM13vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const double &term = itr->second;

            res.setMatrixElement(QCD1, QED1, a, b, res.getOrd(QCD1, QED1)(a, b) + omega * lambda * term * f_g(nnf, i, p, j, -1, 1, eta));
            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * term * (f_g(nnf, i, p, j, -1, 2, eta) - f_g(nnf, i, p, j, -1, 1, eta)));
        }

        for (itr = vM31vi[nnf].begin(); itr != vM31vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const double &term = itr->second;

            res.setMatrixElement(QCD1, QED1, a, b, res.getOrd(QCD1, QED1)(a, b) + omega * lambda * term * f_g(nnf, i, p, j, 1, -1, eta));
            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * term * (f_g(nnf, i, p, j, 2, -1, eta) - f_g(nnf, i, p, j, 1, -1, eta)));
        }

        for (itr = vM34vi[nnf].begin(); itr != vM34vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * itr->second * f_g(nnf, i, p, j, 1, 0, eta));
        }

        for (itr = vM43vi[nnf].begin(); itr != vM43vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * itr->second * f_g(nnf, i, p, j, 0, 1, eta));
        }

        for (itr = vM23vi[nnf].begin(); itr != vM23vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_g(nnf, i, p, j, -2, 1, eta));
        }

        for (itr = vM32vi[nnf].begin(); itr != vM32vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_g(nnf, i, p, j, 1, -2, eta));
        }

        for (itr = vM14vi[nnf].begin(); itr != vM14vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_g(nnf, i, p, j, -1, 0, eta));
        }

        for (itr = vM41vi[nnf].begin(); itr != vM41vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_g(nnf, i, p, j, 0, -1, eta));
        }

        for (itr = vM113vi[nnf].begin(); itr != vM113vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const uint &q = v[5];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_h(nnf, i, p, q, j, -1, -1, 1, eta));
        }

        for (itr = vM131vi[nnf].begin(); itr != vM131vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const uint &q = v[5];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_h(nnf, i, p, q, j, -1, 1, -1, eta));
        }

        for (itr = vM311vi[nnf].begin(); itr != vM311vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const uint &q = v[5];

            res.setMatrixElement(QCD2, QED1, a, b, res.getOrd(QCD2, QED1)(a, b) + omega * omega * lambda * itr->second * f_h(nnf, i, p, q, j, 1, -1, -1, eta));
        }

        for (itr = vM133vi[nnf].begin(); itr != vM133vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const uint &q = v[5];

            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * itr->second * f_h(nnf, i, p, q, j, -1, 1, 1, eta));
        }

        for (itr = vM313vi[nnf].begin(); itr != vM313vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const uint &q = v[5];

            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * itr->second * f_h(nnf, i, p, q, j, 1, -1, 1, eta));
        }

        for (itr = vM331vi[nnf].begin(); itr != vM331vi[nnf].end(); ++itr)
        {
            const std::vector<uint> &v = itr->first;
            const uint &a = v[0];
            const uint &b = v[1];
            const uint &i = v[2];
            const uint &j = v[3];
            const uint &p = v[4];
            const uint &q = v[5];

            res.setMatrixElement(QCD1, QED2, a, b, res.getOrd(QCD1, QED2)(a, b) + omega * lambda * lambda * itr->second * f_h(nnf, i, p, q, j, 1, 1, -1, eta));
        }
        /*
        for (a = 0; a < nops; a++)
            for (b = 0; b < nops; b++)
            {
                for (i = 0; i < nops; i++)
                    for (j = 0; j < nops; j++)
                    {
                        res01(a, b) += (lambda * vM3vi.at(idx(nf, a, b, i, j)) * f_f(nf, i, j, 1, eta)).real();
                        res02(a, b) += (lambda * lambda * vM3vi.at(idx(nf, a, b, i, j)) * (f_f(nf, i, j, 2, eta) - f_f(nf, i, j, 1, eta))).real();
                        res11(a, b) += (omega * lambda * vM4vi.at(idx(nf, a, b, i, j)) * f_f(nf, i, j, 0, eta)).real();
                        res12(a, b) += (omega * lambda * lambda * (-vM4vi.at(idx(nf, a, b, i, j)) * f_f(nf, i, j, 0, eta) + vM6vi.at(idx(nf, a, b, i, j)) * f_f(nf, i, j, 1, eta)) +
                                b5 * vM3vi.at(idx(nf, a, b, i, j)) * f_r(nf, i, j, 1, eta)).real();
                        res21(a, b) += (omega * omega * lambda * vM5vi.at(idx(nf, a, b, i, j)) * f_f(nf, i, j, -1, eta)).real();
                        for (p = 0; p < nops; p++)
                        {
                            res02(a, b) += (lambda * lambda * vM33vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, 1, 1, eta)).real();
                            res11(a, b) += (omega * lambda * (vM13vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, -1, 1, eta) + vM31vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, 1, -1, eta))).real();
                            res12(a, b) += (omega * lambda * lambda * (vM13vi.at(idx(nf, a, b, i, j, p)) * (f_g(nf, i, p, j, -1, 2, eta) - f_g(nf, i, p, j, -1, 1, eta)) +
                                    vM31vi.at(idx(nf, a, b, i, j, p)) * (f_g(nf, i, p, j, 2, -1, eta) - f_g(nf, i, p, j, 1, -1, eta)) + vM34vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, 1, 0, eta) +
                                    vM43vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, 0, 1, eta))).real();
                            res21(a, b) += (omega * omega * lambda * (vM14vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, -1, 0, eta) + vM41vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, 0, -1, eta) +
                                    vM23vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, -2, 1, eta) + vM32vi.at(idx(nf, a, b, i, j, p)) * f_g(nf, i, p, j, 1, -2, eta))).real();

                            for (q = 0; q < nops; q++)
                            {
                                res12(a, b) += (omega * lambda * lambda * (vM133vi.at(idx(nf, a, b, i, j, p, q)) * f_h(nf, i, p, q, j, -1, 1, 1, eta) + vM313vi.at(idx(nf, a, b, i, j, p, q)) *
                                        f_h(nf, i, p, q, j, 1, -1, 1, eta) + vM331vi.at(idx(nf, a, b, i, j, p, q)) * f_h(nf, i, p, q, j, 1, 1, -1, eta))).real();
                                res21(a, b) += (omega * omega * lambda * (vM113vi.at(idx(nf, a, b, i, j, p, q)) * f_h(nf, i, p, q, j, -1, -1, 1, eta) + vM131vi.at(idx(nf, a, b, i, j, p, q)) *
                                        f_h(nf, i, p, q, j, -1, 1, -1, eta) + vM311vi.at(idx(nf, a, b, i, j, p, q)) * f_h(nf, i, p, q, j, 1, -1, -1, eta))).real();
                            }
                        }
                    }
                if (fabs(res01(a, b)) < EPS) res01(a, b) = 0.;
                if (fabs(res02(a, b)) < EPS) res02(a, b) = 0.;
                if (fabs(res11(a, b)) < EPS) res11(a, b) = 0.;
                if (fabs(res12(a, b)) < EPS) res12(a, b) = 0.;
                if (fabs(res21(a, b)) < EPS) res21(a, b) = 0.;
            }
         */
        for (uint a = 0; a < nops; a++)
            for (uint b = 0; b < nops; b++)
            {
                if (fabs(res.getOrd(QCD0, QED0)(a, b)) < EPS) res.setMatrixElement(QCD0, QED0, a, b, 0.);
                if (fabs(res.getOrd(QCD1, QED0)(a, b)) < EPS) res.setMatrixElement(QCD1, QED0, a, b, 0.);
                if (fabs(res.getOrd(QCD2, QED0)(a, b)) < EPS) res.setMatrixElement(QCD2, QED0, a, b, 0.);
                if (fabs(res.getOrd(QCD0, QED1)(a, b)) < EPS) res.setMatrixElement(QCD0, QED1, a, b, 0.);
                if (fabs(res.getOrd(QCD0, QED2)(a, b)) < EPS) res.setMatrixElement(QCD0, QED2, a, b, 0.);
                if (fabs(res.getOrd(QCD1, QED1)(a, b)) < EPS) res.setMatrixElement(QCD1, QED1, a, b, 0.);
                if (fabs(res.getOrd(QCD1, QED2)(a, b)) < EPS) res.setMatrixElement(QCD1, QED2, a, b, 0.);
                if (fabs(res.getOrd(QCD2, QED1)(a, b)) < EPS) res.setMatrixElement(QCD2, QED1, a, b, 0.);
            }
    } else
        for (uint a = 0; a < nops; a++)
            for (uint b = 0; b < nops; b++)
            {
                if (fabs(res.getOrd(QCD0, QED0)(a, b)) < EPS) res.setMatrixElement(QCD0, QED0, a, b, 0.);
                if (fabs(res.getOrd(QCD1, QED0)(a, b)) < EPS) res.setMatrixElement(QCD1, QED0, a, b, 0.);
                if (fabs(res.getOrd(QCD2, QED0)(a, b)) < EPS) res.setMatrixElement(QCD2, QED0, a, b, 0.);
            }

    wilson = res * wilson;
}

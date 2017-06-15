/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <boost/assign.hpp>
#include <initializer_list>
#include "QCD.h"
#include "EvolDF1.h"

double EvolDF1::Beta_s(int i, double nf)
{
    switch(i)
    {
        case 00:
            return(QCD::Beta0(nf));
        case 10:
            return(QCD::Beta1(nf));
        case 20:
            return(QCD::Beta2(nf));
        case 30:
            return(QCD::Beta3(nf));
        case 01:
            if (nf == 5) return(-22./9.);
            else throw std::runtime_error("EvolDF1::Beta_s(01) only known for nf=5.");
        case 11:
            if (nf == 5) return(-308./27.);
            else throw std::runtime_error("EvolDF1::Beta_s(11) only known for nf=5.");
        case 02:
            if (nf == 5) return(4946./243.);
            else throw std::runtime_error("EvolDF1::Beta_s(02) only known for nf=5.");
        default:
            throw std::runtime_error("EvolDF1::Beta_s: case not implemented");
    }
}
double EvolDF1::Beta_e(int i, double nf)
{
    if (nf != 5) throw std::runtime_error("EvolDF1::Beta_e only known for nf=5.");

    switch(i)
    {
        case 00:
            return(80./9.);
        case 10:
            return(464./27.);
        case 01:
            return(176./9.);
        default:
            throw std::runtime_error("EvolDF1::Beta_e: case not implemented");
    }
}

EvolDF1::EvolDF1(unsigned int nops, std::string reqblocks, schemes scheme, orders order, const StandardModel& model) 
:           RGEvolutor(nops, scheme, order), 
            a(boost::extents[4][nops]), b(boost::extents[4][nops][nops][nops]), c(boost::extents[4][nops][nops][nops]), d(boost::extents[4][nops][nops][nops]),
            model(model), blocks(reqblocks),
            v(nops,0.), vi(nops,0.), js(nops,0.), h(nops,0.), gg(nops,0.), s_s(nops,0.),
            jssv(nops,0.), jss(nops,0.), jv(nops,0.), vij(nops,0.), e(nops,0.)
{
    blocks_nops = boost::assign::map_list_of("C",2) ("CP",6) ("CPM",8) ("L",2) ("CPML",10) ("CPQB",11) ("CPMQB",13) ("CPMLQB",15);
//    blocks_nops = {{"C",2},{"CP",6},{"CPM",8},{"L",2},{"CPML",10},{"CPQB",11},{"CPMQB",13},{"CPMLQB",15}};
    blocks_ord = boost::assign::map_list_of("C",NNLO) ("CP",NNLO) ("CPM",NNLO) ("L",NNLO) ("CPML",NNLO) ("CPQB",NLO) ("CPMQB",NLO) ("CPMLQB",NLO);
//    blocks_ord = {{"C",NNLO},{"CP",NNLO},{"CPM",NNLO},{"L",NNLO},{"CPML",NNLO},{"CPQB",NLO},{"CPMQB",NLO},{"CPMLQB",NLO}};
    
    if (blocks_nops[blocks] != nops)
        throw std::runtime_error("EvolDF1(): number of operators does not match block specification");

    unsigned int i, j, k;
    this->nops = nops;
    
    /* magic numbers a & b */ 
    
    for(int L=3; L>-1; L--){
        
    /* L=3 --> u,d,s (nf=3) L=2 --> u,d,s,c (nf=4)  L=1 --> u,d,s,c,b (nf=5) L=0 --> u,d,s,c,b,t (nf=6) */        
    nu = L;  nd = L;
    if(L == 3){nd = 2; nu = 1;} 
    if(L == 1){nd = 3; nu = 2;} 
    if(L == 0){nd = 3; nu = 3;}

    // LO evolutor of the effective Wilson coefficients in the Chetyrkin, Misiak and Munz basis
    
    AnomalousDimension_s(LO,nu,nd).transpose().eigensystem(v,e);
    vi = v.inverse();
    for(i = 0; i < nops; i++){
       a[L][i] = e(i).real();
       for (j = 0; j < nops; j++) {
           for (k = 0; k < nops; k++)  {
                b[L][i][j][k] = v(i, k).real() * vi(k, j).real();
               }
           }
       }
    
    // NLO evolutor of the effective Wilson coefficients in the Chetyrkin, Misiak and Munz basis
    
    gg = vi * AnomalousDimension_s(NLO,nu,nd).transpose() * v;
    double b0 = model.Beta0(nu+nd);
    double b1 = model.Beta1(nu+nd);
    for (i = 0; i < nops; i++){
        for (j = 0; j < nops; j++){
            s_s.assign( i, j, (b1 / b0) * (i==j) * e(i).real() - gg(i,j));    
            if(fabs(e(i).real() - e(j).real() + 2. * b0)>0.00000000001){
                h.assign( i, j, s_s(i,j) / (2. * b0 + e(i) - e(j)));
                }
            else
                throw std::runtime_error("HeffDF1:: to be fixed");
            }
        }
    js = v * h * vi;
    jv = js * v;
    vij = vi * js;
    jss = v * s_s * vi;
    jssv = jss * v;        
    for (i = 0; i < nops; i++){
        for (j = 0; j < nops; j++){
            if(fabs(e(i).real() - e(j).real() + 2. * b0) > 0.00000000001){
                for(k = 0; k < nops; k++){
                        c[L][i][j][k] = jv(i, k).real() * vi(k, j).real();
                        d[L][i][j][k] = -v(i, k).real() * vij(k, j).real();
                        }
                    }
            else{    
                for(k = 0; k < nops; k++){
                   c[L][i][j][k] = (1./(2. * b0)) * jssv(i, k).real() * vi(k, j).real();
                   d[L][i][j][k] = 0.;
                   }   
                }
            }
        }
    }
}
    
EvolDF1::~EvolDF1() 
{}

/* Delta F = 1 anomalous dimension in Misiak basis, in the effective basis (C7eff, C8eff)
       ref. for CC,CP,PP QED NLO, QP,QQ,BP,BB QCD NLO, CL,PL NNLO, PQ,QL,LL,BL NLO: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
       ref. for CC,CP,PP QCD NNLO with nf: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
       ref. for CM,PM QCD NNLO with nf: Czakon, Haisch, Misiak, JHEP 0703, 008, hep-ph/0612329
       ref. for MM QCD NNLO with nf: Gorbahn, Haisch, Misiak, Phys. Rev. Lett. 95, 102004, hep-ph/0504194
       ref. for CM,PM QED LO, QM QCD LO: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
       ref. for MM,QP,QQ,BP,BB QED NLO, BQ NLO: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
                       QM QED,BM??? (in 031209 C7,8 are normalised with alpha_s and not in the effective basis)
*/

gslpp::matrix<double> EvolDF1::GammaCC_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // ref.: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){
            throw std::runtime_error("EvolDF1::GammaCC_s(): wrong number of flavours");
    }
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -4.;
            gammaDF1(0,1) = 8./3.;
            gammaDF1(1,0) = 12.;
            break;
        case NLO:
            gammaDF1(0,0) = -145./3. + nf*16./9.;// -355./9.
            gammaDF1(0,1) = -26. + nf*40./27.;   // -502./27.
            gammaDF1(1,0) = -45. + nf*20./3.;    // -35./3.
            gammaDF1(1,1) = -28./3.;
            break;
        case NNLO:
            gammaDF1(0,0) = -1927./2. + nf*257./9. + nf*nf*40./9. + z3*(224. + nf*160./3.);   // -12773./18. + z3*1472./3.
            gammaDF1(0,1) = 475./9. + nf*362./27. - nf*nf*40./27. - z3*(896./3. + nf*320./9.);// 745./9. - z3*4288./9.
            gammaDF1(1,0) = 307./2. + nf*361./3. - nf*nf*20./3. - z3*(1344. + nf*160.);       // 1177./2. - z3*2144.
            gammaDF1(1,1) = 1298./3. - nf* 76./3. - z3*224.;                                  // 306. + z3*224.
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCC_s(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCP_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // ref.: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
    gslpp::matrix<double> gammaDF1(2, 4, 0.);
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){
            throw std::runtime_error("EvolDF1::GammaCP_s(): wrong number of flavours");
    }
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,1) = -2./9.;
            gammaDF1(1,1) = 4./3.;
            break;
        case NLO:
            gammaDF1(0,0) = -1412./243.;
            gammaDF1(0,1) = -1369./243.;
            gammaDF1(0,2) = 134./243.;
            gammaDF1(0,3) = -35./162.;
            gammaDF1(1,0) = -416./81.;
            gammaDF1(1,1) = 1280./81.;
            gammaDF1(1,2) = 56./81.;
            gammaDF1(1,3) = 35./27.;
            break;
        case NNLO:
            gammaDF1(0,0) = 269107./13122. - nf*2288./729. - z3*1360./81.;   // 63187./13122. - z3*1360./81.
            gammaDF1(0,1) = -2425817./13122. + nf*30815./4374. - z3*776./81.;// -981796./6561. - z3*776./81.
            gammaDF1(0,2) = -343783./52488. + nf*392./729. + z3*124./81.;    // -202663./52488. + z3*124./81.
            gammaDF1(0,3) = -37573./69984. + nf*35./972. + z3*100./27.;      // -24973./69984. + z3*100./27.
            gammaDF1(1,0) = 69797./2187. + nf*904./243. + z3*2720./27.;      // 110477./2187. + z3*2720./.27.
            gammaDF1(1,1) = 1457549./8748. - nf*22067./729. - z3*2768./27.;  // 133529./8748. - z3*2768./27.
            gammaDF1(1,2) = -37889./8748. - nf*28./243. - z3*248./27.;       // -42929./8748. - z3*248./27.
            gammaDF1(1,3) = 366919./11664. - nf*35./162. - z3*110./9.;       // 354319./11664. - z3*110./9.
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCP_s(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCM_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // ref.: Czakon, Haisch, Misiak, JHEP 0703, 008, hep-ph/0612329
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){
            throw std::runtime_error("EvolDF1::GammaCM_s(): wrong number of flavours");
    }
    double Qu = 2./3.;
    double Qd = -1./3.;
    double Qbar = nu*Qu + nd*Qd;
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 8./243. - Qu*4./3.;
            gammaDF1(0,1) = 173./162.;
            gammaDF1(1,0) = -16./81. + Qu*8.;
            gammaDF1(1,1) = 70./27.;
            break;
        case NLO:
            gammaDF1(0,0) = 12614./2187. - nf*64./2187. - Qu*374./27. + nf*Qu*2./27.;
            gammaDF1(0,1) = 65867./5832. + nf*431./5832.;
            gammaDF1(1,0) = -2332./729. + nf*128./729. + Qu*136./9. - nf*Qu*4./9.;
            gammaDF1(1,1) = 10577./486. - nf*917./972.;
            break;
        case NNLO:
            gammaDF1(0,0) = 77506102./531441. - nf*875374./177147. + nf*nf*560./19683. - Qu*9731./162. +
                    nf*Qu*11045./729. + nf*nf*Qu*316./729. + Qbar*3695./486. + z3*(-112216./6561. + nf*728./729. +
                    Qu*25508./81. - nf*Qu*64./81. - Qbar*100./27.);
            gammaDF1(0,1) = -421272953./1417176. - nf*8210077./472392. - nf*nf*1955./6561. + z3*(-953042./2187. -
                    nf*10381./486.);
            gammaDF1(1,0) = -15463055./177147. + nf*242204./59049. - nf*nf*1120./6561. + Qu*55748./27. -
                    nf*Qu*33970./243. - nf*nf*Qu*632./243. - Qbar*3695./81. + z3*(365696./2187. - nf*1168./243. -
                    Qu*51232./27. - nf*Qu*1024./27. + Qbar*200./9.);
            gammaDF1(1,1) = 98548513./472392. - nf*5615165./78732. - nf*nf*2489./2187. + z3*(-607103./729. -
                    nf*1679./81.);
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCM_s(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPP_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // ref.: Gorbahn, Haisch, Nucl. Phys. B 713, 291, hep-ph/0411071
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){
            throw std::runtime_error("EvolDF1::GammaPP_s(): wrong number of flavours");
    }
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,1) = -52./3.;
            gammaDF1(0,3) = 2.;
            gammaDF1(1,0) = -40./9.;
            gammaDF1(1,1) = -160./9. + nf*4./3.;// -100./9.
            gammaDF1(1,2) = 4./9.;
            gammaDF1(1,3) = 5./6.;
            gammaDF1(2,1) = -256./3.;
            gammaDF1(2,3) = 20.;
            gammaDF1(3,0) = -256./9.;
            gammaDF1(3,1) = -544./9. + nf*40./3.;// 56./9.
            gammaDF1(3,2) = 40./9.;
            gammaDF1(3,3) = -2./3.;            
            break;
        case NLO:
            gammaDF1(0,0) = -4468./81.;
            gammaDF1(0,1) = -29129./81. - nf*52./9.;    // -31469./81.
            gammaDF1(0,2) = 400./81.;
            gammaDF1(0,3) = 3493./108. - nf*2./9.;      // 3373./108.
            gammaDF1(1,0) = -13678./243. + nf*368./81;  // -8158./243. 
            gammaDF1(1,1) = -79409./243. + nf*1334./81.;// -59399./243.
            gammaDF1(1,2) = 509./486. - nf*8./81.;      // 269./486.
            gammaDF1(1,3) = 13499./648. - nf*5./27.;    // 12899./648.
            gammaDF1(2,0) = -244480./81. - nf*160./9.;  // -251680./81.
            gammaDF1(2,1) = -29648./81. - nf*2200./9.;  // -128648./81.
            gammaDF1(2,2) = 23116./81. + nf*16./9.;     // 23836./81.
            gammaDF1(2,3) = 3886./27. + nf*148./9.;     // 6106./27.
            gammaDF1(3,0) = 77600./243. - nf*1264./81.; // 58640./243.
            gammaDF1(3,1) = -28808./243. - nf*164./81.; // -26348./243.
            gammaDF1(3,2) = -20324./243. + nf*400./81.; // -14324./243.
            gammaDF1(3,3) = -21211./162. + nf*622./27.; // -2551./162.
            break;
        case NNLO:
            gammaDF1(0,0) = -4203068./2187. + nf*14012./243. - z3*608./27.;                                            // -3572528./2187. - z3*608./27.
            gammaDF1(0,1) = -18422762./2187. + nf*888605./2916. + nf*nf*272./27. + z3*(39824./27. + nf*160.);          // -58158773./8748. + z3*61424./27.
            gammaDF1(0,2) = 674281./4374. - nf*1352./243. - z3*496./27.;                                               // 552601./4374. - z3*496./27.
            gammaDF1(0,3) = 9284531./11664. - nf*2798./81. - nf*nf*26./27. - z3*(1921./9. + nf*20.);                   // 6989171./11664. - z3*2821./9.
            gammaDF1(1,0) = -5875184./6561. + nf*217892./2187. + nf*nf*472./81. + z3*(27520./81. + nf*1360./9.);       // -1651004./6561. + z3*88720./81.
            gammaDF1(1,1) = -70274587./13122. + nf*8860733./17496. - nf*nf*4010./729. + z3*(16592./81. + nf*2512./27.);// -155405353./52488. + z3*54272./81.
            gammaDF1(1,2) = 2951809./52488. - nf*31175./8748. - nf*nf*52./81. - z3*(3154./81. + nf*136./9.);           // 1174159./52488. - z3*9274./81.
            gammaDF1(1,3) = 3227801./8748. - nf*105293./11664. - nf*nf*65./54. + z3*(200./27. - nf*220./9.);           // 10278809./34992. - z3*3100./27.
            gammaDF1(2,0) = -194951552./2187. + nf*358672./81. - nf*nf*2144./81. + z3*87040./27.;                      // -147978032./2187. + z3*87040./27.
            gammaDF1(2,1) = -130500332./2187. - nf*2949616./729. + nf*nf*3088./27. + z3*(238016./27. + nf*640.);       // -168491372./2187. + z3*324416./27.
            gammaDF1(2,2) = 14732222./2187. - nf*27428./81. + nf*nf*272./81. - z3*13984./27.;                          // 11213042./2187. - z3*13984./27.
            gammaDF1(2,3) = 16521659./2916. + nf*8081./54. - nf*nf*316./27. - z3*(22420./9. + nf*200.);                // 17850329./2916. - z3*31420./9.
            gammaDF1(3,0) = 162733912./6561. - nf*2535466./2187. + nf*nf*17920./243. + z3*(174208./81. + nf*12160./9.);// 136797922./6561. + z3*721408./81.
            gammaDF1(3,1) = 13286236./6561. - nf*1826023./4374. - nf*nf*159548./729. - z3*(24832./81. + nf*9440./27.); // -72614473./13122. - z3*166432./81.
            gammaDF1(3,2) = -22191107./13122. + nf*395783./4374. - nf*nf*1720./243. - z3*(33832./81. + nf*1360./9.);   // -9288181./6561. - z3*95032./81.
            gammaDF1(3,3) = -32043361./8748. + nf*3353393./5832. - nf*nf*533./81. + z3*(9248./27. - nf*1120./9.);      // -16664027./17496. - z3*7552./27.
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPP_s(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPM_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // ref.: Czakon, Haisch, Misiak, JHEP 0703, 008, hep-ph/0612329
    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){
            throw std::runtime_error("EvolDF1::GammaPM_s(): wrong number of flavours");
    }
    double Qu = 2./3.;
    double Qd = -1./3.;
    double Qbar = nu*Qu + nd*Qd;
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -176./81.;
            gammaDF1(0,1) = 14./27.;
            gammaDF1(1,0) = 88./243. - nf*16./81.;
            gammaDF1(1,1) = 74./81. - nf*49./54.;
            gammaDF1(2,0) = -6272./81.;
            gammaDF1(2,1) = 1736./27. + nf*36.;
            gammaDF1(3,0) = 3136./243. - nf*160./81. + Qbar*48.;
            gammaDF1(3,1) = 2372./81. + nf*160./27.;
            break;
        case NLO:
            gammaDF1(0,0) = 97876./729. - nf*4352./729. - Qbar*112./3.;
            gammaDF1(0,1) = 42524./243. - nf*2398./243.;
            gammaDF1(1,0) = -70376./2187. - nf*15788./2187. + nf*nf*32./729. - Qbar*140./9.;
            gammaDF1(1,1) = -159718./729. - nf*39719./5832. - nf*nf*253./486.;
            gammaDF1(2,0) = 1764752./729. - nf*65408./729. - Qbar*3136./3.;
            gammaDF1(2,1) = 2281576./243. + nf*140954./243. - nf*nf*14.;
            gammaDF1(3,0) = 4193840./2187. - nf*324128./2187. + nf*nf*896./729. - Qbar*1136./9. - nf*Qbar*56./3.;
            gammaDF1(3,1) = -3031517./729. - nf*15431./1458. - nf*nf*6031./486.;
            break;
        case NNLO:
            gammaDF1(0,0) = 102439553./177147. - nf*12273398/59049. + nf*nf*5824./6561. + Qbar*26639./81. - nf*Qbar*8./27. +
                    z3*(3508864./2187. - nf*1904./243. - Qbar*1984./9. - nf*Qbar*64./9.);
            gammaDF1(0,1) = 3205172129./472392. - nf*108963529./314928. + nf*nf*58903./4374. + z3*(-1597588./729. +
                    nf*13028./81. - nf*nf*20./9.);
            gammaDF1(1,0) = -2493414077./1062882. - nf*9901031./354294. + nf*nf*243872./59049. - nf*nf*nf*1184./6561. -
                    Qbar*49993./972. + nf*Qbar*305./27. + z3*(-1922264./6561. + nf*308648./2187. - nf*nf*1280./243. +
                    Qbar*1010./9. - nf*Qbar*200./27.);
            gammaDF1(1,1) = -6678822461./2834352. + nf*127999025./1889568. + nf*nf*1699073./157464. + nf*nf*nf*505./4374. +
                    z3*(2312684./2187. + nf*128347./729. + nf*nf*920./81.);
            gammaDF1(2,0) = 8808397748./177147. - nf*174839456./59049. + nf*nf*1600./729. - Qbar*669694./81. + nf*Qbar*10672./27. +
                    z3*(123543040./2187. - nf*207712./243. + nf*nf*128./27. - Qbar*24880./9. - nf*Qbar*640./9.);
            gammaDF1(2,1) = 29013624461./118098. - nf*64260772./19683. - nf*nf*230962./243. - nf*nf*nf*148./27. +
                    z3*(-69359224./729. - nf*885356./81. - nf*nf*5080./9.);
            gammaDF1(3,0) = 7684242746./531441. - nf*351775414./177147. - nf*nf*479776./59049. - nf*nf*nf*11456./6561. +
                    Qbar*3950201./243. - nf*Qbar*130538./81. - nf*nf*Qbar*592./81. + z3*(7699264./6561. + nf*2854976./2187. -
                    nf*nf*12320./243. - Qbar*108584./9. - nf*Qbar*1136./27.);
            gammaDF1(3,1) = -72810260309./708588. + nf*2545824851./472392. - nf*nf*33778271./78732. - nf*nf*nf*3988./2187. +
                    z3*(-61384768./2187. - nf*685472./729. + nf*nf*350./81.);
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPM_s(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaMM_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    //ref.: Gorbahn, Haisch, Misiak, Phys. Rev. Lett. 95, 102004, hep-ph/0504194
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){
            throw std::runtime_error("EvolDF1::GammaPM_s(): wrong number of flavours");
    }
    double Qu = 2./3.;
    double Qd = -1./3.;
    double Qbar = nu*Qu + nd*Qd;
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 32./3.;
            gammaDF1(1,0) = Qd*32./3.;
            gammaDF1(1,1) = 28./3.;
            break;
        case NLO:
            gammaDF1(0,0) = 1936./9. - nf*224./27.;
            gammaDF1(1,0) = Qd*368./3. - nf*Qd*224./27.;
            gammaDF1(1,1) = 1456./9. - nf*61./27.;
            break;
        case NNLO:
            gammaDF1(0,0) = 307448./81. - nf*23776./81. - nf*nf*352./81. + z3*(1856./27. - nf*1280./9.);
            gammaDF1(1,0) = -Qbar*1600./27. + Qd*159872./81. - nf*Qd*17108./81. - nf*nf*Qd*352./81. + z3*(Qbar*640./9. -
                    Qd*1856./27. + nf*Qd*1280./9.);
            gammaDF1(1,1) = 268807./81. - nf*4343./27. - nf*nf*461./81. + z3*(-28624./27. - nf*1312./9.);
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaMM_s(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQP_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    //ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,1) = -8./9.;
            gammaDF1(1,1) = 16./27.;
            gammaDF1(2,1) = -128./9.;
            gammaDF1(3,1) = 184./27.;          
            break;
        case NLO:
            gammaDF1(0,0) = 832./243.;
            gammaDF1(0,1) = -4000./243.;
            gammaDF1(0,2) = -112./243.;
            gammaDF1(0,3) = -70./81.;
            gammaDF1(1,0) = 3376./729.;
            gammaDF1(1,1) = 6344./729.;
            gammaDF1(1,2) = -280./729.;
            gammaDF1(1,3) = 55./486.;
            gammaDF1(2,0) = 2272./243.;
            gammaDF1(2,1) = -72088./243.;
            gammaDF1(2,2) = -688./243.;
            gammaDF1(2,3) = -1240./81.;
            gammaDF1(3,0) = 45424./729.;
            gammaDF1(3,1) = 84236./729.;
            gammaDF1(3,2) = -3880./729.;
            gammaDF1(3,3) = 1220./243.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQP_s(): order not implemented");       
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQM_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 176./243.;
            gammaDF1(0,1) = -14./81.;
            gammaDF1(1,0) = -136./729.;
            gammaDF1(1,1) = -295./486.;
            gammaDF1(2,0) = 6272./243.;
            gammaDF1(2,1) = -764./81.;
            gammaDF1(3,0) = 39152./729.;
            gammaDF1(3,1) = -1892./243.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQM_s(): order not implemented");        
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQQ_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,1) = -20.;
            gammaDF1(0,3) = 2.;
            gammaDF1(1,0) = -40./9.;
            gammaDF1(1,1) = -52./3.;
            gammaDF1(1,2) = 4./9.;
            gammaDF1(1,3) = 5./6.;
            gammaDF1(2,1) = -128.;
            gammaDF1(2,3) = 20.;
            gammaDF1(3,0) = -256./9.;
            gammaDF1(3,1) = 160./3.;
            gammaDF1(3,2) = 40./9.;
            gammaDF1(3,3) = -2./3.;        
            break;
        case NLO:
            gammaDF1(0,0) = -404./9.;
            gammaDF1(0,1) = -3077./9.;
            gammaDF1(0,2) = 32./9.;
            gammaDF1(0,3) = 1031./36.;
            gammaDF1(1,0) = -2698./81.;
            gammaDF1(1,1) = -8035./27.;
            gammaDF1(1,2) = -49./162.;
            gammaDF1(1,3) = 4493./216.;
            gammaDF1(2,0) = -19072./9.;
            gammaDF1(2,1) = -14096./9.;
            gammaDF1(2,2) = 1708./9.;
            gammaDF1(2,3) = 1622./9.;
            gammaDF1(3,0) = 32288./81.;
            gammaDF1(3,1) = -15976./27.;
            gammaDF1(3,2) = -6692./81.;
            gammaDF1(3,3) = -2437./54.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQQ_s(): order not implemented");         
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBP_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(1, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,1) = -4./3.;       
            break;
        case NLO:
            gammaDF1(0,0) = -1576./81.;
            gammaDF1(0,1) = 446./27.;
            gammaDF1(0,2) = 172./81.;
            gammaDF1(0,3) = 40./27.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBP_s(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBB_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(1, 1, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 4.;       
            break;
        case NLO:
            gammaDF1(0,0) = 325./9.;
        default:
            throw std::runtime_error("EvolDF1::GammaBB_s(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::AnomalousDimension_s(orders order, unsigned int n_u, unsigned int n_d) const
{
    gslpp::matrix<double> gammaDF1(nops, nops, 0.);

// assign blocks according to user request: "C", "CP", "CPM", "L", "CPML", "CPQB", "CPMQB", "CPMLQB"
    if(blocks[0] == 'C') {
        gammaDF1.assign(0,0,GammaCC_s(order, n_u, n_d));
        if(blocks.size() > 1) {
            if (blocks[1] == 'P') {
                int m = blocks_nops.find("C")->second;
                gammaDF1.assign(0,m,GammaCP_s(order, n_u, n_d));
                gammaDF1.assign(blocks_nops.find("C")->second,blocks_nops.find("C")->second,GammaPP_s(order, n_u, n_d));
                if(blocks.size() > 2) {
                    if (blocks[2] == 'M') {
                    gammaDF1.assign(0,blocks_nops.find("CP")->second,GammaCM_s(order, n_u, n_d));
                    gammaDF1.assign(blocks_nops.find("C")->second,blocks_nops.find("CP")->second,GammaPM_s(order, n_u, n_d));
                    gammaDF1.assign(blocks_nops.find("CP")->second,blocks_nops.find("CP")->second,GammaMM_s(order, n_u, n_d));
                    } else if (blocks[2] == 'Q') {
                        gammaDF1.assign(blocks_nops.find("CP")->second,blocks_nops.find("C")->second,GammaQP_s(order, n_u, n_d));
                        gammaDF1.assign(blocks_nops.find("CP")->second,blocks_nops.find("CP")->second,GammaQQ_s(order, n_u, n_d));
                        if(blocks[3] == 'B') {
                            gammaDF1.assign(blocks_nops.find("CPQ")->second,blocks_nops.find("C")->second,GammaBP_s(order, n_u, n_d));
                            gammaDF1.assign(blocks_nops.find("CPQ")->second,blocks_nops.find("CPQ")->second,GammaBB_s(order, n_u, n_d));
                        }
                    }
                    if(blocks.size() > 3) {
                        if (blocks[3] == 'L') {
                            if(blocks.size() > 4) {
                                if (blocks[4] == 'Q') {
                                    gammaDF1.assign(blocks_nops.find("CPML")->second,blocks_nops.find("C")->second,GammaQP_s(order, n_u, n_d));
                                    gammaDF1.assign(blocks_nops.find("CPML")->second,blocks_nops.find("CPML")->second,GammaQQ_s(order, n_u, n_d));
                                    if(blocks[5] == 'B') {
                                        gammaDF1.assign(blocks_nops.find("CPMLQ")->second,blocks_nops.find("C")->second,GammaBP_s(order, n_u, n_d));
                                        gammaDF1.assign(blocks_nops.find("CPMLQ")->second,blocks_nops.find("CPMLQ")->second,GammaBB_s(order, n_u, n_d));
                                    }
                                }
                            }
                        } else if (blocks[3] == 'Q') {
                            gammaDF1.assign(blocks_nops.find("CPM")->second,blocks_nops.find("C")->second,GammaQP_s(order, n_u, n_d));
                            gammaDF1.assign(blocks_nops.find("CPM")->second,blocks_nops.find("CPM")->second,GammaQQ_s(order, n_u, n_d));
                            if(blocks[4] == 'B') {
                                gammaDF1.assign(blocks_nops.find("CPMQ")->second,blocks_nops.find("C")->second,GammaBP_s(order, n_u, n_d));
                                gammaDF1.assign(blocks_nops.find("CPMQ")->second,blocks_nops.find("CPMQ")->second,GammaBB_s(order, n_u, n_d));
                            }                     
                        }                            
                    }
                }
            }
        }
    }
                    
    return (gammaDF1);
}


gslpp::matrix<double> EvolDF1::GammaCC_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -8./3.;
            gammaDF1(1,1) = -8./3.;
            break;
        case NLO:
            gammaDF1(0,0) = 169./9.;
            gammaDF1(0,1) = 100./27.;
            gammaDF1(1,0) = 50./3.;
            gammaDF1(1,1) = -8./3.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCC_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCP_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(2, 4, 0.);
    
    switch(order)
    {  
        case LO:
            break;
        case NLO:
            gammaDF1(0,1) = 254./729.;
            gammaDF1(1,1) = 1076./243.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCP_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCM_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -832./729.;
            gammaDF1(0,1) = 22./243.;
            gammaDF1(1,0) = -208./243.;
            gammaDF1(1,1) = -116./81.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCM_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCL_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -32./27.;
            gammaDF1(1,0) = -8./9.;
            break;
        case NLO:
            gammaDF1(0,0) = -2272./729.;
            gammaDF1(1,0) = 1952./243.;
            break;
        case NNLO:
            gammaDF1(0,0) = -1359190./19683. + z3*6976./243.;
            gammaDF1(1,0) = -229696./6561. - z3*3584./81.;
        default:
            throw std::runtime_error("EvolDF1::GammaCL_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaCQ_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(2, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 32./27.;
            gammaDF1(1,0) = 8./9.;
            break;
        case NLO:
            gammaDF1(0,0) = 2272./729.;
            gammaDF1(0,1) = 122./81.;
            gammaDF1(0,3) = 49./81.;
            gammaDF1(1,0) = -1952./243.;
            gammaDF1(1,1) = -748./27.;
            gammaDF1(1,3) = 82./27.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaCQ_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPP_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    
    switch(order)
    {  
        case LO:
            break;
        case NLO:
            gammaDF1(0,1) = 11116./243.;
            gammaDF1(0,3) = -14./3.;
            gammaDF1(1,0) = 280./27.;
            gammaDF1(1,1) = 18763./729.;
            gammaDF1(1,2) = -28./27.;
            gammaDF1(1,3) = -35./18.;
            gammaDF1(2,1) = 111136./243.;
            gammaDF1(2,3) = -140./3.;
            gammaDF1(3,0) = 2944./27.;
            gammaDF1(3,1) = 193312./729.;
            gammaDF1(3,2) = -280./27.;
            gammaDF1(3,3) = -175./9.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPP_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPM_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Baranowski, Misiak, Phys. Lett. B 483, 410, hep-ph/9907427
    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -20./243.;
            gammaDF1(0,1) = 20./81.;
            gammaDF1(1,0) = -176./729.;
            gammaDF1(1,1) = 14./243.;
            gammaDF1(2,0) = -22712./243.;
            gammaDF1(2,1) = 1328./81.;
            gammaDF1(3,0) = -6272./729.;
            gammaDF1(3,1) = -1180./243.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPM_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPL_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    double z3 = gslpp_special_functions::zeta(3);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -16./9.;
            gammaDF1(1,0) = 32./27.;
            gammaDF1(2,0) = -112./9.;
            gammaDF1(3,0) = 512./27.;
            break;
        case NLO:
            gammaDF1(0,0) = -6752./243.;
            gammaDF1(1,0) = -2192./729.;
            gammaDF1(2,0) = -84032./243.;
            gammaDF1(3,0) = -37856./729.;
            break;
        case NNLO:
            gammaDF1(0,0) = -1290092./6561. + z3*3200./81.;
            gammaDF1(1,0) = -819971./19673. - z3*19936./243.;
            gammaDF1(2,0) = -16821944./6561. + z3*30464./81.;
            gammaDF1(3,0) = -17787268./19683. - z3*286720./243.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPL_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaPQ_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 76./9.;
            gammaDF1(0,2) = -2./3.;
            gammaDF1(1,0) = -32./27.;
            gammaDF1(1,1) = 20./3.;
            gammaDF1(1,3) = -2./3.;
            gammaDF1(2,0) = 496./9.;
            gammaDF1(2,2) = -20./3.;
            gammaDF1(3,0) = -512./27.;
            gammaDF1(3,1) = 128./3.;
            gammaDF1(3,3) = -20./3.;
            break;
        case NLO:
            gammaDF1(0,0) = -23488./243.;
            gammaDF1(0,1) = 6280./27.;
            gammaDF1(0,2) = 112./9.;
            gammaDF1(0,3) = -538./27.;
            gammaDF1(1,0) = 31568./729.;
            gammaDF1(1,1) = 9481./81.;
            gammaDF1(1,2) = -92./27.;
            gammaDF1(1,3) = -1012./81.;
            gammaDF1(2,0) = -233920./243.;
            gammaDF1(2,1) = 68848./27.;
            gammaDF1(2,2) = 1120./9.;
            gammaDF1(2,3) = -5056./27.;
            gammaDF1(3,0) = 352352./729.;
            gammaDF1(3,1) = 116680./81.;
            gammaDF1(3,2) = -752./27.;
            gammaDF1(3,3) = -10147./81.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaPQ_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaMM_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 16./9.;
            gammaDF1(0,1) = -8./3.;
            gammaDF1(1,1) = 8./9.;
            break;
        case NLO:
            gammaDF1(0,0) = -256./27.;
            gammaDF1(0,1) = -52./9.;
            gammaDF1(1,0) = 128./81.;
            gammaDF1(1,1) = -40./27.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaMM_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaLL_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(2, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 8.;
            gammaDF1(0,1) = -4.;
            gammaDF1(1,0) = -4.;
            break;
        case NLO:
            gammaDF1(0,1) = 16.;
            gammaDF1(1,0) = 16.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaLL_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQP_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 40./27.;
            gammaDF1(0,2) = -4./27.;
            gammaDF1(1,1) = 40./27.;
            gammaDF1(1,3) = -4./27.;
            gammaDF1(2,0) = 256./27.;
            gammaDF1(2,2) = -40./27.;
            gammaDF1(3,1) = 256./27.;
            gammaDF1(3,3) = -40./27.;
            break;
        case NLO:
            gammaDF1(0,0) = -2240./81.;
            gammaDF1(0,1) = 39392./729.;
            gammaDF1(0,2) = 224./81.;
            gammaDF1(0,3) = -92./27.;
            gammaDF1(1,0) = 2176./243.;
            gammaDF1(1,1) = 84890./2187.;
            gammaDF1(1,2) = -184./243.;
            gammaDF1(1,3) = -224./81.;
            gammaDF1(2,0) = -23552./81.;
            gammaDF1(2,1) = 399776./729.;
            gammaDF1(2,2) = 2240./81.;
            gammaDF1(2,3) = -752./27.;
            gammaDF1(3,0) = 23296./243.;
            gammaDF1(3,1) = 933776./2187.;
            gammaDF1(3,2) = -1504./243.;
            gammaDF1(3,3) = -2030./81.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQP_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQL_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(4, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -272./27.;
            gammaDF1(1,0) = -32./81.;
            gammaDF1(2,0) = -2768./27.;
            gammaDF1(3,0) = -512./81.;
            break;
        case NLO:
            gammaDF1(0,0) = -24352./729.;
            gammaDF1(1,0) = 54608./2187.;
            gammaDF1(2,0) = -227008./729.;
            gammaDF1(3,0) = 551648./2187.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQL_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaQQ_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
    gslpp::matrix<double> gammaDF1(4, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 332./27.;
            gammaDF1(0,2) = -2./9.;
            gammaDF1(1,0) = 32./81.;
            gammaDF1(1,1) = 20./9.;
            gammaDF1(1,3) = -2./9.;
            gammaDF1(2,0) = 3152./27.;
            gammaDF1(2,2) = -20./9.;
            gammaDF1(3,0) = 512./81.;
            gammaDF1(3,1) = 128./9.;
            gammaDF1(3,3) = -20./9.;
            break;
        case NLO:
            gammaDF1(0,0) = -5888./729.;
            gammaDF1(0,1) = 13916./81.;
            gammaDF1(0,2) = 112./27.;
            gammaDF1(0,3) = -812./81.;
            gammaDF1(1,0) = -2552./2187.;
            gammaDF1(1,1) = 15638./243.;
            gammaDF1(1,2) = -176./81.;
            gammaDF1(1,3) = -2881./486.;
            gammaDF1(2,0) = -90944./729.;
            gammaDF1(2,1) = 90128./81.;
            gammaDF1(2,2) = 1120./27.;
            gammaDF1(2,3) = -1748./81.;
            gammaDF1(3,0) = 1312./2187.;
            gammaDF1(3,1) = 102488./243.;
            gammaDF1(3,2) = -1592./81.;
            gammaDF1(3,3) = -6008./243.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaQQ_e(): order not implemented");
    }
     return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBP_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
    gslpp::matrix<double> gammaDF1(1, 4, 0.);
    
    switch(order)
    {  
        case LO:     
            break;
        case NLO:
            gammaDF1(0,1) = -232./81.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBP_e(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBL_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Huber, Lunghi, Misiak, Wyler, Nucl. Phys. B 740, 105, hep-ph/0512066
    gslpp::matrix<double> gammaDF1(1, 2, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 16./9.;
            break;
        case NLO:
            gammaDF1(0,0) = -8./9.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBL_e(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBQ_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
    gslpp::matrix<double> gammaDF1(1, 4, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = -16./9.;
            break;
        case NLO:
            gammaDF1(0,1) = 580./27.;
            gammaDF1(0,3) = -94./27.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBQ_e(): order not implemented");
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDF1::GammaBB_e(orders order, unsigned int n_u, unsigned int n_d) const
{
    // only available for nf = 5
    // ref.: Bobeth, Gambino, Gorbahn, Haisch, JHEP 0404, 071, hep-ph/0312090
    gslpp::matrix<double> gammaDF1(1, 1, 0.);
    
    switch(order)
    {  
        case LO:
            gammaDF1(0,0) = 4./3.;
            break;
        case NLO:
            gammaDF1(0,0) = -388./9.;
            break;
        default:
            throw std::runtime_error("EvolDF1::GammaBB_e(): order not implemented");
    }
    return (gammaDF1);
}


//gslpp::matrix<double> EvolDF1::ToRescaleBasis(orders order, unsigned int n_u, unsigned int n_d) const
//{
//    
//    /* matrix entries for the anomalous dimension in the Chetyrkin, Misiak and Munz basis,
//       ref. hep-ph/9711280v1, hep-ph/0504194 */
//    
//    gslpp::matrix<double> mat(nops, 0.);
//    gslpp::matrix<double> mat1(nops, 0.);
//    unsigned int nf = n_u + n_d;
//    double z3 = gsl_sf_zeta_int(3);
//    
//    mat1(0,6) = - 13454./2187. + 44./2187.*nf;
//    mat1(1,6) = 20644./729. - 88./729.*nf;
//    mat1(2,6) = 119456./729. + 5440./729.*n_d -21776./729.*n_u;
//    mat1(3,6) = - 202990./2187. + 32./729.*n_d*n_d + n_d*(16888./2187. + 64./729.*n_u)
//                - 17132./2187.*n_u + 32./729.*n_u*n_u;
//    mat1(4,6) = 530240./243. + 300928./729.*n_d - 461120./729.*n_u;
//    mat1(5,6) = - 1112344./729. + 5432./729.*n_d*n_d + n_d*(419440./2187. - 
//                2744./729.*n_u) + 143392./2187.*n_u - 8176./729.*n_u*n_u;
//    
//    mat1(0,7) = 25759./5832. + 431./5832.*nf;
//    mat1(1,7) = 9733./486. - 917./972.*nf;
//    mat1(2,7) = 82873./243. - 3361./243.*nf;
//    mat1(3,7) = - 570773./2916. - 253./486.*n_d*n_d +n_d*(-40091./5832. - 
//                253./243.*n_u) - 40091./5832.*n_u - 253./486.*n_u*n_u;
//    mat1(4,7) = 838684./81. - 14.*n_d*n_d + n_d*(129074./243. - 28.*n_u) + 
//                129074./243.*n_u - 14.*n_u*n_u;
//    mat1(5,7) = - 923522./243. - 6031./486.*n_d*n_d + n_d*(-13247./1458. - 6031./243.*n_u)
//                -13247./1458.*n_u - 6031./486.*n_u*n_u;
//    
//    mat1(0,8) = - 2357278./19683. + 14440./6561.*n_d + 144688./6561.*n_u + 6976./243.*z3;
//    mat1(1,8) = - 200848./6561. - 23696./2187.*n_d + 30736./2187.*n_u - 3584./81.*z3;
//    mat1(2,8) = - 1524104./6561. - 176./27.*n_d*n_d + 352./27.*n_u*n_u +
//                n_d*(257564./2187. + 176./27.*n_u - 128./3.*z3) - 256./81.*z3 + 
//                n_u*(-382984./2187. + 256./3.*z3); 
//    mat1(3,8) = 1535926./19683. + 1984./2187.*n_d*n_d - 5792./2187.*n_u*n_u +
//                n_d*(-256901./6561. - 3808./2187.*n_u - 2720./81.*z3) - 
//                5056./243.*z3 + n_u*(34942./6561. + 1600./81.*z3);
//    mat1(4,8) = - 31433600./6561. - 2912./27.*n_d*n_d + 5824./27.*n_u*n_u +
//                n_d*(- 3786616./2187. + 2912./27.*n_u - 1280./3.*z3) -
//                4096./81.*z3 + n_u*(7525520./2187. + 2560./3.*z3);
//    mat1(5,8) = 48510784./19683. -51296./2187.*n_d*n_d + 54976./2187.*n_u*n_u +
//                n_u*(-11231648./6561. - 22016./81.*z3) + n_d*(340984./6561. + 
//                3680./2187.*n_u - 8192./81.*z3) - 80896./243.*z3;
//     
//    
//    switch(order){
//        case(NLO): 
//            mat = AnomalousDimension_M(NLO, n_u, n_d);
//            for (unsigned int i=0; i<6; i++){
//                for (unsigned int j=6; j<nops; j++){
//                    mat(i,j) = mat1(i,j);
//                }
//            }
//            for (unsigned int i=6; i<nops; i++){
//                for (unsigned int j=6; j<nops; j++){
//                    mat(i,j) = mat(i,j) + 2. * (i==j) * model.Beta1(nf);
//                }
//            }
//            return (mat);
//        case(LO):
//            mat = AnomalousDimension_M(LO, n_u, n_d);
//            for (unsigned int i=0; i<6; i++){
//                for (unsigned int j=6; j<nops; j++){
//                    mat(i,j) = AnomalousDimension_M(NLO, n_u, n_d)(i,j);
//                }
//            }
//            for (unsigned int i=6; i<nops; i++){
//                for (unsigned int j=6; j<nops; j++){
//                    mat(i,j) = mat(i,j) + 2. * (i==j) * model.Beta0(nf);
//                }
//            }
//            return (mat);
//        default:
//            throw std::runtime_error("change to rescaled operator basis: order not implemented"); 
//    }
//    
//}

//gslpp::matrix<double> EvolDF1::ToEffectiveBasis(gslpp::matrix<double> mat) const
//{
//    
//    gslpp::matrix<double> y(nops, 0.);
//    
//    y(0,0) = 1.;
//    y(1,1) = 1.;
//    y(2,2) = 1.;
//    y(3,3) = 1.;
//    y(4,4) = 1.;
//    y(5,5) = 1.;
//    y(6,6) = 1.;
//    y(7,7) = 1.;
//    y(8,8) = 1.;
//    y(9,9) = 1.;
//    y(10,10) = 1.;
//    y(11,11) = 1.;
//    y(12,12) = 1.;
//    
//    y(6,2) = -1./3.;
//    y(6,3) = -4./9.;
//    y(6,4) = -20./3.;
//    y(6,5) = -80./9.;
//    
//    y(7,2) = 1.;
//    y(7,3) = -1./6.;
//    y(7,4) = 20.;
//    y(7,5) = -10./3.;
//    
//    y(8,2) = 4./3.;
//    y(8,4) = 64./9.;
//    y(8,5) = 64./27.; // Add terms proportional to Log(mb/mub))
//    
//    return( (y.inverse()).transpose() * mat * y.transpose() );
//    
//}

gslpp::matrix<double>& EvolDF1::DF1Evol(double mu, double M, orders order, schemes scheme) 
{
    
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1::Df1Evol(): scheme " + out.str() + " not implemented "); 
    }
    
    double alsMZ = model.getAlsMz();
    double Mz = model.getMz();
    if(alsMZ == alsMZ_cache && Mz == Mz_cache) {
        if (mu == this->mu && M == this->M && scheme == this->scheme)
            return (*Evol(order));        
    }
    alsMZ_cache = alsMZ;
    Mz_cache = Mz;
        
    if (M < mu) {
        std::stringstream out;
        out << "M = " << M << " < mu = " << mu;
        throw out.str();
    }

    setScales(mu, M); // also assign evol to identity

    double m_down = mu;
    double m_up = model.AboveTh(m_down);
    double nf = model.Nf(m_down);
    
    while (m_up < M) {
        DF1Evol(m_down, m_up, nf, scheme);
        m_down = m_up;
        m_up = model.AboveTh(m_down);
        nf += 1.;
    }
    DF1Evol(m_down, M, nf, scheme);
    
    return (*Evol(order));
    
    }
    
 void EvolDF1::DF1Evol(double mu, double M, double nf, schemes scheme) 
 {

    gslpp::matrix<double> resLO(nops, 0.), resNLO(nops, 0.), resNNLO(nops, 0.);

    int L = 6 - (int) nf;
    double alsM = model.Als(M) / 4. / M_PI;
    double alsmu = model.Als(mu) / 4. / M_PI;
    
    double eta = alsM / alsmu;
    
    for (unsigned int k = 0; k < nops; k++) {
        double etap = pow(eta, a[L][k] / 2. / model.Beta0(nf));
        for (unsigned int i = 0; i < nops; i++){
            for (unsigned int j = 0; j < nops; j++) {
                resNNLO(i, j) += 0.;
                
                if(fabs(e(i).real() - e(j).real() + 2. * model.Beta0(nf))>0.000000000001)  {
                    resNLO(i, j) += c[L][i][j][k] * etap * alsmu;
                    resNLO(i, j) += d[L][i][j][k] * etap * alsM;
                }
                else{
                    resNLO(i, j) += - c[L][i][j][k] * etap * alsmu * log(eta);    
                }        
                resLO(i, j) += b[L][i][j][k] * etap;
                if (fabs(resLO(i, j)) < 1.e-12) {resLO(i, j) = 0.;}
                if (fabs(resNLO(i, j)) < 1.e-12) {resNLO(i, j) = 0.;}
            }
        }
    }
    
    switch(order) {
        case NNLO:
            *elem[NNLO] = 0.;
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
            break;
        case FULLNNLO:
        case FULLNLO:
        default:
            throw std::runtime_error("Error in EvolDF1bsg::Df1Evolbsg()");
    }   
  }
 


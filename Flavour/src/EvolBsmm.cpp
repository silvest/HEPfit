/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EvolBsmm.h"
#include "StandardModel.h"
#include <gsl/gsl_sf.h>
#include <algorithm>

EvolBsmm::EvolBsmm(unsigned int dim_i, schemes scheme, orders order, orders_qed order_qed, const StandardModel & model)
:   RGEvolutor(dim_i, scheme, order, order_qed), model(model), V(dim_i,0.), Vi(dim_i,0.),
    AA(dim_i,0.), BB(dim_i,0.), CC(dim_i,0.), DD(dim_i,0.), EE(dim_i,0.), FF(dim_i,0.),
    RR(dim_i,0.), e(dim_i,0.), vavi(0,0.), vbvi(0,0.), vcvi(0,0.), vdvi(0,0.),
        vevi(0,0.), vfvi(0,0.), vrvi(0,0.),vaevi(0,0.), vbbvi(0,0.), vbdvi(0,0.), vbevi(0,0.),
        vdbvi(0,0.), vdevi(0,0.), veavi(0,0.), vebvi(0,0.), vedvi(0,0.), veevi(0,0.), vbeevi(0,0.), 
        vebevi(0,0.), veebvi(0,0.), vbbevi(0,0.), vbebvi(0,0.), vebbvi(0,0.), dim(dim_i)
{
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int l = 0;
    unsigned int m = 0;
    unsigned int p = 0;
    unsigned int q = 0;
    int L = 1;
    double b0 = 0.;

  
    vavi.clear();
    vbvi.clear();
    vcvi.clear();
    vdvi.clear();
    vfvi.clear();
    vevi.clear();
    vrvi.clear();
    vaevi.clear();
    vbbvi.clear(); 
    vbdvi.clear();
    vbevi.clear();
    vdbvi.clear();
    vdevi.clear();
    veavi.clear();
    vebvi.clear();
    vedvi.clear();
    veevi.clear();
    vbeevi.clear();
    vebevi.clear(); 
    veebvi.clear();
    vbbevi.clear(); 
    vbebvi.clear(); 
    vebbvi.clear();
    
    
 
/* define L, nu, nd */
    if(L == 1){nd = 3; nu = 2;} 
    b0 = model.Beta0(6-L);
	    
    AnomalousDimension(10,nu,nd).transpose().eigensystem(V,e); 
   
    Vi = V.inverse();
	      
		/* magic numbers of U0 */
    for(unsigned int i = 0; i < dim; i++){
        a[L][i] = e(i).real()/2./b0;
	for (unsigned int j = 0; j < dim; j++){
            for (unsigned int k = 0; k < dim; k++){
		b[L][i][j][k] = V(i, k).real() * Vi(k, j).real();
            }
	}
    }
	    


    AA = BuiltB('A', nu, nd);  
    BB = BuiltB('B', nu, nd);
    CC = BuiltB('C', nu, nd);
    DD = BuiltB('D', nu, nd);
    EE = BuiltB('E', nu, nd);
    FF = BuiltB('F', nu, nd);
    RR = BuiltB('R', nu, nd);
        
    /* The products below only need the real parts of V, V^-1 and of the B matrices: copy them once to plain
     * arrays (element access through gslpp::matrix<gslpp::complex> dominated the construction time). The loop
     * order and the order of the factors are unchanged, so the stored lists are identical to the previous ones. */
    const unsigned int n = dim;
    auto realParts = [n](const gslpp::matrix<gslpp::complex>& matrix) {
        std::vector<double> r(n * n);
        for (unsigned int row = 0; row < n; row++)
            for (unsigned int col = 0; col < n; col++)
                r[row * n + col] = matrix(row, col).real();
        return r;
    };
    const std::vector<double> rV = realParts(V), rVi = realParts(Vi), rA = realParts(AA), rB = realParts(BB),
            rC = realParts(CC), rD = realParts(DD), rE = realParts(EE), rF = realParts(FF), rR = realParts(RR);

    const double cutoff = 0.000000001;
    auto store = [cutoff](std::vector<double>& list, std::initializer_list<double> indices, double value) {
        if (fabs(value) > cutoff) {
            list.insert(list.end(), indices);
            list.push_back(value);
        }
    };

    for (l = 0; l < n; l++) {
        for (i = 0; i < n; i++) {
            const double Vli = rV[l * n + i];
            for (j = 0; j < n; j++) {
                for (m = 0; m < n; m++) {
                    const double Vijm = rVi[j * n + m];
                    store(vavi, {(double) l, (double) i, (double) j, (double) m}, Vli * rA[i * n + j] * Vijm);
                    store(vbvi, {(double) l, (double) i, (double) j, (double) m}, Vli * rB[i * n + j] * Vijm);
                    store(vcvi, {(double) l, (double) i, (double) j, (double) m}, Vli * rC[i * n + j] * Vijm);
                    store(vdvi, {(double) l, (double) i, (double) j, (double) m}, Vli * rD[i * n + j] * Vijm);
                    store(vevi, {(double) l, (double) i, (double) j, (double) m}, Vli * rE[i * n + j] * Vijm);
                    store(vfvi, {(double) l, (double) i, (double) j, (double) m}, Vli * rF[i * n + j] * Vijm);
                    store(vrvi, {(double) l, (double) i, (double) j, (double) m}, Vli * rR[i * n + j] * Vijm);

                    for (p = 0; p < n; p++) {
                        const double Aip = rA[i * n + p], Bip = rB[i * n + p], Dip = rD[i * n + p], Eip = rE[i * n + p];
                        const double Epj = rE[p * n + j], Bpj = rB[p * n + j], Dpj = rD[p * n + j], Apj = rA[p * n + j];
                        store(vaevi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Aip * Epj * Vijm);
                        store(vbbvi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Bip * Bpj * Vijm);
                        store(vbdvi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Bip * Dpj * Vijm);
                        store(vbevi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Bip * Epj * Vijm);
                        store(vdbvi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Dip * Bpj * Vijm);
                        store(vdevi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Dip * Epj * Vijm);
                        store(veavi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Eip * Apj * Vijm);
                        store(vebvi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Eip * Bpj * Vijm);
                        store(vedvi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Eip * Dpj * Vijm);
                        store(veevi, {(double) l, (double) i, (double) p, (double) j, (double) m}, Vli * Eip * Epj * Vijm);

                        for (q = 0; q < n; q++) {
                            const double Epq = rE[p * n + q], Bpq = rB[p * n + q];
                            const double Eqj = rE[q * n + j], Bqj = rB[q * n + j];
                            store(vbeevi, {(double) l, (double) i, (double) p, (double) q, (double) j, (double) m}, Vli * Bip * Epq * Eqj * Vijm);
                            store(vebevi, {(double) l, (double) i, (double) p, (double) q, (double) j, (double) m}, Vli * Eip * Bpq * Eqj * Vijm);
                            store(veebvi, {(double) l, (double) i, (double) p, (double) q, (double) j, (double) m}, Vli * Eip * Epq * Bqj * Vijm);
                            store(vbbevi, {(double) l, (double) i, (double) p, (double) q, (double) j, (double) m}, Vli * Bip * Bpq * Eqj * Vijm);
                            store(vbebvi, {(double) l, (double) i, (double) p, (double) q, (double) j, (double) m}, Vli * Bip * Epq * Bqj * Vijm);
                            store(vebbvi, {(double) l, (double) i, (double) p, (double) q, (double) j, (double) m}, Vli * Eip * Bpq * Bqj * Vijm);
                        }
                    }
                }
            }
        }
    }

    for (std::vector<double>* list : {&vavi, &vbvi, &vcvi, &vdvi, &vevi, &vfvi, &vrvi})
        groupByArguments(*list, 2);
    for (std::vector<double>* list : {&vaevi, &vbbvi, &vbdvi, &vbevi, &vdbvi, &vdevi, &veavi, &vebvi, &vedvi, &veevi})
        groupByArguments(*list, 3);
    for (std::vector<double>* list : {&vbeevi, &vebevi, &veebvi, &vbbevi, &vbebvi, &vebbvi})
        groupByArguments(*list, 4);
}


template<typename Function>
void EvolBsmm::accumulate(const std::vector<double>& list, unsigned int nargs, std::vector<double>& target,
        Function function)
{
    // entries (l, arguments, m, coefficient), grouped by their arguments (see groupByArguments())
    const std::size_t stride = nargs + 3;
    const double* group = nullptr;
    double value = 0.;
    for (std::size_t k = 0; k + stride <= list.size(); k += stride) {
        const double* entry = &list[k];
        if (group == nullptr || !std::equal(entry + 1, entry + 1 + nargs, group + 1)) {
            value = function(entry + 1);
            group = entry;
        }
        target[(unsigned int) entry[0] * dim + (unsigned int) entry[nargs + 1]] += entry[nargs + 2] * value;
    }
}

void EvolBsmm::groupByArguments(std::vector<double>& list, unsigned int nargs)
{
    // Reorder the entries (l, arguments, m, coefficient) so that equal arguments are contiguous; the stable sort keeps
    // the original order within a group.
    const std::size_t stride = nargs + 3, size = list.size() / stride;
    std::vector<std::size_t> order(size);
    for (std::size_t k = 0; k < size; k++) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return std::lexicographical_compare(list.begin() + a * stride + 1, list.begin() + a * stride + 1 + nargs,
                list.begin() + b * stride + 1, list.begin() + b * stride + 1 + nargs);
    });
    std::vector<double> grouped;
    grouped.reserve(list.size());
    for (std::size_t k : order) grouped.insert(grouped.end(), list.begin() + k * stride, list.begin() + (k + 1) * stride);
    list.swap(grouped);
}

EvolBsmm::~EvolBsmm()
{}

gslpp::matrix<double> EvolBsmm::AnomalousDimension(int gam, unsigned int n_u, unsigned int n_d) const

{
	       
	   /* Delta F = 1 anomalous dimension in Misiak basis, 
	        ref.: ref. hep-ph/1311.1348v2, hep-ph/0512066 */   
	        
	/* gamma(row, coloumn) at the LO */
	    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    double zeta3 = 0;
        
    gslpp::matrix<double> gammaDF1(dim, 0.);
	
    zeta3 = gsl_sf_zeta_int(3);
        
    switch(gam){
	       
        case 10:
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            gammaDF1(0,0) = -4. ;
            gammaDF1(0,1) = 8./3. ;
            gammaDF1(0,3) = -2./9.;

            gammaDF1(1,0) = 12.;
            gammaDF1(1,3) = 4./3.;

            gammaDF1(2,3) = -52./3.;
            gammaDF1(2,5) = 2.;

            gammaDF1(3,2) = -40./9.;
            gammaDF1(3,3) = -100./9.;
            gammaDF1(3,4) = 4./9.;
            gammaDF1(3,5) = 5./6.;

            gammaDF1(4,3) = -256./3.;
            gammaDF1(4,5) = 20.;

            gammaDF1(5,2) = -256./9.;
            gammaDF1(5,3) = 56./9.;
            gammaDF1(5,4) = 40./9.;
            gammaDF1(5,5) = -2./3.;   
	  
	break;    
	case 20:
	   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            /* gamma(row, coloumn) at the NLO */

            gammaDF1(0,0) = -355./9.;
            gammaDF1(0,1) = -502./27.;
            gammaDF1(0,2) = -1412./243.;
            gammaDF1(0,3) = -1369./243.;
            gammaDF1(0,4) = 134./243.;
            gammaDF1(0,5) = -35./162.;

            gammaDF1(1,0) = -35./3.; 
            gammaDF1(1,1) = -28./3.;
            gammaDF1(1,2) = -416./81.;
            gammaDF1(1,3) = 1280./81.;
            gammaDF1(1,4) = 56./81.;
            gammaDF1(1,5) = 35./27.;

            gammaDF1(2,2) = -4468./81.;
            gammaDF1(2,3) = -31469./81.;
            gammaDF1(2,4) = 400./81.;
            gammaDF1(2,5) = 3373./108.;

            gammaDF1(3,2) = -8158./243.;       
            gammaDF1(3,3) = -59399./243.;
            gammaDF1(3,4) = 269./486.;
            gammaDF1(3,5) = 12899./648.;

            gammaDF1(4,2) = -251680./81.;
            gammaDF1(4,3) = -128648./81.;
            gammaDF1(4,4) = 23836./81.;
            gammaDF1(4,5) = 6106./27.;

            gammaDF1(5,2) = 58640./243.;
            gammaDF1(5,3) = -26348./243.;
            gammaDF1(5,4) = -14324./243.;
            gammaDF1(5,5) = -2551./162.;
	 
        break;

	case 30:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            /* gamma(row, coloumn) at the NNLO */

            gammaDF1(0,0) = -12773./18. +  zeta3 * 1472./3.;
            gammaDF1(0,1) = 745./9. - zeta3 *4288./9.;
            gammaDF1(0,2) = 63187./13122. - zeta3 * 1360./81.;
            gammaDF1(0,3) = -981796./6561. - zeta3 * 776./81.;
            gammaDF1(0,4) = -202663./52488. + zeta3 * 124./81.;
            gammaDF1(0,5) = -24973./69984. + zeta3 * 100./27.;

            gammaDF1(1,0) = 1177./2. - zeta3 * 2144.; 
            gammaDF1(1,1) = 306. - zeta3 * 224.;
            gammaDF1(1,2) = 110477./2187. + zeta3 * 2720./27.;
            gammaDF1(1,3) = 133529./8748. - zeta3 * 2768./27.;
            gammaDF1(1,4) = -42929./8748. - zeta3 * 248./27.;
            gammaDF1(1,5) = 354319./11664. - zeta3 * 110./9.;

            gammaDF1(2,2) = -3572528./2187. - zeta3 * 608./27.;
            gammaDF1(2,3) = -58158773./8748. + zeta3 * 61424./27.;
            gammaDF1(2,4) = 552601./4374. - zeta3 * 496./27.;
            gammaDF1(2,5) = 6989171./11664. - zeta3 * 2821./9.;

            gammaDF1(3,2) = -1651004./6561. + zeta3 * 88720./81.;       
            gammaDF1(3,3) = -155405353./52488 + zeta3 * 54272./81.;
            gammaDF1(3,4) = 1174159./52488. - zeta3 * 9274./81.;
            gammaDF1(3,5) = 10278809./34992. - zeta3 * 3100./27.;

            gammaDF1(4,2) = -147978032./2187. + zeta3 * 87040./27.;
            gammaDF1(4,3) = -168491372./2187. + zeta3 * 324416./27.;
            gammaDF1(4,4) = 11213042./2187. - zeta3 * 13984./27.;
            gammaDF1(4,5) = 17850329./2916. - zeta3 * 31420./9.;

            gammaDF1(5,2) = 136797922./6561. + zeta3 * 721408./81.;
            gammaDF1(5,3) = -72614473./13122. - zeta3 * 166432./81.;
            gammaDF1(5,4) = -9288181./6561. - zeta3 * 95032./81.;
            gammaDF1(5,5) = -16664027./17496. - zeta3 * 7552./27.;
	
        break;	
	case 01:

            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }    

            gammaDF1(0,0) = -8./3.;
            gammaDF1(0,6) = -32./27.;

            gammaDF1(1,1) = -8./3.;
            gammaDF1(1,6) = -8./9.;

            gammaDF1(2,6) = -16./9.;

            gammaDF1(3,6) = 32./27.;

            gammaDF1(4,6) = -112./9.;

            gammaDF1(5,6) = 512./27.;  

            gammaDF1(6,6) = 8.;
            gammaDF1(6,7) = -4.; 

            gammaDF1(7,6) = -4.; 
  	
        break;  
	case 11:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }


            gammaDF1(0,0) = 169./9.;
            gammaDF1(0,1) = 100./27.;
            gammaDF1(0,3) = 254./729.;
            gammaDF1(0,6) = -2272./729.;

            gammaDF1(1,0) = 50./3.; 
            gammaDF1(1,1) = -8./3.;
            gammaDF1(1,3) = 1076./243.;
            gammaDF1(1,6) = 1952./243.;

            gammaDF1(2,3) = 11116./243.;
            gammaDF1(2,5) = -14./3.;
            gammaDF1(2,6) = -6752./243.;

            gammaDF1(3,2) = 280./27.;       
            gammaDF1(3,3) = 18763./729.;
            gammaDF1(3,4) = -28./27.;
            gammaDF1(3,5) = -35./18.;
            gammaDF1(3,6) = -2192./729.;  

            gammaDF1(4,3) = 111136./243.;
            gammaDF1(4,5) = -140./3.;
            gammaDF1(4,6) = -84032./243.;

            gammaDF1(5,2) = 2944./27.;
            gammaDF1(5,3) = 193312./729.;
            gammaDF1(5,4) = -280./27.;
            gammaDF1(5,5) = -175./9.;
            gammaDF1(5,6) = -37856./729.;

            gammaDF1(6,7) = 16.;

            gammaDF1(7,6) = 16.; 

        break;
	case 02:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                    throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

             gammaDF1(0,6)= -11680./2187.;
             gammaDF1(0,7)= -416./81.;

             gammaDF1(1,6)= -2920./729.;
             gammaDF1(1,7)= -104./27.;

             gammaDF1(2,6)= -39752./729.;
             gammaDF1(2,7)= -136./27.;

             gammaDF1(3,6)= 1024./2187.;
             gammaDF1(3,7)= -448./81.;

             gammaDF1(4,6)= -381344./729.;
             gammaDF1(4,7)= -15616./27.;

             gammaDF1(5,6)= 24832./2187.;
             gammaDF1(5,7)= -7936./81.;


        break;
	case 21:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            gammaDF1(0,6)= -1359190./19683. + zeta3 * 6976./243.;

            gammaDF1(1,6)= -229696./6561 - zeta3 * 3584./81.;

            gammaDF1(2,6)= -1290092./6561 + zeta3 * 3200./81.;

            gammaDF1(3,6)= -819971./19683. - zeta3 * 19936./243;

            gammaDF1(4,6)= -16821944./6561 + zeta3 * 30464./81.;

            gammaDF1(5,6)= -17787368./19683. - zeta3 * 286720./243.;


        break;

	default:
            throw std::runtime_error("EvolBsmm::AnomalousDimension(): order not implemented"); 
    }
    return (gammaDF1);
}


gslpp::matrix<double> EvolBsmm::BuiltB(char letter, unsigned int n_u, unsigned int n_d)

{

    unsigned int nf = 5; //all the class works for nf = 5

    double B00S = model.Beta0(nf), B10S = model.Beta1(nf), B20S = model.Beta2(nf), 
            B01S = -22./9., B11S = -308./27.; 

    double B00E = 80./9., B01E = 176./9.;

    double b1 = B10S/(2. * B00S * B00S), b2 = B20S/(4. * B00S * B00S * B00S) - b1 * b1 , 
            b3 = B01S/(2. * B00S * B00E), b4 = B11S /(4. * B00S * B00S * B00E) - 2 * b1 * b3, 
            b5 = B01E/(2. * B00S * B00E) - b1;



    gslpp::matrix<double> B(dim, 0.);
    gslpp::matrix<double> W10T(dim, 0.);
    gslpp::matrix<double> W20T(dim, 0.);
    gslpp::matrix<double> W30T(dim, 0.);
    gslpp::matrix<double> W01T(dim, 0.);
    gslpp::matrix<double> W02T(dim, 0.);
    gslpp::matrix<double> W11T(dim, 0.);
    gslpp::matrix<double> W21T(dim, 0.);


    W10T = (AnomalousDimension(10, n_u, n_d).transpose())/2./B00S;
    W20T = (AnomalousDimension(20, n_u, n_d).transpose())/4./B00S/B00S;
    W30T = (AnomalousDimension(30, n_u, n_d).transpose())/8./B00S/B00S/B00S;
    W01T = (AnomalousDimension(01, n_u, n_d).transpose())/2./B00E;
    W02T = (AnomalousDimension(02, n_u, n_d).transpose())/4./B00E/B00E;
    W11T = (AnomalousDimension(11, n_u, n_d).transpose())/4./B00S/B00E;
    W21T = (AnomalousDimension(21, n_u, n_d).transpose())/8./B00S/B00S/B00E;



    switch(letter){

        case 'A':

            B = Vi.real() * (W30T - b1 * W20T - b2 * W10T) * V.real();

        break;    
        case 'B':

            B = Vi.real() * (W20T - b1 * W10T) * V.real();

        break;
        case 'C':

            B = Vi.real() * (W21T - b1 * W11T - b2 * W01T - b3 * W20T - b4 * W10T) * V.real();

        break;
        case 'D':

            B = Vi.real() * (W11T - b1 * W01T - b3 * W10T) * V.real();

        break;
        case 'E':

            B = Vi.real() * (W01T) * V.real();

        break;
        case 'F':

            B = Vi.real() * (W02T + W11T - (b1 + b3) * W01T - b3 * W10T) * V.real();

        break;
        case 'R':

            B = b5 * Vi.real() * (W01T) * V.real();

        break;
        default:
            throw std::runtime_error("EvolBsmm::BuiltB(): order not implemented"); 
    }
    return (B);
}



gslpp::matrix< double > & EvolBsmm::Df1Evol(double mu, double M, orders order, orders_qed order_qed, schemes scheme) 

{
    switch (scheme) {             /*  complete this method */
        case NDR:
        break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1nlep::Df1Evol(): scheme " + out.str()  + " not implemented ");
    }

    // the evolution also depends on alpha_s(M_Z), alpha(M_Z) and M_Z: include them in the cache check
    // (it used to be keyed on the scales only, so the first evaluated point fixed the evolution for all others)
    bool same = (mu == this->mu && M == this->M && scheme == this->scheme && model.getAlsMz() == cacheAlsMz
                 && model.alphaMz() == cacheAleMz && model.getMz() == cacheMz);

    if (same && order_qed == NO_QED)
        return (*Evol(order));

    if (same && order_qed != NO_QED)
        return (*Evol(order_qed));

    cacheAlsMz = model.getAlsMz();
    cacheAleMz = model.alphaMz();
    cacheMz = model.getMz();


    if (M < mu) {
        std::stringstream out;
        out << "M = " << M << " < mu = " << mu;
        throw out.str();
    }

    setScales(mu, M); // also assign evol to identity
    if (M != mu) {
        double m_down = mu;
        //double m_up = model.AboveTh(m_down);
        //double nf = model.Nf(m_down);
        double nf = 5; //all the process in implemented for nf = 5
        alsM = alphatilde_s(M);
        alsmu = alphatilde_s(mu);

        eta = alsM / alsmu;
        logeta = log(eta);
        for (unsigned int j = 0; j < dim; j++)
            for (int shift = -8; shift <= 8; shift++)
                etaPowTable[j][shift + 8] = pow(eta, a[1][j] + shift); // a[L] with L = 6 - nf = 1

        /*while (m_up < M) {                         //there is no thresholds
                Df1Evol(m_down, m_up, nf, scheme);
                //Df1threshold_nlep(m_up, nf+1.);
                m_down = m_up;
                m_up = model.AboveTh(m_down);
                nf += 1.;
        } */

        Df1Evol(m_down, M, nf, scheme);
    }

    if(order_qed != NO_QED) return (*Evol(order_qed));
    else return (*Evol(order)); 
    
}	

void EvolBsmm::Df1Evol(double mu, double M, double nf, schemes scheme)

{ 
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;


    int L = 6 - (int) nf;

    double B00S = model.Beta0(nf);// B10S = model.Beta1(nf); /*inizializza i B*/

    double B00E = 80./9.;// B01E = 176./9.; /*inizializza i B*/

    double omega = 2. * B00S , lambda = B00E  /B00S ;    /*   E ale dipende da mu?*/



    /* Each list entry is (l, arguments, m, coefficient) and adds coefficient * f(arguments) to the (l, m) element
     * of one of the matrices below. The lists are grouped by their arguments in the constructor, so that every
     * F, R, G and H value is computed once per group instead of once per (l, m) pair. */
    const unsigned int n = dim;
    std::vector<double> us(n * n, 0.), us2(n * n, 0.), ue(n * n, 0.), ueos(n * n, 0.), ue2os(n * n, 0.),
            ue2os2(n * n, 0.), ues(n * n, 0.);

    accumulate(vbeevi, 4, ue2os, [&](const double* k) { return H(k[0], k[1], k[2], k[3], 2, 4, 4, mu, M, nf); });
    accumulate(vebevi, 4, ue2os, [&](const double* k) { return H(k[0], k[1], k[2], k[3], 4, 2, 4, mu, M, nf); });
    accumulate(veebvi, 4, ue2os, [&](const double* k) { return H(k[0], k[1], k[2], k[3], 4, 4, 2, mu, M, nf); });
    accumulate(vbbevi, 4, ues, [&](const double* k) { return H(k[0], k[1], k[2], k[3], 2, 2, 4, mu, M, nf); });
    accumulate(vbebvi, 4, ues, [&](const double* k) { return H(k[0], k[1], k[2], k[3], 2, 4, 2, mu, M, nf); });
    accumulate(vebbvi, 4, ues, [&](const double* k) { return H(k[0], k[1], k[2], k[3], 4, 2, 2, mu, M, nf); });

    accumulate(vaevi, 3, ues, [&](const double* k) { return G(k[0], k[1], k[2], 1, 4, mu, M, nf); });
    accumulate(vbbvi, 3, us2, [&](const double* k) { return G(k[0], k[1], k[2], 2, 2, mu, M, nf); });
    accumulate(vbdvi, 3, ues, [&](const double* k) { return G(k[0], k[1], k[2], 2, 3, mu, M, nf); });
    accumulate(vbevi, 3, ue, [&](const double* k) { return G(k[0], k[1], k[2], 2, 4, mu, M, nf); });
    accumulate(vbevi, 3, ue2os, [&](const double* k) {
        return -G(k[0], k[1], k[2], 2, 4, mu, M, nf) + G(k[0], k[1], k[2], 2, 5, mu, M, nf); });
    accumulate(vdbvi, 3, ues, [&](const double* k) { return G(k[0], k[1], k[2], 3, 2, mu, M, nf); });
    accumulate(vdevi, 3, ue2os, [&](const double* k) { return G(k[0], k[1], k[2], 3, 4, mu, M, nf); });
    accumulate(veavi, 3, ues, [&](const double* k) { return G(k[0], k[1], k[2], 4, 1, mu, M, nf); });
    accumulate(vebvi, 3, ue, [&](const double* k) { return G(k[0], k[1], k[2], 4, 2, mu, M, nf); });
    accumulate(vebvi, 3, ue2os, [&](const double* k) {
        return -G(k[0], k[1], k[2], 4, 2, mu, M, nf) + G(k[0], k[1], k[2], 5, 2, mu, M, nf); });
    accumulate(vedvi, 3, ue2os, [&](const double* k) { return G(k[0], k[1], k[2], 4, 3, mu, M, nf); });
    accumulate(veevi, 3, ue2os2, [&](const double* k) { return G(k[0], k[1], k[2], 4, 4, mu, M, nf); });

    accumulate(vavi, 2, us2, [&](const double* k) { return F(k[0], k[1], 1, mu, M, nf); });
    accumulate(vbvi, 2, us, [&](const double* k) { return F(k[0], k[1], 2, mu, M, nf); });
    accumulate(vcvi, 2, ues, [&](const double* k) { return F(k[0], k[1], 2, mu, M, nf); });
    accumulate(vdvi, 2, ue, [&](const double* k) { return F(k[0], k[1], 3, mu, M, nf); });
    accumulate(vdvi, 2, ue2os, [&](const double* k) { return -F(k[0], k[1], 3, mu, M, nf); });
    accumulate(vevi, 2, ueos, [&](const double* k) { return F(k[0], k[1], 4, mu, M, nf); });
    accumulate(vevi, 2, ue2os2, [&](const double* k) {
        return F(k[0], k[1], 5, mu, M, nf) - F(k[0], k[1], 4, mu, M, nf); });
    accumulate(vfvi, 2, ue2os, [&](const double* k) { return F(k[0], k[1], 4, mu, M, nf); });
    accumulate(vrvi, 2, ue2os, [&](const double* k) { return R(k[0], k[1], 4, mu, M, nf); });

    /* Combine with the evolution matrices of the previous step (*elem[...], all taken before any update, as in the
     * fall-through switches below) using plain arrays: the 8x8 products through gslpp element access dominated the
     * cost otherwise. Sums are taken left to right, as before. */
    typedef std::vector<double> Dense;
    auto scaled = [n](const Dense& values, double factor) {
        Dense r(n * n);
        for (unsigned int k = 0; k < n * n; k++) r[k] = factor * values[k];
        return r;
    };
    auto dense = [n](const gslpp::matrix<double>& matrix) {
        Dense r(n * n);
        for (unsigned int row = 0; row < n; row++)
            for (unsigned int col = 0; col < n; col++)
                r[row * n + col] = matrix(row, col);
        return r;
    };
    auto mul = [n](const Dense& a, const Dense& b) {
        Dense c(n * n);
        for (unsigned int row = 0; row < n; row++)
            for (unsigned int col = 0; col < n; col++) {
                double sum = 0.;
                for (unsigned int k = 0; k < n; k++) sum += a[row * n + k] * b[k * n + col];
                c[row * n + col] = sum;
            }
        return c;
    };
    auto add = [n](std::initializer_list<Dense> terms) {
        Dense r(n * n, 0.);
        for (const Dense& term : terms)
            for (unsigned int k = 0; k < n * n; k++) r[k] += term[k];
        return r;
    };
    auto assign = [n](gslpp::matrix<double>& matrix, const Dense& values) {
        for (unsigned int row = 0; row < n; row++)
            for (unsigned int col = 0; col < n; col++)
                matrix(row, col) = values[row * n + col];
    };

    Dense rLO(n * n, 0.), rUe2(n * n, 0.);
    for (k = 0; k < n; k++) {
        const double etap = etaPow(k, 0);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                rLO[i * n + j] += b[L][i][j][k] * etap;
    }
    for (i = 0; i < n; i++) rUe2[i * n + i] = 1.;
    const Dense rUs = scaled(us, omega), rUs2 = scaled(us2, omega * omega), rUeos = scaled(ueos, lambda),
            rUe = scaled(ue, lambda * omega), rUe2os2 = scaled(ue2os2, lambda * lambda),
            rUes = scaled(ues, omega * omega * lambda), rUe2os = scaled(ue2os, omega * lambda * lambda);

    const Dense eLO = dense(*elem[LO]), eNLO = dense(*elem[NLO]), eNNLO = dense(*elem[NNLO]),
            e11 = dense(*elem[NLO_QED11]), e21 = dense(*elem[NLO_QED21]), e22 = dense(*elem[NLO_QED22]);

    switch(order_qed) {
        case NLO_QED22:
            assign(*elem[NLO_QED22], add({mul(e22, rLO), mul(e11, rUe), mul(eLO, rUe2), mul(e21, rUeos),
                    mul(eNNLO, rUe2os2), mul(eNLO, rUe2os)}));
        case NLO_QED12:
            assign(*elem[NLO_QED12], add({mul(e11, rUeos), mul(eNLO, rUe2os2), mul(eLO, rUe2os)}));
        case NLO_QED21:
            assign(*elem[NLO_QED21], add({mul(e21, rLO), mul(e11, rUs), mul(eNLO, rUe), mul(eLO, rUes), mul(eNNLO, rUeos)}));
        case NLO_QED02:
            assign(*elem[NLO_QED02], mul(eLO, rUe2os2));
        case NLO_QED11:
            assign(*elem[NLO_QED11], add({mul(e11, rLO), mul(eLO, rUe), mul(eNLO, rUeos)}));
        case LO_QED:
            assign(*elem[LO_QED], mul(eLO, rUeos));
            break;
        default:
            throw std::runtime_error("Error in EvolBsmm::Df1Evol");
    }

    switch(order) {
        case NNLO:
            assign(*elem[NNLO], add({mul(eLO, rUs2), mul(eNNLO, rLO), mul(eNLO, rUs)}));
        case NLO:
            assign(*elem[NLO], add({mul(eLO, rUs), mul(eNLO, rLO)}));
        case LO:
            assign(*elem[LO], mul(eLO, rLO));
            break;
        default:
            throw std::runtime_error("Error in EvolBsmm::Df1Evol");
    }
}

double EvolBsmm::F(unsigned int i, unsigned int j, int x, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;

    double etai = etaPow(i, 0);
    double etajx = etaPow(j, x - 3);
        
    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][j] + x - 3. - a[L][i]) < cut ) value = 0;
    else value = 1;	
    
    switch(value) {
        case 0:
            result = etai * logeta;
        break;
        case 1:
            result = (etajx - etai)/(a[L][j] + x - 3. - a[L][i]);
        break;
    }
    return (result);
}


double EvolBsmm::R(unsigned int i, unsigned int j, int x, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;

    double etai = etaPow(i, 0);
    double etajx = etaPow(j, x - 3);
        
    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][j] + x - 3. - a[L][i]) < cut ) value = 0;
    else value = 1;	
    
    switch(value) {
        case 0:
            result = etai * logeta * logeta * 0.5;
        break;
        case 1:
            result = (etajx * logeta - ((etajx - etai)/(a[L][j] + x - 3. - a[L][i])))/(a[L][j] + x - 3. - a[L][i]);
        break;
    }
return (result);
}

double EvolBsmm::G(unsigned int i, unsigned int p, unsigned int j, int x, int y, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;
    
    double etai = etaPow(i, 0);
    double etapx = etaPow(p, x - 3);
    double etajxy = etaPow(j, x - 3 + y - 3);
        
    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][j] + y - 3. - a[L][p]) < cut && fabs(a[L][p] + x - 3. - a[L][i]) < cut ) value = 0;
    else if(fabs(a[L][j] + y - 3. - a[L][p]) < cut && fabs(a[L][p] + x - 3. - a[L][i]) > cut) value = 1;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) < cut &&  fabs(a[L][p] + x - 3. - a[L][i]) < cut ) value = 2;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) > cut &&  fabs(a[L][p] + x - 3. - a[L][i]) > cut ) value = 3;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) < cut &&  fabs(a[L][p] + x - 3. - a[L][i]) > cut ) value = 4;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) > cut &&  fabs(a[L][p] + x - 3. - a[L][i]) < cut ) value = 5;	

    switch(value) {
        case 0:
            result = etai * logeta * logeta * 0.5;
        break;
        case 1:
            result = ((etapx * logeta - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][p] + x - 3. - a[L][i]));
        break;
        case 2:
            result = (etai * logeta - etai * logeta)/(a[L][j] + y - 3. - a[L][p]);
        break;
        case 3:
            result = (((etajxy - etai)/(a[L][j] + x - 3. + y - 3. - a[L][i])) - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][j] + y - 3. - a[L][p]);
        break;
        case 4:
            result = (etai * logeta - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][j] + y - 3. - a[L][p]);
        break;
        case 5:
            result = (((etajxy - etai)/(a[L][j] + x - 3. + y - 3. - a[L][i])) - etai * logeta)/(a[L][j] + y - 3. - a[L][p]);
        break;
        }
    return (result);
}

double EvolBsmm::H(unsigned int i, unsigned int p, unsigned int q, unsigned int j, int x, int y, int z, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;

    double etai = etaPow(i, 0);
    double etapx = etaPow(p, x - 3);

    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][p] + x - 3. - a[L][i]) < cut && fabs(a[L][q] + y - 3. - a[L][p]) < cut && fabs(a[L][j] + z - 3. - a[L][q]) < cut) value = 0;
    else if(fabs(a[L][p] + x - 3. - a[L][i]) > cut  && fabs(a[L][q] + y - 3. - a[L][p]) < cut && fabs(a[L][j] + z - 3. - a[L][q]) < cut) value = 1;
    /* H is the third divided difference of eta^t at t = a_i, b_p = a_p + x - 3, c_q = a_q + x - 3 + y - 3 and
     * c_j = a_j + x - 3 + y - 3 + z - 3 (G and F are the second and first ones); degenerate points give its limits.
     * c_j = c_q != b_p: the limit of case 6 is (R(i,q,x+y-3) - G(i,p,q,x,y))/(c_q - b_p), where R and G also handle
     * c_q = a_i and b_p = a_i (as in EvolDF1::f_h) */
    else if(fabs(a[L][j] + z - 3. - a[L][q]) < cut) value = 2;
    else if(fabs(a[L][j] + z - 3. - a[L][q]) > cut) value = 6;

    switch(value) {
        case 0:
            result = etai * logeta * logeta * logeta /6.;
        break;
        case 1:
            result = (0.5 * etapx * logeta * logeta - ((etapx * logeta 
                    - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][p] + x - 3. - a[L][i])))/(a[L][p] + x - 3. - a[L][i]);
        break;
        case 2:
            result = (R(i, q, x + y - 3, mu, M, nf) - G(i, p, q, x, y, mu, M, nf))/(a[L][q] + y - 3. - a[L][p]);
        break;
        case 6:
            result = (G(i, p, j, x , y + z - 3., mu, M, nf) - G(i, p, q, x , y, mu, M, nf))/(a[L][j] + z - 3. - a[L][q]);		
        break;
    }
return (result);
}

double EvolBsmm::alphatilde_e(double mu)

{     // also the running is only for nf = 5
    
    //double mu_0 = 91.1876;
    double mu_0 = model.getMz();
    //double alphatilde_e = 1./(127.751 * 4. * M_PI); // alpha_e at mu_0 = 91.1876 Gev
    double alphatilde_e = model.alphaMz()/4./M_PI;
    //double alphatilde_s = 0.1184/(4.* M_PI); // alpha_s at mu_0 = 91.1876 Gev
    double alphatilde_s = model.getAlsMz()/4./M_PI;
    unsigned int nf = 5;

    double B00S = model.Beta0(nf), B10S = model.Beta1(nf);//B20S = model.Beta2(nf), 
    double B01S = -22./9.; //B11S = -308./27.; 

    double B00E = 80./9., B01E = 176./9., B10E = 464./27.;

    //double b1 = B10S/(2. * B00S * B00S); //b2 = B20S/(4. * B00S * B00S * B00S) - b1 * b1 , 
    //double b3 = B01S/(2. * B00S * B00E);// b4 = B11S /(4. * B00S * B00S * B00E) - 2 * b1 * b3, 
            //b5 = B01E/(2. * B00S * B00E) - b1;

    double vs= 1. + 2. * B00S * alphatilde_s * log(mu/ mu_0);
    double ve= 1. - 2. * B00E * alphatilde_e * log(mu/ mu_0);
    //double ps= B00S * alphatilde_s /(B00S * alphatilde_s + B00E * alphatilde_e);
    double pe= B00E * alphatilde_e /(B00S * alphatilde_s + B00E * alphatilde_e);

    double logve = log(ve);
    double logvs = log(vs);
    double logeos = log(ve/vs);
    double asovs = alphatilde_s/vs;
    double aeove = alphatilde_e/ve;

    double result = 0;

    result = aeove - pow(aeove, 2) * (logve * B10E/ B00E - logvs * B01E/B00S) 
            + pow(aeove, 2) * (asovs) * ((logvs - vs + 1.) * B01E * B10S/(B00S * B00S)
            + logve * vs * pe * B01E * B10E/(B00E * B00E) +(logeos * vs * pe - logve) * B01E * B01S/(B00S * B00E));
    return (result);
}

double EvolBsmm::alphatilde_s(double mu)

{  // also the running is only for nf = 5
    
    //double mu_0 = 91.1876;
    double mu_0 = model.getMz();
    //double alphatilde_e = 1./(127.751 * 4. * M_PI); // alpha_e at mu_0 = 91.1876 Gev
    double alphatilde_e = model.alphaMz()/4./M_PI;
    //double alphatilde_s = 0.1184/(4.* M_PI); // alpha_s at mu_0 = 91.1876 Gev
    double alphatilde_s = model.getAlsMz()/4./M_PI;
    unsigned int nf = 5;

    double B00S = model.Beta0(nf), B10S = model.Beta1(nf), B20S = model.Beta2(nf), B30S = gsl_sf_zeta_int(3) * 352864./81. - 598391./1458,
            B01S = -22./9., B11S = -308./27., B02S = 4945./243.; 

    double B00E = 80./9., B01E = 176./9., B10E = 464./27.; 

    //double B00S2 = B00S * B00S;
    double B10soB00s = B10S / B00S;
    double B01soB00e = B01S/B00E;

    //double b1 = B10soB00s/(2. * B00S), b2 = B20S/(4. * B00S2 * B00S) - b1 * b1 , 
    //        b3 = B01soB00e/(2. * B00S ), b4 = B11S /(4. * B00S2 * B00E) - 2 * b1 * b3, 
    //        b5 = B01E/(2. * B00S * B00E) - b1;

    double vs= 1. + 2. * B00S * alphatilde_s * log(mu/ mu_0);
    double ve= 1. - 2. * B00E * alphatilde_e * log(mu/ mu_0);
    double ps= B00S * alphatilde_s /(B00S * alphatilde_s + B00E * alphatilde_e);
    //double pe= B00E * alphatilde_e /(B00S * alphatilde_s + B00E * alphatilde_e);

    double logve = log(ve);
    double logvs = log(vs);
    double logeos = log(ve/vs);
    double logsoe = log(vs/ve);
    double asovs = alphatilde_s/vs;
    double aeove = alphatilde_e/ve;

    double result = 0;

    result = asovs - pow(asovs, 2) * (logvs * B10soB00s - logve * B01soB00e) 
            +  pow(asovs, 3) * ((1. - vs) * B20S / B00S + B10soB00s * B10soB00s * (logvs * logvs - logvs
            + vs - 1.) + B01soB00e * B01soB00e * logve * logve + (-2. * logvs * logve 
            + ps * ve * logve) * B01S * B10S/(B00E * B00S)) 
            +  pow(asovs, 4) * (0.5 * B30S *(1. - vs * vs)/ B00S + ((2. * vs - 3.) * logvs + vs * vs 
            - vs) * B20S * B10soB00s /(B00S) + B10soB00s * B10soB00s * B10soB00s * (- pow(logvs,3) 
            + 5. * pow(logvs,2) / 2. + 2. * (1. - vs) * logvs - (vs - 1.) * (vs - 1.)* 0.5)) 
            + pow(asovs, 2) * (aeove) * ((ve - 1.) * B02S / B00E 
            + ps * ve * logeos * B11S /B00S +(logve - ve + 1.) * B01soB00e * B10E/(B00E) 
            + logvs * ps * B01S * B10soB00s/(B00S) +(logsoe * ve * ps - logvs) * B01soB00e * B01E/( B00S));
    return (result);
}
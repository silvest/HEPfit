/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "StandardModel.h"
using namespace std;

void setSMparameters(StandardModel& SM_i) {
    std::map<std::string, double> Parameters;
    // 18+5 parameters defined in StandardModel
    Parameters["GF"] = 1.16637E-5;
    Parameters["mneutrino_1"] = 0.0;
    Parameters["mneutrino_2"] = 0.0;
    Parameters["mneutrino_3"] = 0.0;
    Parameters["melectron"] = 0.54857990943e-3;
    Parameters["mmu"] = 0.1134289256;
    Parameters["mtau"] = 1.77682;
    Parameters["lambda"] = 0.2253;
    Parameters["A"] = 0.808;
    Parameters["rhob"] = 0.132;
    Parameters["etab"] = 0.341;
    Parameters["ale"] = 1.0/137.035999679;
    Parameters["dAle5Mz"] = 0.02758;
    Parameters["mHl"] = 130.0;
    Parameters["muw"] = 80.0;
    //
    Parameters["phiEpsK"] = 0.0;
    Parameters["DeltaMK"] = 0.0;
    Parameters["KbarEpsK"] = 0.0;
    Parameters["Dmk"] = 0.0;
    Parameters["SM_M12D" ] = 0.0;
    
    // 26+16+1 parameters defined in QCD    
    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    Parameters["mtop"] = 174.0;
    Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 175.0;
    Parameters["mub"] = 4.7;
    Parameters["muc"] = 1.5;
    Parameters["MBd"] = 0.0;
    Parameters["MBs"] = 0.0;
    Parameters["MBp"] = 0.0;
    Parameters["MK0"] = 0.0;
    Parameters["MKp"] = 0.0;
    Parameters["FBs"] = 0.0;
    Parameters["FBsoFBd"] = 0.0;
    Parameters["BBsoBBd"] = 0.0;
    Parameters["BBs1"] = 0.0;
    Parameters["BBs2"] = 0.0;  
    Parameters["BBs3"] = 0.0;
    Parameters["BBs4"] = 0.0;
    Parameters["BBs5"] = 0.0;
    Parameters["BBsscale"] = 0.0;
    Parameters["BBsscheme"] = 0.0;
    //
    Parameters["MD"] = 0.0;
    Parameters["FD"] = 0.0;
    Parameters["BD1"] = 0.0;
    Parameters["BD2"] = 0.0;
    Parameters["BD3"] = 0.0;
    Parameters["BD4"] = 0.0;
    Parameters["BD5"] = 0.0;
    Parameters["BDscale"] = 0.0;
    Parameters["BDscheme"] = 0.0;
    Parameters["BK1"] = 0.0;
    Parameters["BK2"] = 0.0;
    Parameters["BK3"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK5"] = 0.0;
    Parameters["BKscale"] = 0.0;
    Parameters["BKscheme"] = 0.0; 
    Parameters["FK"] = 0.0; 
    
    /** Test for alpha_lep **/
    //Parameters["melectron"] = 0.00051099907;
    //Parameters["mmu"] = 0.105658389;
    //Parameters["mtau"] = 1.777;    
    //Parameters["AlsMz"] = 1.0/137.0359895;
    //Parameters["Mz"] = 91.187;    
    
    /** To make comparisons with ZFITTER codes **/
    Parameters["GF"] = 1.16637E-5;  // for GFER=2
    Parameters["melectron"] = 0.51099907e-3;
    Parameters["mmu"] = 0.105658389;
    Parameters["mtau"] = 1.77705;
    Parameters["mup"] = 0.062;
    Parameters["mdown"] = 0.083;
    Parameters["mcharm"] = 1.50; // In ZFitter, 1.5 is the pole mass.
    Parameters["mstrange"] = 0.215;
    Parameters["mbottom"] = 4.70; // In ZFitter, 4.7 is the pole mass.
    Parameters["ale"] = 1.0/137.0359895;

    /* TEST for Table 6.1 in hep-ph/0507146*/
    /* flags: AMT4=6, ALEM=2 */
    Parameters["Mz"] = 91.1875;
    Parameters["mtop"] = 175.0;    
    Parameters["mHl"] = 150.0;
    Parameters["AlsMz"] = 0.118;
    Parameters["dAle5Mz"] = 0.02758;    

    /* to make mb(Mz) and mc(Mz) similar to ZFitter ones */
    /* mcMz = 0.56381685, mbMz = 2.8194352 */
    /* muMz = 0.062, mdMz = 0.083, msMz = 0.215 */
    Parameters["mbottom"] = 4.122;
    Parameters["mub"] = 4.122;    
    Parameters["mcharm"] = 1.171;
    Parameters["muc"] = 1.171;   
    
    SM_i.Init(Parameters);
}

int main(int argc, char** argv) {
    
    try {    
        StandardModel* mySM;
        mySM = new StandardModel(true);
        mySM->InitializeModel();
        setSMparameters(*mySM);

        cout.precision(8);
        cout.setf(ios::floatfield);
        cout.setf(ios::left);

        cout << "[Inputs]" << endl
             << "   M_Z = " << mySM->getMz()
             << "   alpha_s(M_Z) = " << mySM->getAlsMz() << endl
             << "   M_t = " << mySM->getMtpole() 
             << "   m_b(m_b) = " << mySM->getQuarks(mySM->BOTTOM).getMass()
             << "   m_c(m_c) = " << mySM->getQuarks(mySM->CHARM).getMass() << endl
             << "   mut = " << mySM->getMut()
             << "   mub = " << mySM->getMub()
             << "   muc = " << mySM->getMuc() << endl << endl;
        
        cout << setw(30) << " " << setw(18) << "LO"
             << setw(18) << "FULLNLO"
             << setw(18) << "FULLNNLO" << endl;

        ////////////////////////////////////////////////////////////////////////        
        
        cout << "[Test for Lambda_QCD]" << endl;
        cout << setw(30) << "Lambda(6):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(200.0, (orders)order));
        cout << endl;
        cout << setw(30) << "Lambda(5):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda5((orders)order));
        cout << endl;
        cout << setw(30) << "Lambda(4):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(3.0, (orders)order));
        cout << endl;
        cout << setw(30) << "Lambda(3):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(1.0, (orders)order));
        cout << endl << endl;

        ////////////////////////////////////////////////////////////////////////
        
        cout << "[Test for QCD::Als(mu, order)]" << std::endl;
        cout << setw(30) << "alpha_s(M_t):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getMtpole(), (orders)order);
        cout << endl;
        cout << setw(30) << "alpha_s(M_Z):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getMz(), (orders)order);
        cout << endl;            
        cout << setw(30) << "alpha_s(m_b(m_b)):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getQuarks(mySM->BOTTOM).getMass(), (orders)order);
        cout << endl;            
        cout << setw(30) << "alpha_s(m_c(m_c)):";
        for (int order=LO; order<=FULLNNLO; order++) 
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getQuarks(mySM->CHARM).getMass(), (orders)order);
        cout << endl << endl;      

        ////////////////////////////////////////////////////////////////////////
        
        cout << "[Test for QCD::AlsWithLambda(mu, logLambda, nf, order)]" << endl;
        cout << setw(30) << "alpha_s(M_t) for Nf=6:";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getMtpole();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 6., (orders)order);
            }
        }
        cout << endl; 
        cout << setw(30) << "alpha_s(M_t) for Nf=5:";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getMtpole();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 5., (orders)order);
            }
        }
        cout << endl; 
        cout << setw(30) << "alpha_s(M_Z):";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getMz();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 5., (orders)order);
            }
        }
        cout << endl;         
        cout << setw(30) << "alpha_s(m_b(m_b)) for Nf=5:";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getQuarks(mySM->BOTTOM).getMass();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 5., (orders)order);
            }
        }
        cout << endl;         
        cout << setw(30) << "alpha_s(m_b(m_b)) for Nf=4:";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getQuarks(mySM->BOTTOM).getMass();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 4., (orders)order);
            }
        }
        cout << endl;         
        cout << setw(30) << "alpha_s(m_c(m_c)) for Nf=4:";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getQuarks(mySM->CHARM).getMass();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 4., (orders)order);
            }
        }
        cout << endl;        
        cout << setw(30) << "alpha_s(m_c(m_c)) for Nf=3:";
        for (int order=LO; order<=FULLNNLO; order++) {
            if (order!=NLO && order!=NNLO) {
                double mu = mySM->getQuarks(mySM->CHARM).getMass();
                double logLambda = mySM->logLambda(mu,(orders)order);
                cout << setw(18) << mySM->AlsWithLambda(mu, logLambda, 3., (orders)order);
            }
        }
        cout << endl << endl;        

        ////////////////////////////////////////////////////////////////////////
        
        cout << "1/alpha(M_Z) = " << 1.0/mySM->alphaMz() << endl;
        cout << "1/ale_OS(100.0, LO) = " << 1.0/mySM->ale_OS(100.0, LO) << endl;
        cout << "1/ale_OS(100.0, FULLNLO) = " << 1.0/mySM->ale_OS(100.0, FULLNLO) << endl;
        cout << "1/ale_OS(150.0, LO) = " << 1.0/mySM->ale_OS(150.0, LO) << endl;
        cout << "1/ale_OS(150.0, FULLNLO) = " << 1.0/mySM->ale_OS(150.0, FULLNLO) << endl;
        cout << "1/ale_OS(200.0, LO) = " << 1.0/mySM->ale_OS(200.0, LO) << endl;
        cout << "1/ale_OS(200.0, FULLNLO) = " << 1.0/mySM->ale_OS(200.0, FULLNLO) << endl;
        
        ////////////////////////////////////////////////////////////////////////
        
        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}


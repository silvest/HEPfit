/* 
 * File:   RunningMasses.cpp
 * Author: mishima
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
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
    
    // 26+16 parameters defined in QCD    
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

        cout << setw(20) << " " << setw(18) << "LO"
             << setw(18) << "FULLNLO"
             << setw(18) << "FULLNNLO" << endl;        
        
        ////////////////////////////////////////////////////////////////////////
        
        cout << setw(20) << "m_b(M_Z):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mb = mySM->getQuarks(mySM->BOTTOM).getMass();
                cout << setw(18) << mySM->Mrun(mySM->getMz(), mb, (orders)order);
            }
        }
        cout << endl;
        cout << setw(20) << "m_c(M_Z):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mc = mySM->getQuarks(mySM->CHARM).getMass();
                cout << setw(18) << mySM->Mrun(mySM->getMz(), mc, (orders)order);
            }
        }
        cout << endl << endl;

        cout << setw(20) << "m_b(m_b):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mb = mySM->getQuarks(mySM->BOTTOM).getMass();
                cout << setw(18) << mySM->Mrun(mb, mb, (orders)order);
            }
        }
        cout << endl;
        cout << setw(20) << "m_c(m_b):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mc = mySM->getQuarks(mySM->CHARM).getMass();
                double mb = mySM->getQuarks(mySM->BOTTOM).getMass();
                cout << setw(18) << mySM->Mrun(mb, mc, (orders)order);
            }
        }
        cout << endl << endl;;
        
        cout << setw(20) << "m_b(2.0):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mb = mySM->getQuarks(mySM->BOTTOM).getMass();
                cout << setw(18) << mySM->Mrun(2.0, mb, (orders)order);
            }
        }
        cout << endl;
        cout << setw(20) << "m_c(2.0):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mc = mySM->getQuarks(mySM->CHARM).getMass();
                cout << setw(18) << mySM->Mrun(2.0, mc, (orders)order);
            }
        }
        cout << endl << endl;
        
        cout << setw(20) << "m_b(m_c):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mb = mySM->getQuarks(mySM->BOTTOM).getMass();
                double mc = mySM->getQuarks(mySM->CHARM).getMass();
                cout << setw(18) << mySM->Mrun(mc, mb, (orders)order);
            }
        }
        cout << endl;
        cout << setw(20) << "m_c(m_c):";
        for (int order=LO; order<=FULLNNLO; order++) { 
            if (order!=NLO && order!=NNLO) { 
                double mc = mySM->getQuarks(mySM->CHARM).getMass();
                cout << setw(18) << mySM->Mrun(mc, mc, (orders)order);
            }
        }
        cout << endl << endl;;
       
        ////////////////////////////////////////////////////////////////////////        
        
        return EXIT_SUCCESS;
    } catch (const char* c) {
        cerr << c << endl;
        return EXIT_FAILURE;
    }        
}


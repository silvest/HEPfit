/* 
 * Copyright (C) 2012 HEPfit Collaboration
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

    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1875;
    Parameters["mtop"] = 173.2;
    //
    Parameters["mup"] = 0.0023;
    Parameters["mdown"] = 0.0048;
    Parameters["mstrange"] = 0.095;
    Parameters["mcharm"] = 1.275;
    Parameters["mbottom"] = 4.18;
    //
    Parameters["muc"] = 1.275;
    Parameters["mub"] = 4.18;
    Parameters["mut"] = 164.0;
    //Parameters["mut"] = 173.2;

    Parameters["mHl"] = 125.6;
    
    SM_i.Update(Parameters);    
}


int main(int argc, char** argv) {
    
    try {    
        StandardModel* mySM;
        mySM = new StandardModel();
        mySM->InitializeModel();
        setSMparameters(*mySM);

        cout.precision(8);
        //cout.setf(ios::floatfield);
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
        
        cout << "M_t_pole:               " << mySM->getMtpole() << endl;
        cout << "M_t(M_t) from M_t_pole: " << mySM->Mp2Mbar(mySM->getMtpole(), FULLNNLO) << endl;
        cout << "M_t_pole from M_t(M_t): " << mySM->Mbar2Mp(mySM->Mp2Mbar(mySM->getMtpole(), FULLNNLO)) 
             << endl << endl;

        cout << "M_t_pole from m_t=163.3: " << mySM->Mbar2Mp(163.3) 
             << endl << endl;

        
        cout << "m_b(m_b):               " << mySM->getQuarks(mySM->BOTTOM).getMass() << endl;
        cout << "M_b_pole from m_b(m_b): " << mySM->Mbar2Mp(mySM->getQuarks(mySM->BOTTOM).getMass(), FULLNNLO) << endl;
        cout << "m_b(m_b) from M_b_pole: " << mySM->Mp2Mbar(mySM->Mbar2Mp(mySM->getQuarks(mySM->BOTTOM).getMass(), FULLNNLO)) 
             << endl << endl;
        
        ////////////////////////////////////////////////////////////////////////        
        
        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }   
}


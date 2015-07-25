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

#define EPS 1.e-10


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

        cout << "Inputs: " << endl
             << "   M_Z = " << mySM->getMz()
             << "   alpha_s(M_Z) = " << mySM->getAlsMz() << endl
             << "   M_t = " << mySM->getMtpole()
             << "   m_b(m_b) = " << mySM->getQuarks(mySM->BOTTOM).getMass()
             << "   m_c(m_c) = " << mySM->getQuarks(mySM->CHARM).getMass()
             << "   m_s(2) = " << mySM->getQuarks(mySM->STRANGE).getMass() << endl
             << "   mut = " << mySM->getMut()
             << "   mub = " << mySM->getMub()
             << "   muc = " << mySM->getMuc() << endl << endl;

        cout << "Decoupling scales (inputs): " << endl
             << "  mut = " << mySM->getMut()
             << "  mub = " << mySM->getMub()
             << "  muc = " << mySM->getMuc()
             << endl;

        cout << "Pole masses: " << endl
             << "  M_t = " << mySM->getMtpole() << endl
             << "  M_b = " << mySM->Mbar2Mp(mySM->getQuarks(mySM->BOTTOM).getMass(), LO) << "[LO], "
             << mySM->Mbar2Mp(mySM->getQuarks(mySM->BOTTOM).getMass(), FULLNLO) << "[NLO], "
             << mySM->Mbar2Mp(mySM->getQuarks(mySM->BOTTOM).getMass(), FULLNNLO) << "[NNLO]"
             //<< "  M_c = " << mySM->Mbar2Mp(mySM->getQuarks(mySM->CHARM).getMass())
             << endl;

        cout << "MSbar invariant masses: " << endl
             << "  m_t(m_t) = "
             << mySM->Mp2Mbar(mySM->getMtpole(), LO) << "[LO], "
             << mySM->Mp2Mbar(mySM->getMtpole(), FULLNLO) << "[NLO], "
             << mySM->Mp2Mbar(mySM->getMtpole(), FULLNNLO) << "[NNLO], "
             << mySM->getQuarks(mySM->TOP).getMass() << "[getMass]" << endl
             << "  m_b(m_b) = " << mySM->getQuarks(mySM->BOTTOM).getMass() << endl
             << "  m_c(m_c) = " << mySM->getQuarks(mySM->CHARM).getMass()
             << endl << endl;

        ////////////////////////////////////////////////////////////////////////

        cout << "[Test for Lambda_QCD]" << endl;
        cout << setw(30) << " " << setw(18) << "LO"
             << setw(18) << "FULLNLO"
             << setw(18) << "FULLNNLO" << endl;
        cout << setw(30) << "Lambda(6):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(6., (orders)order));
        cout << endl;
        cout << setw(30) << "Lambda(5):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(5., (orders)order));
        cout << endl;
        cout << setw(30) << "Lambda(4):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(4., (orders)order));
        cout << endl;
        cout << setw(30) << "Lambda(3):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << exp(mySM->logLambda(3., (orders)order));
        cout << endl << endl;

        ////////////////////////////////////////////////////////////////////////

        cout << "[Test for QCD::AlsWithLambda(mu, order)]" << endl;
        cout << setw(30) << " " << setw(18) << "LO"
             << setw(18) << "NLO"
             << setw(18) << "NNLO"
             << setw(18) << "FULLNLO"
             << setw(18) << "FULLNNLO" << endl;
        cout << setw(30) << "alpha_s(M_t):";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMtpole();
            cout << setw(18) << mySM->AlsWithLambda(mu, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(mut) for Nf=6:";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMut();
            cout << setw(18) << mySM->AlsWithLambda(mu+EPS, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(mut) for Nf=5:";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMut();
            cout << setw(18) << mySM->AlsWithLambda(mu-EPS, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(M_Z):";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMz();
            cout << setw(18) << mySM->AlsWithLambda(mu, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(m_b(m_b)):";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getQuarks(mySM->BOTTOM).getMass();
            cout << setw(18) << mySM->AlsWithLambda(mu, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(mub) for Nf=5:";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMub();
            cout << setw(18) << mySM->AlsWithLambda(mu+EPS, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(mub) for Nf=4:";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMub();
            cout << setw(18) << mySM->AlsWithLambda(mu-EPS, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(m_c(m_c)):";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getQuarks(mySM->CHARM).getMass();
            cout << setw(18) << mySM->AlsWithLambda(mu, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(muc) for Nf=4:";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMuc();
            cout << setw(18) << mySM->AlsWithLambda(mu+EPS, (orders)order);
        }
        cout << endl;
        cout << setw(30) << "alpha_s(muc) for Nf=3:";
        for (int order=LO; order<=FULLNNLO; order++) {
            double mu = mySM->getMuc();
            cout << setw(18) << mySM->AlsWithLambda(mu-EPS, (orders)order);
        }
        cout << endl << endl;

        ////////////////////////////////////////////////////////////////////////

        cout << "[Test for QCD::Als(mu, order)]" << std::endl;
        cout << setw(30) << " " << setw(18) << "LO"
             << setw(18) << "FULLNLO"
             << setw(18) << "FULLNNLO" << endl;
        cout << setw(30) << "alpha_s(M_t):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getMtpole(), (orders)order);
        cout << endl;
        cout << setw(30) << "alpha_s(mut):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getMut(), (orders)order);
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
        cout << setw(30) << "alpha_s(mub):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getMub(), (orders)order);
        cout << endl;
        cout << setw(30) << "alpha_s(m_c(m_c)):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getQuarks(mySM->CHARM).getMass(), (orders)order);
        cout << endl;
        cout << setw(30) << "alpha_s(muc):";
        for (int order=LO; order<=FULLNNLO; order++)
            if (order!=NLO && order!=NNLO)
                cout << setw(18) << mySM->Als(mySM->getMuc(), (orders)order);
        cout << endl << endl;

        ////////////////////////////////////////////////////////////////////////

        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}


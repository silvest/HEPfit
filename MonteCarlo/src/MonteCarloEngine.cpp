/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MonteCarloEngine.h"
#include <BAT/BCParameter.h>
#include <BAT/BCMath.h>
#ifdef _MPI
#include <mpi.h>
#endif
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <fstream>
#include <stdexcept>
#include <iomanip>

MonteCarloEngine::MonteCarloEngine(
        const std::vector<ModelParameter>& ModPars_i,
        boost::ptr_vector<Observable>& Obs_i,
        std::vector<Observable2D>& Obs2D_i,
        std::vector<CorrelatedGaussianObservables>& CGO_i,
        std::vector<CorrelatedGaussianParameters>& CGP_i)
: BCModel(""), ModPars(ModPars_i), CGP(CGP_i), Obs_ALL(Obs_i), Obs2D_ALL(Obs2D_i),
  CGO(CGO_i), NumOfUsedEvents(0), NumOfDiscardedEvents(0) {
    obval = NULL;
    obweight = NULL;
    Mod = NULL;
    TH1::StatOverflows(kTRUE);
#ifdef _MPI
    rank = MPI::COMM_WORLD.Get_rank();
#else
    rank = 0;
#endif
};

void MonteCarloEngine::Initialize(Model* Mod_i) {
    Mod = Mod_i;
    int k = 0, kweight = 0;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            k++;
            if (it->getDistr().compare("noweight") != 0)
                kweight++;
        }
        std::string HistName = it->getName();
        if (Histo1D.find(HistName) == Histo1D.end()) {
            TH1D * histo = new TH1D(HistName.c_str(), it->getLabel().c_str(),
                    NBINS1D, it->getMin(), it->getMax());
            histo->GetXaxis()->SetTitle(it->getLabel().c_str());
            BCH1D * bchisto = new BCH1D(histo);
            Histo1D[HistName] = bchisto;
            thMin[it->getName()] = std::numeric_limits<double>::max();
            thMax[it->getName()] = -std::numeric_limits<double>::max();
        }
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            if (!it->isTMCMC())
                throw std::runtime_error("ERROR: cannot handle noMCMC for Observable2D file yet!");
        } else if (it->getDistr().compare("weight") == 0)
            throw std::runtime_error("ERROR: do not use Observable2D for analytic 2D weights!");
        std::string HistName = it->getName();
        if (Histo2D.find(HistName) == Histo2D.end()) {
            TH2D * histo2 = new TH2D(HistName.c_str(),
                    (it->getLabel() + " vs " + it->getLabel2()).c_str(),
                    NBINS2D, it->getMin(), it->getMax(),
                    NBINS2D, it->getMin2(), it->getMax2());
            histo2->GetXaxis()->SetTitle(it->getLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel2().c_str());
            BCH2D * bchisto2 = new BCH2D(histo2);
            Histo2D[HistName] = bchisto2;
        }
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 != CGO.end(); ++it1) {
        std::vector<Observable> pino(it1->getObs());
        for (std::vector<Observable>::iterator it = pino.begin();
                it != pino.end(); ++it) {
            //for (int i = 0; i < it1->getObs().size(); i++) {
            //Observable * it = &(it1->getObs().at(i));
            if ((it->getDistr()).compare("file") == 0)
                throw std::runtime_error("Cannot use file in CorrelatedGaussianObservables!");
            if (!(it->isTMCMC())) {
                k++;
                if (it->getDistr().compare("noweight") != 0)
                    throw std::runtime_error("Cannot use weight in CorrelatedGaussianObservables!");
            }
            std::string HistName = it->getName();
            if (Histo1D.find(HistName) == Histo1D.end()) {
                TH1D * histo = new TH1D(HistName.c_str(), it->getLabel().c_str(),
                        NBINS1D, it->getMin(), it->getMax());
                histo->GetXaxis()->SetTitle(it->getLabel().c_str());
                BCH1D * bchisto = new BCH1D(histo);
                Histo1D[HistName] = bchisto;
                thMin[HistName] = std::numeric_limits<double>::max();
                thMax[HistName] = -std::numeric_limits<double>::max();
            }
        }
    }
    kmax = k;
    kwmax = kweight;

    DefineParameters();

};

void MonteCarloEngine::setNChains(unsigned int i) {
    MCMCSetNChains(i);
    obval = new double[fMCMCNChains * kmax];
    obweight = new double[fMCMCNChains * kwmax];
}

// ---------------------------------------------------------

MonteCarloEngine::~MonteCarloEngine()
// default destructor
{
    delete [] obval;
    delete [] obweight;
    /* The following code has been commented out pending further review.
       It is causing crashes at the termination of the code if the histograms
       are accessed from the main program.*/
    //    for (std::map<std::string, BCH1D *>::iterator it = Histo1D.begin();
    //            it != Histo1D.end(); it++)
    //        delete it->second;
    //    for (std::map<std::string, BCH2D *>::iterator it = Histo2D.begin();
    //            it != Histo2D.end(); it++)
    //        delete it->second;
};

// ---------------------------------------------------------

void MonteCarloEngine::DefineParameters() {
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    if (rank == 0) std::cout << "\nParameters varied in Monte Carlo:" << std::endl;
    int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin();
            it < ModPars.end(); it++) {
        if (it->geterrf() == 0. && it->geterrg() == 0.)
            continue;
        AddParameter(it->getname().c_str(), it->getmin(), it->getmax());
        if (rank == 0) std::cout << "  " << it->getname() << " " << k << std::endl;
        if (it->IsCorrelated()) {
            for (int i = 0; i < CGP.size(); i++) {
                if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                    std::string index = it->getname().substr(CGP[i].getName().size());
                    long int lindex = strtol(index.c_str(),NULL,10);
                    if (lindex > 0)
                        DPars[CGP[i].getPar(lindex - 1).getname()] = 0.;
                    else {
                        std::stringstream out;
                        out << it->getname();
                        throw std::runtime_error("MonteCarloEngine::DefineParameters(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                    }
                }
            }
        } else
            DPars[it->getname()] = 0.;
        if (it->geterrf() == 0.) SetPriorGauss(k, it->getave(), it->geterrg());
        else if (it->geterrg() == 0.) SetPriorConstant(k);
        else {
            TF1 * combined = new TF1(it->getname().c_str(),
                    "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                    it->getmin(), it->getmax());
            combined->SetParameter(0, it->getave());
            combined->SetParameter(1, it->geterrg());
            combined->SetParameter(2, it->geterrf());
            SetPrior(k, combined);
            delete combined;
        }
        k++;
    }
}

void MonteCarloEngine::setDParsFromParameters(const std::vector<double>& parameters, 
        std::map<std::string,double>& DPars_i) 
{
    std::map<std::string, std::vector<double> > cgpmap;

    unsigned int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it < ModPars.end(); it++){
        if(it->IsFixed())
            continue;
        if(it->getname().compare(GetParameter(k)->GetName()) != 0)
            {
                        std::stringstream out;
                        out << it->getname();
                        throw std::runtime_error("MonteCarloEngine::setDParsFromParameters(): " + out.str() + "is sitting at the wrong position in the BAT parameters vector");
                    }
        if (it->IsCorrelated()) {
            std::string index = it->getname().substr(it->getCgp_name().size());
            long int lindex = strtol(index.c_str(),NULL,10);
            if (lindex - 1 == cgpmap[it->getCgp_name()].size())
                cgpmap[it->getCgp_name()].push_back(parameters[k]);
            else {
                std::stringstream out;
                out << it->getname() << " " << lindex;
                throw std::runtime_error("MonteCarloEngine::setDParsFromParameters(): " + out.str() + "seems to be a CorrelatedGaussianParameters object but the corresponding parameters are missing or not in the right order");
            }

        } else
            DPars_i[it->getname()] = parameters[k];
        k++;
    }

    for (int j = 0; j < CGP.size(); j++) {
        std::vector<double> current = cgpmap.at(CGP[j].getName());
        if (current.size() != CGP[j].getPars().size()) {
            std::stringstream out;
            out << CGP[j].getName();
            throw std::runtime_error("MonteCarloEngine::setDParsFromParameters(): " + out.str() + " appears to be represented in cgpmap with a wrong size");
        }
        
        std::vector<double> porig = CGP[j].getOrigParsValue(current);

        for(int l = 0; l < porig.size(); l++) {
            DPars_i[CGP[j].getPar(l).getname()] = porig[l];
        }
        

    }

}

// ---------------------------------------------------------

double MonteCarloEngine::LogLikelihood(const std::vector<double>& parameters) {
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.

    double logprob = 0.;
    
    setDParsFromParameters(parameters, DPars);

    // if update false set probability equal zero
    if (!Mod->Update(DPars)) {
#ifdef _MCDEBUG
        std::cout << "event discarded" << std::endl;

        /* Debug */
        //for (int k = 0; k < parameters.size(); k++)
        //    std::cout << "  " << GetParameter(k)->GetName() << " = "
        //              << DPars[GetParameter(k)->GetName()] << std::endl;
#endif
        NumOfDiscardedEvents++;
        return (log(0.));
    }
#ifdef _MCDEBUG
    //std::cout << "event used in MC" << std::endl;
#endif

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it != Obs_ALL.end(); it++) {
        if (it->isTMCMC())
            logprob += it->computeWeight();
    }

    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it != Obs2D_ALL.end(); it++) {
        if (it->isTMCMC())
            logprob += it->computeWeight();
    }

    for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++) {
        logprob += it->computeWeight();
    }
    if (isnan(logprob)) {
        NumOfDiscardedEvents++;
        std::cout << "Event discarded since logprob evaluated to NAN.\n";
        return (log(0.));
    }
    NumOfUsedEvents++;
    return logprob;
}

void MonteCarloEngine::MCMCIterationInterface() {
#ifdef _MPI
    unsigned mychain = 0;
    int iproc = 0;
    unsigned npars = GetNParameters();
    int buffsize = npars + 1;
    int index_chain[procnum];
    double *recvbuff = new double[buffsize];
    std::vector<std::vector<double> > fMCMCxvect;
    std::vector<double> pars;
    double **buff;


    buff = new double*[procnum];
    int obsbuffsize = Obs_ALL.size() + 2 * Obs2D_ALL.size();
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++)
        obsbuffsize += it1->getObs().size();
    buff[0] = new double[procnum * obsbuffsize];
    for (int i = 1; i < procnum; i++) {
        buff[i] = buff[i - 1] + obsbuffsize;
        index_chain[i] = -1;
    }

    double ** sendbuff = new double *[procnum];
    sendbuff[0] = new double[procnum * buffsize];
    for (int il = 1; il < procnum; il++)
        sendbuff[il] = sendbuff[il - 1] + buffsize;

    while (mychain < fMCMCNChains) {
        pars.clear();
        for (unsigned int k = 0; k < npars; k++)
            pars.push_back(fMCMCx.at(mychain * npars + k));

        fMCMCxvect.push_back(pars);
        index_chain[iproc] = mychain;
        iproc++;
        mychain++;
        if (iproc < procnum && mychain < fMCMCNChains)
            continue;

        for (int il = 0; il < iproc; il++) {
            //The first entry of the array specifies the task to be executed.

            sendbuff[il][0] = 2.; // 2 = observables calculation
            for (int im = 1; im < buffsize; im++)
                sendbuff[il][im] = fMCMCxvect[il][im - 1];
        }
        for (int il = iproc; il < procnum; il++) {
            sendbuff[il][0] = 3.; // 3 = nothing to execute, but return a buffer of observables
            index_chain[il] = -1;
        }
        //       double inittime = MPI::Wtime();
        MPI::COMM_WORLD.Scatter(sendbuff[0], buffsize, MPI::DOUBLE,
                recvbuff, buffsize, MPI::DOUBLE,
                0);

        if (recvbuff[0] == 2.) { // compute observables
            double sbuff[obsbuffsize];
            std::map<std::string, double> DPars;
            pars.assign(recvbuff + 1, recvbuff + buffsize);
            setDParsFromParameters(pars,DPars);
            Mod->Update(DPars);

            int k = 0;
            for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
                sbuff[k++] = it->computeTheoryValue();
            }
            for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end(); it++) {
                sbuff[k++] = it->computeTheoryValue();
                sbuff[k++] = it->computeTheoryValue2();
            }

            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                std::vector<Observable> pino(it1->getObs());
                for (std::vector<Observable>::iterator it = pino.begin(); it != pino.end(); ++it)
                    sbuff[k++] = it->computeTheoryValue();
            }
            MPI::COMM_WORLD.Gather(sbuff, obsbuffsize, MPI::DOUBLE, buff[0], obsbuffsize, MPI::DOUBLE, 0);
        } else if (recvbuff[0] == 3.) { // do not compute observables, but gather the buffer
            double sbuff[obsbuffsize];
            MPI::COMM_WORLD.Gather(sbuff, obsbuffsize, MPI::DOUBLE, buff[0], obsbuffsize, MPI::DOUBLE, 0);
        }

        for (int il = 0; il < procnum; il++) {
            if (index_chain[il] >= 0) {
                int k = 0;
                // fill the histograms for observables
                int ko = 0, kweight = 0;
                for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                        it < Obs_ALL.end(); it++) {
                    double th = buff[il][k++];
                    /* set the min and max of theory values */
                    if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                    if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                    Histo1D[it->getName()]->GetHistogram()->Fill(th);
                    if (!it->isTMCMC()) {
                        obval[index_chain[il] * kmax + ko++] = th;
                        if (it->getDistr().compare("noweight") != 0) {
                            obweight[index_chain[il] * kwmax + kweight++] = it->computeWeight(th);
                        }
                    }
                }

                // fill the 2D histograms for observables
                for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
                        it < Obs2D_ALL.end(); it++) {
                    double th1 = buff[il][k++];
                    double th2 = buff[il][k++];
                    Histo2D[it->getName()]->GetHistogram()->Fill(th1, th2);
                }

                // fill the histograms for correlated observables
                for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                        it1 < CGO.end(); it1++) {
                    std::vector<Observable> pino(it1->getObs());
                    for (std::vector<Observable>::iterator it = pino.begin();
                            it != pino.end(); ++it) {
                        double th = buff[il][k++];
                        /* set the min and max of theory values */
                        if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                        if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                        Histo1D[it->getName()]->GetHistogram()->Fill(th);
                    }
                }

            }
        }
        iproc = 0;
        fMCMCxvect.clear();
    }
    delete sendbuff[0];
    delete [] sendbuff;
    delete [] recvbuff;
    delete buff[0];
    delete [] buff;
#else
    for (unsigned int i = 0; i < fMCMCNChains; ++i) {
        std::vector<double>::const_iterator first = fMCMCx.begin() + i * GetNParameters();
        std::vector<double>::const_iterator last = first + GetNParameters();
        std::vector<double> currvec(first, last);
        setDParsFromParameters(currvec,DPars);

        Mod->Update(DPars);

        // fill the histograms for observables
        int k = 0, kweight = 0;
        for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) {
            double th = it->computeTheoryValue();
            /* set the min and max of theory values */
            if (th < thMin[it->getName()]) thMin[it->getName()] = th;
            if (th > thMax[it->getName()]) thMax[it->getName()] = th;
            Histo1D[it->getName()]->GetHistogram()->Fill(th);
            if (!it->isTMCMC()) {
                obval[i * kmax + k] = th;
                k++;
                if (it->getDistr().compare("noweight") != 0) {
                    obweight[i * kwmax + kweight] = it->computeWeight();
                    kweight++;
                }
            }
        }

        // fill the 2D histograms for observables
        for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
                it < Obs2D_ALL.end(); it++) {
            double th1 = it->computeTheoryValue();
            double th2 = it->computeTheoryValue2();
            Histo2D[it->getName()]->GetHistogram()->Fill(th1, th2);
        }

        // fill the histograms for correlated observables
        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                it1 < CGO.end(); it1++) {
            std::vector<Observable> pino(it1->getObs());
            for (std::vector<Observable>::iterator it = pino.begin();
                    it != pino.end(); ++it) {
                double th = it->computeTheoryValue();
                /* set the min and max of theory values */
                if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                Histo1D[it->getName()]->GetHistogram()->Fill(th);
            }
        }
    }
#endif
}

void MonteCarloEngine::CheckHistogram(const TH1D& hist, const std::string name) {
    // output the portions of underflow and overflow bins
    double UnderFlowContent = hist.GetBinContent(0);
    double OverFlowContent = hist.GetBinContent(NBINS1D + 1);
    double Integral = hist.Integral();
    double TotalContent = 0.0;
    for (int n = 0; n <= NBINS1D + 1; n++)
        TotalContent += hist.GetBinContent(n);
    HistoLog << name << ": "
            << Integral / TotalContent * 100. << "% within the range, "
            << UnderFlowContent / TotalContent * 100. << "% underflow, "
            << OverFlowContent / TotalContent * 100. << "% overflow"
            << std::endl;
}

void MonteCarloEngine::CheckHistogram(const TH2D& hist, const std::string name) {
    double Integral = hist.Integral();
    double TotalContent = 0.0;
    for (int m = 0; m <= NBINS2D + 1; m++)
        for (int n = 0; n <= NBINS2D + 1; n++)
            TotalContent += hist.GetBinContent(m, n);
    HistoLog << name << ": "
            << Integral / TotalContent * 100. << "% within the ranges"
            << std::endl;
}

void MonteCarloEngine::PrintHistogram(BCModelOutput & out, Observable& it,
        const std::string OutputDir) {
    std::string HistName = it.getName();
    double min = thMin[it.getName()];
    double max = thMax[it.getName()];

    if (Histo1D[HistName]->GetHistogram()->Integral() > 0.0) {
        std::string fname = OutputDir + "/" + HistName + ".pdf";
        //        BCH1D* pippo =  Histo1D[HistName];
        //        double x = pippo->GetMean();
        //        pippo->Print("Dmd1.pdf");
        Histo1D[HistName]->SetGlobalMode(it.computeTheoryValue());
        Histo1D[HistName]->Print(fname.c_str());
        std::cout << fname << " has been created." << std::endl;
        out.Write(Histo1D[HistName]->GetHistogram());
        CheckHistogram(*Histo1D[HistName]->GetHistogram(), it.getName());
    } else
        HistoLog << "WARNING: The histogram of "
            << it.getName() << " is empty!" << std::endl;

    HistoLog.precision(10);
    HistoLog << "  [min, max]=[" << min << ", " << max << "]" << std::endl;
    HistoLog.precision(6);
}

void MonteCarloEngine::PrintHistogram(BCModelOutput & out, const std::string OutputDir) {
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) {
        if (rank == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables histogram cannot be generated.\n");
        else sleep(2);
    }
    setDParsFromParameters(mode,DPars);

    Mod->Update(DPars);

    // print the histograms to pdf files
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end();
            it++) {
        PrintHistogram(out, *it, OutputDir);
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 < CGO.end(); it1++) {
        std::vector<Observable> pino(it1->getObs());
        for (std::vector<Observable>::iterator it = pino.begin();
                it != pino.end(); ++it)
            PrintHistogram(out, *it, OutputDir);
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        std::string HistName = it->getName();
        if (Histo2D[HistName]->GetHistogram()->Integral() > 0.0) {
            std::string fname = OutputDir + "/" + HistName + ".pdf";
            double th[2];
            th[0] = it->computeTheoryValue();
            th[1] = it->computeTheoryValue2();
            Histo2D[HistName]->SetGlobalMode(th);
            Histo2D[HistName]->Print(fname.c_str());
            std::cout << fname << " has been created." << std::endl;
            out.Write(Histo2D[HistName]->GetHistogram());
            CheckHistogram(*Histo2D[HistName]->GetHistogram(), HistName);
        } else
            HistoLog << "WARNING: The histogram of "
                << HistName << " is empty!" << std::endl;
    }
}

void MonteCarloEngine::AddChains() {
    fMCMCFlagWritePreRunToFile = false;
    int k = 0, kweight = 0;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            for (unsigned int i = 0; i < fMCMCNChains; ++i)
                fMCMCTrees[i]->Branch(it->getName().c_str(), &obval[i * kmax + k],
                    (it->getName() + "/D").c_str());
            k++;
            if (it->getDistr().compare("noweight") != 0) {
                for (unsigned int i = 0; i < fMCMCNChains; ++i)
                    fMCMCTrees[i]->Branch((it->getName() + "_weight").c_str(),
                        &obweight[i * kwmax + kweight],
                        (it->getName() + "_weight/D").c_str());
                kweight++;
            }
        }
    }
}

void MonteCarloEngine::PrintCorrelationMatrix(const std::string filename) {
    std::ofstream out;
    out.open(filename.c_str(), std::ios::out);

    int npar = GetNParameters();

    for (int i = 0; i < npar; ++i)
        out << " & " << GetParameter(i)->GetName();
    out << " \\\\" << std::endl;

    for (int i = 0; i < npar; ++i) {
        out << GetParameter(i)->GetName() << " & $";
        for (int j = 0; j < npar; ++j) {
            if (i != j) {
                BCH2D* bch2d_temp = GetMarginalized(GetParameter(i), GetParameter(j));
                if (bch2d_temp != NULL)
                    out << bch2d_temp->GetHistogram()->GetCorrelationFactor();
                else
                    out << 0.;
            } else
                out << 1.;
            if (j == npar - 1) out << "$ \\\\" << std::endl;
            else out << "$ & $";
        }
    }

    out.close();
}

std::string MonteCarloEngine::computeStatistics() {
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) {
        if(rank == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables statistics cannot be generated.\n");
        else sleep (2);
    }
    std::ostringstream StatsLog;
    int i = 0;
    StatsLog << "Statistics file for Observables, Binned Observables and Corellated Gaussian Observables.\n" << std::endl;
    if (Obs_ALL.size() > 0) StatsLog << "Observables:\n" << std::endl;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {

        if (it->getObsType().compare("BinnedObservable") == 0) {
            StatsLog << "  (" << ++i << ") Binned Observable \"";
            StatsLog << it->getName() << "[" << it->getTho()->getBinMin() << ", " << it->getTho()->getBinMax() << "]" << "\":";
        } else if (it->getObsType().compare("FunctionObservable") == 0) {
            StatsLog << "  (" << ++i << ") Function Observable \"";
            StatsLog << it->getName() << "[" << it->getTho()->getBinMin() << "]" << "\":";
        } else if (it->getObsType().compare("HiggsObservable") == 0) {
            StatsLog << "  (" << ++i << ") Higgs Observable \"";
            StatsLog << it->getName() << "\":";
        } else {
            StatsLog << "  (" << ++i << ") " << it->getObsType() << " \"";
            StatsLog << it->getName() << "\":";
        }

        StatsLog << std::endl;

        BCH1D * bch1d = Histo1D[it->getName()];
        if (bch1d->GetHistogram()->Integral() > 0.0) {
            StatsLog << "      Mean +- sqrt(V):                " << std::setprecision(4)
                    << bch1d->GetMean() << " +- " << std::setprecision(4)
                    << bch1d->GetRMS() << std::endl

                    << "      (Marginalized) mode:            " << bch1d->GetMode() << std::endl;

            std::vector<double> v1;
            std::vector<double> v2;
            std::vector<double> v3;
            v1 = bch1d->GetSmallestIntervals(0.6827);
            v2 = bch1d->GetSmallestIntervals(0.9545);
            v3 = bch1d->GetSmallestIntervals(0.9973);
            StatsLog << "      Smallest interval(s) containing at least 68.27% and local mode(s):"
                    << std::endl;
            for (unsigned j = 0; j < v1.size(); j += 5)
                StatsLog << "       (" << v1[j] << ", " << v1[j + 1]
                    << ") (local mode at " << v1[j + 3] << " with rel. height "
                    << v1[j + 2] << "; rel. area " << v1[j + 4] << ")"
                    << std::endl;
            StatsLog << std::endl;

            StatsLog << "      Smallest interval(s) containing at least 95.45% and local mode(s):"
                    << std::endl;
            for (unsigned j = 0; j < v2.size(); j += 5)
                StatsLog << "       (" << v2[j] << ", " << v2[j + 1]
                    << ") (local mode at " << v2[j + 3] << " with rel. height "
                    << v2[j + 2] << "; rel. area " << v2[j + 4] << ")"
                    << std::endl;
            StatsLog << std::endl;

            StatsLog << "      Smallest interval(s) containing at least 99.73% and local mode(s):"
                    << std::endl;
            for (unsigned j = 0; j < v3.size(); j += 5)
                StatsLog << "       (" << v3[j] << ", " << v3[j + 1]
                    << ") (local mode at " << v3[j + 3] << " with rel. height "
                    << v3[j + 2] << "; rel. area " << v3[j + 4] << ")"
                    << std::endl;
            StatsLog << std::endl;
        } else {
            StatsLog << "\nWARNING: The histogram of " << it->getName() << " is empty! Statistics cannot be generated\n" << std::endl;
        }
    }

    if (CGO.size() > 0) StatsLog << "\nCorellated Gaussian Observables:\n" << std::endl;
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 < CGO.end(); it1++) {
        StatsLog << "\n" << it1->getName() << ":\n" << std::endl;
        i = 0;
        std::vector<Observable> CGObs(it1->getObs());
        for (std::vector<Observable>::iterator it2 = CGObs.begin(); it2 < CGObs.end(); it2++) {

            if (it2->getObsType().compare("BinnedObservable") == 0) {
                StatsLog << "  (" << ++i << ") Binned Observable \"";
                StatsLog << it2->getName() << "[" << it2->getTho()->getBinMin() << ", " << it2->getTho()->getBinMax() << "]" << "\":";
            } else if (it2->getObsType().compare("FunctionObservable") == 0) {
                StatsLog << "  (" << ++i << ") Function Observable \"";
                StatsLog << it2->getName() << "[" << it2->getTho()->getBinMin() << "]" << "\":";
            } else if (it2->getObsType().compare("HiggsObservable") == 0) {
                StatsLog << "  (" << ++i << ") Higgs Observable \"";
                StatsLog << it2->getName() << "\":";
            } else {
                StatsLog << "  (" << ++i << ") " << it2->getObsType() << " \"";
                StatsLog << it2->getName() << "\":";
            }

            StatsLog << std::endl;
            BCH1D * bch1d = Histo1D[it2->getName()];
            if (bch1d->GetHistogram()->Integral() > 0.0) {
                StatsLog << "      Mean +- sqrt(V):                " << std::setprecision(4)
                        << bch1d->GetMean() << " +- " << std::setprecision(4)
                        << bch1d->GetRMS() << std::endl

                        << "      (Marginalized) mode:            " << bch1d->GetMode() << std::endl;

                std::vector<double> v1;
                std::vector<double> v2;
                std::vector<double> v3;
                v1 = bch1d->GetSmallestIntervals(0.6827);
                v2 = bch1d->GetSmallestIntervals(0.9545);
                v3 = bch1d->GetSmallestIntervals(0.9973);
                StatsLog << "      Smallest interval(s) containing at least 68.27% and local mode(s):"
                        << std::endl;
                for (unsigned j = 0; j < v1.size(); j += 5)
                    StatsLog << "       (" << v1[j] << ", " << v1[j + 1]
                        << ") (local mode at " << v1[j + 3] << " with rel. height "
                        << v1[j + 2] << "; rel. area " << v1[j + 4] << ")"
                        << std::endl;
                StatsLog << std::endl;

                StatsLog << "      Smallest interval(s) containing at least 95.45% and local mode(s):"
                        << std::endl;
                for (unsigned j = 0; j < v2.size(); j += 5)
                    StatsLog << "       (" << v2[j] << ", " << v2[j + 1]
                        << ") (local mode at " << v2[j + 3] << " with rel. height "
                        << v2[j + 2] << "; rel. area " << v2[j + 4] << ")"
                        << std::endl;
                StatsLog << std::endl;

                StatsLog << "      Smallest interval(s) containing at least 99.73% and local mode(s):"
                        << std::endl;
                for (unsigned j = 0; j < v3.size(); j += 5)
                    StatsLog << "       (" << v3[j] << ", " << v3[j + 1]
                        << ") (local mode at " << v3[j + 3] << " with rel. height "
                        << v3[j + 2] << "; rel. area " << v3[j + 4] << ")"
                        << std::endl;
                StatsLog << std::endl;
            } else {
                StatsLog << "\nWARNING: The histogram of " << it2->getName() << " is empty! Statistics cannot be generated\n" << std::endl;
            }
        }
    }

    return StatsLog.str().c_str();
}

std::string MonteCarloEngine::writePreRunData() 
{
    std::vector<double> mode(GetBestFitParameters());
    std::vector<double> scales(MCMCGetTrialFunctionScaleFactor(0));
    std::ostringstream StatsLog;
    for (int i = 0; i < mode.size(); i++)
        StatsLog << GetParameter(i)->GetName() << " " << mode.at(i) << " " << scales.at(i) << std::endl;
    return StatsLog.str().c_str();
}

double MonteCarloEngine::computeNormalization() {

    unsigned int Npars = GetNParameters();
    std::vector<double> mode(GetBestFitParameters());
    gslpp::matrix<double> Hessian(Npars, Npars, 0.);

    for (unsigned int i = 0; i < Npars; i++)
        for (unsigned int j = 0; j < Npars; j++) {
            // calculate Hessian matrix element
            Hessian.assign(i, j, -SecondDerivative(GetParameter(i), GetParameter(j), GetBestFitParameters()));
        }

    double det_Hessian = Hessian.determinant();

    return exp(Npars / 2. * log(2. * M_PI) + 0.5 * log(1. / det_Hessian) + LogLikelihood(mode) + LogAPrioriProbability(mode));
}

double MonteCarloEngine::SecondDerivative(const BCParameter * par1, const BCParameter * par2, std::vector<double> point) {

    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MCMCENgine::SecondDerivative : Invalid number of entries in the vector.");
    }

    // define steps
    const double dy1 = par2->GetRangeWidth() / NSTEPS;
    const double dy2 = dy1 * 2.;
    const double dy3 = dy1 * 3.;

    // define points at which to evaluate
    std::vector<double> y1p = point;
    std::vector<double> y1m = point;
    std::vector<double> y2p = point;
    std::vector<double> y2m = point;
    std::vector<double> y3p = point;
    std::vector<double> y3m = point;

    unsigned idy = GetParameters().Index(par2->GetName());

    y1p[idy] += dy1;
    y1m[idy] -= dy1;
    y2p[idy] += dy2;
    y2m[idy] -= dy2;
    y3p[idy] += dy3;
    y3m[idy] -= dy3;

    const double m1 = (FirstDerivative(par1, y1p) - FirstDerivative(par1, y1m)) / 2. / dy1;
    const double m2 = (FirstDerivative(par1, y2p) - FirstDerivative(par1, y2m)) / 4. / dy1;
    const double m3 = (FirstDerivative(par1, y3p) - FirstDerivative(par1, y3m)) / 6. / dy1;

    return 3. / 2. * m1 - 3. / 5. * m2 + 1. / 10. * m3;
}

double MonteCarloEngine::FirstDerivative(const BCParameter * par, std::vector<double> point) {

    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MCMCENgine::FirstDerivative : Invalid number of entries in the vector.");
    }

    // define steps
    const double dx1 = par->GetRangeWidth() / NSTEPS;
    const double dx2 = dx1 * 2.;
    const double dx3 = dx1 * 3.;

    // define points at which to evaluate
    std::vector<double> x1p = point;
    std::vector<double> x1m = point;
    std::vector<double> x2p = point;
    std::vector<double> x2m = point;
    std::vector<double> x3p = point;
    std::vector<double> x3m = point;

    unsigned idx = GetParameters().Index(par->GetName());

    x1p[idx] += dx1;
    x1m[idx] -= dx1;
    x2p[idx] += dx2;
    x2m[idx] -= dx2;
    x3p[idx] += dx3;
    x3m[idx] -= dx3;

    const double m1 = (Function_h(x1p) - Function_h(x1m)) / 2. / dx1;
    const double m2 = (Function_h(x2p) - Function_h(x2m)) / 4. / dx1;
    const double m3 = (Function_h(x3p) - Function_h(x3m)) / 6. / dx1;

    return 3. / 2. * m1 - 3. / 5. * m2 + 1. / 10. * m3;
}

double MonteCarloEngine::Function_h(std::vector<double> point) {
    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MCMCENgine::Function_h : Invalid number of entries in the vector.");
    }
    return LogLikelihood(point) + LogAPrioriProbability(point);
}
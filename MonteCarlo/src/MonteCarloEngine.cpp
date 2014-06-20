/* 
 * Copyright (C) 2012 SusyFit Collaboration
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

MonteCarloEngine::MonteCarloEngine(
        const std::vector<ModelParameter>& ModPars_i,
        boost::ptr_vector<Observable>& Obs_i,
        std::vector<Observable2D>& Obs2D_i,
        std::vector<CorrelatedGaussianObservables>& CGO_i)
: BCModel(""), ModPars(ModPars_i), Obs_ALL(Obs_i), Obs2D_ALL(Obs2D_i),
CGO(CGO_i), NumOfUsedEvents(0), NumOfDiscardedEvents(0)
{
    obval = NULL;
    obweight = NULL;
    Mod = NULL;
#ifdef _MPI
    rank = MPI::COMM_WORLD.Get_rank();
#else
    rank = 0;
#endif
};

void MonteCarloEngine::Initialize(Model* Mod_i)
{
    Mod = Mod_i;
    int k = 0, kweight = 0;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            k++;
            if (it->getDistr().compare("noweight") != 0)
                kweight++;
        }
        if (Histo1D.find(it->getThname()) == Histo1D.end()) {
            TH1D * histo = new TH1D(it->getThname().c_str(), it->getLabel().c_str(), NBINS1D,
                    it->getMin(), it->getMax());
            histo->GetXaxis()->SetTitle(it->getLabel().c_str());
            BCH1D * bchisto = new BCH1D(histo);
            Histo1D[it->getThname()] = bchisto;
            thMin[it->getThname()] = std::numeric_limits<double>::max();
            thMax[it->getThname()] = -std::numeric_limits<double>::max();
        }
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            if (!it->isTMCMC())
                throw std::runtime_error("ERROR: cannot handle noMCMC for Observable2D file yet!");
        } else if (it->getDistr().compare("weight") == 0)
            throw std::runtime_error("ERROR: do not use Observable2D for analytic 2D weights!");
        if (Histo2D.find(it->getThname() + "_vs_" + it->getThname2()) == Histo2D.end()) {
            TH2D * histo2 = new TH2D((it->getThname() + "_vs_" + it->getThname2()).c_str(),
                    (it->getLabel() + " vs " + it->getLabel2()).c_str(),
                    NBINS2D, it->getMin(), it->getMax(),
                    NBINS2D, it->getMin2(), it->getMax2());
            histo2->GetXaxis()->SetTitle(it->getLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel2().c_str());
            BCH2D * bchisto2 = new BCH2D(histo2);
            Histo2D[it->getThname() + "_vs_" + it->getThname2()] = bchisto2;
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
            if (Histo1D.find(it->getThname()) == Histo1D.end()) {
                TH1D * histo = new TH1D(it->getThname().c_str(), it->getLabel().c_str(),
                        NBINS1D, it->getMin(), it->getMax());
                histo->GetXaxis()->SetTitle(it->getLabel().c_str());
                BCH1D * bchisto = new BCH1D(histo);
                Histo1D[it->getThname()] = bchisto;
                thMin[it->getThname()] = std::numeric_limits<double>::max();
                thMax[it->getThname()] = -std::numeric_limits<double>::max();
            }
        }
    }
    kmax = k;
    kwmax = kweight;

    DefineParameters();

};

void MonteCarloEngine::setNChains(unsigned int i)
{
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
    for (std::map<std::string, BCH1D *>::iterator it = Histo1D.begin();
            it != Histo1D.end(); it++)
        delete it->second;
    for (std::map<std::string, BCH2D *>::iterator it = Histo2D.begin();
            it != Histo2D.end(); it++)
        delete it->second;
};

// ---------------------------------------------------------

void MonteCarloEngine::DefineParameters()
{
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    if (rank == 0) std::cout << "Parameters varied in Monte Carlo:" << std::endl;
    int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin();
            it < ModPars.end(); it++) {
        if (it->errf == 0. && it->errg == 0.)
            continue;
        AddParameter(it->name.c_str(), it->min, it->max);
        if (rank == 0) std::cout << "  " << it->name << " " << k << std::endl;
        DPars[it->name] = 0.;
        if (it->errf == 0.) SetPriorGauss(k, it->ave, it->errg);
        else if (it->errg == 0.) SetPriorConstant(k);
        else {
            TF1 * combined = new TF1(it->name.c_str(),
                    "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                    it->min, it->max);
            combined->SetParameter(0, it->ave);
            combined->SetParameter(1, it->errg);
            combined->SetParameter(2, it->errf);
            SetPrior(k, combined);
            delete combined;
        }
        k++;
    }
}

// ---------------------------------------------------------

double MonteCarloEngine::LogLikelihood(const std::vector<double>& parameters)
{
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.

    double logprob = 0.;

    for (unsigned int k = 0; k < parameters.size(); k++) {
        //        std::string pippo = GetParameter(k)->GetName();
        //        double pluto = parameters[k];
        //        DPars[pippo]=pluto;
        DPars[GetParameter(k)->GetName()] = parameters[k];
    }

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
    NumOfUsedEvents++;
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
    //std::cout << "logprob " << logprob <<std::endl;    
    return logprob;
}

void MonteCarloEngine::MCMCIterationInterface()
{
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
            for (unsigned int k = 0; k < pars.size(); k++) {
                DPars[GetParameter(k)->GetName()] = pars[k];
            }
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
                    if (th < thMin[it->getThname()]) thMin[it->getThname()] = th;
                    if (th > thMax[it->getThname()]) thMax[it->getThname()] = th;
                    Histo1D[it->getThname()]->GetHistogram()->Fill(th);
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
                    Histo2D[it->getThname() + "_vs_" + it->getThname2()]->GetHistogram()->Fill(th1, th2);
                }

                // fill the histograms for correlated observables
                for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                        it1 < CGO.end(); it1++) {
                    std::vector<Observable> pino(it1->getObs());
                    for (std::vector<Observable>::iterator it = pino.begin();
                            it != pino.end(); ++it) {
                        double th = buff[il][k++];

                        /* set the min and max of theory values */
                        if (th < thMin[it->getThname()]) thMin[it->getThname()] = th;
                        if (th > thMax[it->getThname()]) thMax[it->getThname()] = th;

                        Histo1D[it->getThname()]->GetHistogram()->Fill(th);
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
        for (unsigned int k = 0; k < GetNParameters(); k++) {
            //        std::string pippo = GetParameter(k)->GetName();
            //        double pluto = parameters[k];
            //        DPars[pippo]=pluto;
            DPars[GetParameter(k)->GetName()] = fMCMCx.at(i * GetNParameters() + k);
        }

        Mod->Update(DPars);

        // fill the histograms for observables
        int k = 0, kweight = 0;
        for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) {
            double th = it->computeTheoryValue();

            /* set the min and max of theory values */
            if (th < thMin[it->getThname()]) thMin[it->getThname()] = th;
            if (th > thMax[it->getThname()]) thMax[it->getThname()] = th;

            Histo1D[it->getThname()]->GetHistogram()->Fill(th);
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
            Histo2D[it->getThname() + "_vs_" + it->getThname2()]->GetHistogram()->Fill(th1, th2);
        }

        // fill the histograms for correlated observables
        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                it1 < CGO.end(); it1++) {
            std::vector<Observable> pino(it1->getObs());
            for (std::vector<Observable>::iterator it = pino.begin();
                    it != pino.end(); ++it) {
                double th = it->computeTheoryValue();

                /* set the min and max of theory values */
                if (th < thMin[it->getThname()]) thMin[it->getThname()] = th;
                if (th > thMax[it->getThname()]) thMax[it->getThname()] = th;

                Histo1D[it->getThname()]->GetHistogram()->Fill(th);
            }
        }
    }
#endif
}

void MonteCarloEngine::CheckHistogram(const TH1D& hist, const std::string name)
{
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

void MonteCarloEngine::CheckHistogram(const TH2D& hist, const std::string name)
{
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
        const std::string OutputDir)
{
    if (Histo1D[it.getThname()]->GetHistogram()->Integral() > 0.0) {
        std::string fname = OutputDir + "/" + it.getThname() + ".pdf";
        //        BCH1D* pippo =  Histo1D[it.getThname()];
        //        double x = pippo->GetMean();
        //        pippo->Print("Dmd1.pdf");
        Histo1D[it.getThname()]->SetGlobalMode(it.computeTheoryValue());
        Histo1D[it.getThname()]->Print(fname.c_str());
        std::cout << fname << " has been created." << std::endl;
        out.Write(Histo1D[it.getThname()]->GetHistogram());
        CheckHistogram(*Histo1D[it.getThname()]->GetHistogram(), it.getThname());
    } else
        HistoLog << "WARNING: The histogram of "
            << it.getThname() << " is empty!" << std::endl;

    double min = thMin[it.getThname()];
    double max = thMax[it.getThname()];
    double range = max - min;
    HistoLog.precision(10);
    HistoLog << "  [min, max]=[" << min << ", " << max << "]" << std::endl;
    HistoLog.precision(6);
}

void MonteCarloEngine::PrintHistogram(BCModelOutput & out, const std::string OutputDir)
{
    std::vector<double> mode(GetBestFitParameters());
    for (unsigned int k = 0; k < GetNParameters(); k++)
        DPars[GetParameter(k)->GetName()] = mode[k];
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
        std::string HistName = it->getThname() + "_vs_" + it->getThname2();
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

void MonteCarloEngine::AddChains()
{
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

void MonteCarloEngine::PrintCorrelationMatrix(const std::string filename)
{
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


/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBSERVABLE_H
#define	OBSERVABLE_H

#include <string>
#include <iostream>
#include "ThObservable.h"

class Observable {
public:
    Observable(const std::string name_i, const std::string thname_i,
            const std::string label_i, const bool tMCMC_i, const double min_i,
            const double max_i, ThObservable * tho_i);
    Observable(const Observable& orig);
    double getTheoryValue();
    virtual ~Observable();

    double getAve() const {
        return ave;
    }

    void setAve(double ave) {
        this->ave = ave;
    }

    std::string getDistr() const {
        return distr;
    }

    void setDistr(std::string distr) {
        this->distr = distr;
    }

    double getErrf() const {
        return errf;
    }

    void setErrf(double errf) {
        this->errf = errf;
    }

    double getErrg() const {
        return errg;
    }

    void setErrg(double errg) {
        this->errg = errg;
    }

    std::string getFilename() const {
        return filename;
    }

    void setFilename(std::string filename) {
        this->filename = filename;
    }

    std::string getHistoname() const {
        return histoname;
    }

    void setHistoname(std::string histoname) {
        this->histoname = histoname;
    }

    std::string getLabel() const {
        return label;
    }

    void setLabel(std::string label) {
        this->label = label;
    }

    double getMax() const {
        return max;
    }

    void setMax(double max) {
        this->max = max;
    }

    double getMin() const {
        return min;
    }

    void setMin(double min) {
        this->min = min;
    }

    std::string getName() const {
        return name;
    }

    void setName(std::string name) {
        this->name = name;
    }

    bool isTMCMC() const {
        return tMCMC;
    }

    void setTMCMC(bool tMCMC) {
        this->tMCMC = tMCMC;
    }

    std::string getThname() const {
        return thname;
    }

    void setThname(std::string thname) {
        this->thname = thname;
    }

    ThObservable* getTho() const {
        return tho;
    }

    void setTho(ThObservable* tho) {
        this->tho = tho;
    }

 
    friend std::ostream& operator<<(std::ostream& output, const Observable& o);

protected:
    ThObservable * tho;
    std::string name, thname, label, distr, filename, histoname;
    double ave, errg, errf, min, max;
    bool tMCMC;
};


#endif	/* OBSERVABLE_H */


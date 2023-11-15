/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OTHERTHOBSERVABLES_H
#define	OTHERTHOBSERVABLES_H

#include "ThObservable.h"

class NPbase;

//-----  Collider observables: LHC dilepton events  ----------

/**
 * @class NevLHCee13
 * @brief An observable class for the number of events at the LHC for a given process
 * @f$N_{events}^{LHC}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
 *
 */
class NevLHCee13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] i_bin_i the bin number
     */
    NevLHCee13(const StandardModel& SM_i, const int i_bin_i);
      
    /**
     * @brief Destructor of the NevLHCee13 class.
     */
    virtual ~NevLHCee13();

    /**
     * @brief The number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
     * @return @f$N_{events}^{LHC}@f$
     */
    double computeThValue();
      
private:
    const NPbase* myNPbase;
    const int i_bin;
};

//---------------------------------------------------------------

/**
 * @class NevLHCmumu13
 * @brief An observable class for the number of events at the LHC for a given process
 * @f$N_{events}^{LHC}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
 *
 */
class NevLHCmumu13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] i_bin_i the bin number
     */
    NevLHCmumu13(const StandardModel& SM_i, const int i_bin_i);
      
    /**
     * @brief Destructor of the NevLHCmumu13 class.
     */
    virtual ~NevLHCmumu13();

    /**
     * @brief The number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
     * @return @f$N_{events}^{LHC}@f$
     */
    double computeThValue();
      
private:
    const NPbase* myNPbase;
    const int i_bin;
};

//---------------------------------------------------------------

/**
 * @class NevLHCtautau13
 * @brief An observable class for the number of events at the LHC for a given process
 * @f$N_{events}^{LHC}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
 *
 */
class NevLHCtautau13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] i_bin_i the bin number
     */
    NevLHCtautau13(const StandardModel& SM_i, const int i_bin_i);
      
    /**
     * @brief Destructor of the NevLHCtautau13 class.
     */
    virtual ~NevLHCtautau13();

    /**
     * @brief The number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
     * @return @f$N_{events}^{LHC}@f$
     */
    double computeThValue();
      
private:
    const NPbase* myNPbase;
    const int i_bin;
};

//---------------------------------------------------------------


//-----  Collider observables: LHC mono-lepton events  ----------

/**
 * @class NevLHCenu13
 * @brief An observable class for the number of events at the LHC for a given process
 * @f$N_{events}^{LHC}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
 *
 */
class NevLHCenu13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] i_bin_i the bin number
     */
    NevLHCenu13(const StandardModel& SM_i, const int i_bin_i);
      
    /**
     * @brief Destructor of the NevLHCenu13 class.
     */
    virtual ~NevLHCenu13();

    /**
     * @brief The number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
     * @return @f$N_{events}^{LHC}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const int i_bin;
};

//---------------------------------------------------------------

/**
 * @class NevLHCmunu13
 * @brief An observable class for the number of events at the LHC for a given process
 * @f$N_{events}^{LHC}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
 *
 */
class NevLHCmunu13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] i_bin_i the bin number
     */
    NevLHCmunu13(const StandardModel& SM_i, const int i_bin_i);
      
    /**
     * @brief Destructor of the NevLHCmunu13 class.
     */
    virtual ~NevLHCmunu13();

    /**
     * @brief The number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
     * @return @f$N_{events}^{LHC}@f$
     */
    double computeThValue();
      
private:
    const NPbase* myNPbase;
    const int i_bin;
};

//---------------------------------------------------------------

/**
 * @class NevLHCtaunu13
 * @brief An observable class for the number of events at the LHC for a given process
 * @f$N_{events}^{LHC}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
 *
 */
class NevLHCtaunu13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] i_bin_i the bin number
     */
    NevLHCtaunu13(const StandardModel& SM_i, const int i_bin_i);
      
    /**
     * @brief Destructor of the NevLHCtaunu13 class.
     */
    virtual ~NevLHCtaunu13();

    /**
     * @brief The number of events at the LHC for a given process @f$N_{events}^{LHC}@f$.
     * @return @f$N_{events}^{LHC}@f$
     */
    double computeThValue();
      
private:
    const NPbase* myNPbase;
    const int i_bin;
};

//---------------------------------------------------------------




#endif	/* OTHERTHOBSERVABLES_H */


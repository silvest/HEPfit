/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ORDERSCHEME_H
#define	ORDERSCHEME_H

#define MAXORDER NNLO
#define MAXORDER_EW NLO_ew

/**
 * @enum schemes
 * @ingroup StandardModel
 * @brief An enum type for regularization schemes.
 */
enum schemes
{
    NDR=0, /**< Naive dimensional regularization (NDR) scheme */
    HV, /**< 't Hooft-Veltman (HV) scheme */
    LRI /**< Regularization-Independent (RI) renormalization schemes with the Landau gauge */
};

/**
 * @enum orders
 * @ingroup StandardModel
 * @brief An enum type for orders in %QCD.
 */
enum orders
{
    LO=0, /**< Leading order */
    NLO,  /**< Next-to-leading order */
    NNLO,  /**< Next-to-next-to-leading order */
    FULLNLO, /**< Full NLO = LO + NLO */
    FULLNNLO /**< Full NNLO = LO + NLO + NNLO */
};

/**
 * @enum orders_ew
 * @ingroup StandardModel
 * @brief An enum type for orders in electroweak.
 */
enum orders_ew
{
    NULL_ew = orders(MAXORDER+1), /**< An auxiliary enumerator */
    LO_ew, /**< Leading order */
    NLO_ew  /**< Next-to-leading order */
};

#endif	/* ORDERSCHEME_H */
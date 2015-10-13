/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ORDERSCHEME_H
#define	ORDERSCHEME_H

#define MAXORDER FULLNNLO
#define MAXORDER_EW FULLNLO_ew

/**
 * @enum schemes
 * @ingroup StandardModel
 * @brief An enum type for regularization schemes.
 */
enum schemes
{
    NDR = 0, /**< Naive dimensional regularization (NDR) scheme */
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
    LO = 0, /**< Leading order */
    NLO, /**< Next-to-leading order */
    NNLO, /**< Next-to-next-to-leading order */
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
    NULL_ew = orders(MAXORDER + 1), /**< An auxiliary enumerator */
    LO_ew, /**< Leading order */
    NLO_ew, /**< Next-to-leading order */
    NLO_ewt1,   /* e^2/s^2 */
    NLO_ewt2,     /* es */
    NLO_ewt3,    /* e^2/s */
    NLO_ewt4,     /* e^2 */
    FULLNLO_ew /**< Full NLO_ew = LO + NLO + LO_ew + NLO_ew +...+ NLO_ewt4 */
};

#endif	/* ORDERSCHEME_H */
/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ORDERSCHEME_H
#define	ORDERSCHEME_H

#define MAXORDER FULLNNNLO
#define MAXORDER_QED FULLNLO_QED

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
    NNNLO, /**< Next-to-next-to-next-to-leading order */
    FULLNLO, /**< Full NLO = LO + NLO */
    FULLNNLO, /**< Full NNLO = LO + NLO + NNLO */
    FULLNNNLO /**< Full NNNLO = LO + NLO + NNLO + NNNLO */        
};

/**
 * @enum orders_qed
 * @ingroup StandardModel
 * @brief An enum type for orders in electroweak.
 */
enum orders_qed // WARNING: don't change the ordering, it matters in HeffDF1
{
    NO_QED = orders(MAXORDER) + 1, /**< An auxiliary enumerator */
    LO_QED, /**< Leading order e/s */
    NLO_QED11, /**< Next-to-leading order e */
    NLO_QED21,     /* e*s */
    NLO_QED02,   /* e^2/s^2 */
    NLO_QED12,    /* e^2/s */
    NLO_QED22,     /* e^2 */
    FULLNLO_QED /**< Full NLO_QED = LO + NLO + LO_QED + NLO_QED +...+ NLO_QED22 */
};

// New enum for orders introduced with Expanded
/**
 * @enum qcd_orders
 * @ingroup StandardModel
 * @brief An enum type for qcd_orders in %QCD.
 */
enum qcd_orders
{
    QCD0 = 0, /**< Leading order */
    QCD1, /**< Next-to-leading order */
    QCD2, /**< Next-to-next-to-leading order */
    QCD3, /**< Next-to-next-to-next-to-leading order */
    FULLQCD1, /**< Full NLO = LO + NLO */
    FULLQCD2, /**< Full NNLO = LO + NLO + NNLO */
    FULLQCD3 /**< Full NNLO = LO + NLO + NNLO + NNNLO */        
};

/**
 * @enum qed_orders
 * @ingroup StandardModel
 * @brief An enum type for qed_orders in electroweak.
 */
enum qed_orders // WARNING: don't change the ordering, it matters in HeffDF1
{
    QED0=0, /* Leading order e/s */
    QED1, /* */
    QED2, /**< Next-to-leading order e */
    FULLQED1,
    FULLQED2 /* all terms up to QED2 included */
};

#endif	/* ORDERSCHEME_H */
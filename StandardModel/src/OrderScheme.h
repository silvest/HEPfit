/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
 * @brief Available schemes are NDR, HV and LRI.
 */
enum schemes {NDR=0, HV, LRI};

/**
 * @enum orders
 * @ingroup StandardModel
 * @brief Available orders in QCD are LO, NLO, NNLO, FULLNLO = LO + NLO , FULLNNLO = LO + NLO + NNLO.
 */
enum orders {LO=0, NLO, NNLO, FULLNLO, FULLNNLO};

/**
 * @enum orders_ew
 * @ingroup StandardModel
 * @brief Available orders in Electroweak are NULL_ew, LO_ew and NLO_ew.
 */
enum orders_ew {NULL_ew = orders(MAXORDER+1), LO_ew, NLO_ew};

#endif	/* ORDERSCHEME_H */
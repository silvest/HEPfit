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
 */
enum schemes {NDR=0, HV, LRI};

/**
 * @enum orders
 * @ingroup StandardModel
 */
enum orders {LO=0, NLO, NNLO, FULLNLO, FULLNNLO};

/**
 * @enum orders_ew
 * @ingroup StandardModel
 */
enum orders_ew {NULL_ew = orders(MAXORDER+1), LO_ew, NLO_ew};

#endif	/* ORDERSCHEME_H */


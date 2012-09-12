/* 
 * File:   orderscheme.h
 * Author: enrico
 *
 * Created on May 12, 2011, 11:19 AM
 */

#ifndef ORDERSCHEME_H
#define	ORDERSCHEME_H

#define MAXORDER NNLO
#define MAXORDER_EW NLO_ew

enum schemes {NDR=0, HV, LRI};
enum orders {LO=0, NLO, FULLNLO, NNLO, FULLNNLO};
enum orders_ew {NULL_ew = orders(MAXORDER+1), LO_ew, NLO_ew};

#endif	/* ORDERSCHEME_H */


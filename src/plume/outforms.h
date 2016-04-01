/*
 * outforms.h
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#ifndef OUTFORMS_H
#define OUTFORMS_H
#include <hashtab.h>
#include <buckhead.h>
#include <properties.h>

void write_h5part (int, int, THashTable *, TimeProps *);

void write_h5part_show (int, int, THashTable *, TimeProps *);

void write_matlab (int, THashTable *, HashTable *, TimeProps *, vector<BucketHead> &);


#endif /* OUTFORMS_H_ */

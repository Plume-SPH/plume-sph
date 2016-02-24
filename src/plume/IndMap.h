/*
 * IndMap.h
 *
 *  Created on: Mar 31, 2015
 *      Author: zhixuanc
 */

#ifndef INDMAP_H
#define INDMAP_H


#include <hashtab.h>
#include <thashtab.h>
#include "constant.h"

struct IndMap
{
// friend bool sortByIndex(const IndMap&, const IndMap&);
// friend bool MyidIs::operator()(const Initial_Index_Map&);
// private:
	unsigned key[TKEYLENGTH];//particle key
	int hash;
	int index;
	int p_id; //processor id to which the particle belong to
// public:
	//! Constructors
	IndMap (TKey *keyin, int hashin, int idin)
	{
		int i;
		for (i = 0; i < TKEYLENGTH; i++)
			key [i] = keyin->key[i];

		hash=hashin;
		p_id=idin;
		index=-1;
	}

	//default constructor
	IndMap ()
	{
	 int i;
	 for (i = 0; i < TKEYLENGTH; i++)
			key [i] = 0;

		hash=0;
		p_id=0;
		index=-1;
	}
    //
	void put_key (TKey* keyin)
	{
		int i;
		for (i = 0; i <TKEYLENGTH; i++)
			key [i] = keyin->key[i];
	}
    //
	void put_key (const unsigned* keyin)
	{
		int i;
		for (i = 0; i <TKEYLENGTH; i++)
			key [i] = *(keyin + i);
	}
	//
	void put_key (unsigned keyin, int i)
	{
		key [i] = keyin;
	}
	//
	void put_hash (int hashin)
	{
		hash=hashin;
	}
	//
	void put_index (int indexin)
	{
		index=indexin;
	}
	//
	void put_pid (int pidin)
	{
		p_id=pidin;
	}
	//
	const unsigned* get_key() const
	{
		return key;
	}

	int get_hash () const
	{
		return hash;
	}
	//
	int get_index () const
	{
		return index;
	}

	int get_pid() const
	{
		return p_id;
	}
};

struct Initial_Index_Map
{
	int id;
	int iidex;

	Initial_Index_Map (int idin, int indexin)
	{
		id =idin;
		iidex=indexin;
	}

	Initial_Index_Map ()
	{
		return;
	}

	int get_id () const
	{
		return id;
	}

	int get_iidex () const
	{
		return iidex;
	}

	void put_id (int idin)
	{
		id=idin;
	}

	void put_iidex (int idexin)
	{
		iidex=idexin;
	}
};

#endif /* INDMAP_H_ */

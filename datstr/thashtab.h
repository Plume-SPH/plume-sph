/*
 * thashtab.h
 *
 *  Created on: May 5, 2015
 *      Author: zhixuanc
 */

#ifndef THASHTAB_H
#define THASHTAB_H

#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

#include <constant.h>

struct TKey
{
  unsigned key[TKEYLENGTH];

  // default constructor
  TKey ()
  {
    int i;
    for (i=0; i<TKEYLENGTH; i++)
    	key[i]=0 ;
  }

  // constructor
  TKey (unsigned keyin[])
  {
    for (int i = 0; i < TKEYLENGTH;  i++)
      key[i] = keyin[i];
  }

  // copy constructor
  TKey (const TKey & rhs)
  {
    int i;
    for (i = 0; i < TKEYLENGTH; i++)
        key[i] = rhs.key[i];
  }

  // Key to unsigned-array
  void fill_key (unsigned k[])
  {
    for (int i = 0; i < TKEYLENGTH; i++)
      k[i] = key[i];
  }

  // overload assignment
  TKey & operator= (const TKey & rhs)
  {
    for (int i = 0; i < TKEYLENGTH; i++)
      key[i] = rhs.key[i];
    return * this;
  }

  // overload equality
  bool operator== (const TKey & rhs) const
  {
    for (int i = 0; i < TKEYLENGTH; i++)
      if ( key[i] != rhs.key[i] )
        return false;
    return true;
  }

  // overload less than TKEYLENGTH = 3;
  bool operator< (const TKey & rhs) const
  {
    if (key[0] < rhs.key[0])
      return true;
    else if (key[0] > rhs.key[0])
      return false;
    else if (key[1] < rhs.key[1])
      return true;
    else if (key[1] > rhs.key[1])
      return false;
    else if (key[2] < rhs.key[2])
      return true;
    else
      return false;
  }
};

// KEYLENGTH = 3;
// return true if k1<k2
inline bool
compare_keys (const TKey & K1, const TKey & K2)
{

  if (K1.key[0] < K2.key[0])
    return true;
  else if (K1.key[0] > K2.key[0])
    return false;
  else if (K1.key[1] < K2.key[1])
	return true;
  else if (K1.key[1] > K2.key[1])
    return false;
  else if (K1.key[2] < K2.key[2])
	return true;
  else
    return false;
}

struct THashEntry
{
  unsigned key[TKEYLENGTH];      // key: object key word
  void * value;                 // value: poiter to record
  THashEntry * pre;              // pre, next: objects with same entry
  THashEntry * next;             // will be stored in a two-way linked-list
  //int N;                        //The total number of particles that was added 

  THashEntry (unsigned * keyi)
  {
    int i;
    for (i = 0; i < TKEYLENGTH; i++)
        key[i] = keyi[i];
      next = NULL;
      pre = NULL;
    //N=0;
  }

  THashEntry ()
  {
    value = NULL;
    next = NULL;
    pre = NULL;
    //N = 0;
  }

  ~THashEntry ()                 //keep the follower when deleting an object
  {
    if (next)
      next->pre = pre;
    if (pre)
      pre->next = next;
  }
};

typedef THashEntry *THashEntryPtr;

class THashTable
{

protected:
  unsigned Range;
  double umax;
  double minDom[DIMENSION];
  double maxDom[DIMENSION];

  THashEntryPtr *bucket;//The name is misleading, it should be particles, too lazy to change.
  int NBUCKETS; //also too lazy to change
  int PRIME;
//  int ENTRIES;
  int SIZE01;
  int SIZE02;
//  int MAX_ADD;
  THashEntryPtr addElement (int entry, unsigned *key);
  THashEntryPtr searchBucket (THashEntryPtr p, unsigned *key);

public:
  THashTable (int,  int, double *, double *); //constructor that I will use.
//  THashTable (int,  int, double *, double *, int); //overloading constructor with maximum number of material adding included.
  ~THashTable ();

  void add (unsigned * key, void * value);
  void *lookup (unsigned * key);
  void *lookup (TKey); // not sure about this one? The syntax is incorrect??
  void remove (unsigned * key);
  int hash (unsigned * key);

  // remove entery
  void remove (TKey k)
  {
    unsigned key[TKEYLENGTH];
      k.fill_key (key);
      remove (key);
  }

  //! add entry to the hash-table
  void add (TKey inkey, void * value)
  {
    unsigned key[TKEYLENGTH];
    inkey.fill_key (key);
    add (key, value);
  }

  //! get size of the hash-table
  int get_no_of_buckets ()
  {
    return NBUCKETS;
  }

  //! get first entry
  THashEntryPtr * getbucketptr ()
  {
    return bucket;
  }

  //! get i-th bucket
  THashEntryPtr getBucket (int entry)
  {
    return bucket[entry];
  }

  //! get min of domain
  double *get_minDom ()
  {
    return minDom;
  }

  //! get max of domain
  double *get_maxDom ()
  {
    return maxDom;
  }

};


/*********************************
 *  Hash-Table iterator
 ********************************/
class THTIterator
{
private:
  THashTable * table;
  THashEntryPtr current;
  int index;
  int size;

public:
   ~THTIterator ()
  {
  };

  THTIterator (THashTable * ht);
  THashEntryPtr getNextBucket ();
  void *next ();
  void reset ()
  {
    current = table->getBucket (0);
    index = 0;
  };
};

#endif /* THASHTAB_H_ */

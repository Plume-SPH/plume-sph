/*
 * hashtab.h
 *
 *  Created on: Mar 5, 2015
 *      Author: zhixuanc
 */

#ifndef HASHTAB_H
#define HASHTAB_H

#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

#include <constant.h>

struct Key
{
  unsigned key[KEYLENGTH];

  // default constructor
  Key ()
  {
    int i;
    for (i=0; i<KEYLENGTH; i++)
    	key[i]=0 ;
  }

  // constructor
  Key (unsigned keyin[])
  {
    for (int i = 0; i < KEYLENGTH;  i++)
      key[i] = keyin[i];
  }

  // copy constructor
  Key (const Key & rhs)
  {
    int i;
    for (i = 0; i < KEYLENGTH; i++)
        key[i] = rhs.key[i];
  }

  // Key to unsigned-array
  void fill_key (unsigned k[])
  {
    for (int i = 0; i < KEYLENGTH; i++)
      k[i] = key[i];
  }

  // overload assignment
  Key & operator= (const Key & rhs)
  {
    for (int i = 0; i < KEYLENGTH; i++)
      key[i] = rhs.key[i];
    return * this;
  }

  // overload equality
  bool operator== (const Key & rhs) const
  {
    for (int i = 0; i < KEYLENGTH; i++)
      if ( key[i] != rhs.key[i] )
        return false;
    return true;
  }

  // overload less than
  bool operator< (const Key & rhs) const
  {
    if (key[0] < rhs.key[0])
      return true;
    else if (key[0] > rhs.key[0])
      return false;
    else if (key[1] < rhs.key[1])
      return true;
    else
      return false;
  }
};

inline bool
compare_keys (const Key & K1, const Key & K2)
{

  if (K1.key[0] < K2.key[0])
    return true;
  else if (K1.key[0] > K2.key[0])
    return false;
  else if (K1.key[1] < K2.key[1])
    return true;
  else
    return false;
}

struct HashEntry
{
  unsigned key[KEYLENGTH];      // key: object key word
  void * value;                 // value: poiter to record
  HashEntry * pre;              // pre, next: objects with same entry
  HashEntry * next;             // will be stored in a two-way linked-list

  HashEntry (unsigned * keyi)
  {
    int i;
    for (i = 0; i < KEYLENGTH; i++)
        key[i] = keyi[i];
      next = NULL;
      pre = NULL;
  }

  HashEntry ()
  {
    value = NULL;
    next = NULL;
    pre = NULL;
  }

  ~HashEntry ()                 //keep the follower when deleting an object
  {
    if (next)
      next->pre = pre;
    if (pre)
      pre->next = next;
  }
};

typedef HashEntry *HashEntryPtr;

class HashTable
{

protected:
  unsigned Range;
  double umax;
  double minDom[DIMENSION];
  double maxDom[DIMENSION];

  HashEntryPtr *bucket;
  int NBUCKETS;
  int PRIME;
//  int ENTRIES;
  int SIZE01;
  int SIZE02;

  HashEntryPtr addElement (int entry, unsigned *key);
  HashEntryPtr searchBucket (HashEntryPtr p, unsigned *key);

public:
  HashTable (int,  int, double *, double *);
  ~HashTable ();

  void add (unsigned * key, void * value);
  void *lookup (unsigned * key);
  void *lookup (Key);
  void remove (unsigned * key);
  int hash (unsigned * key);

  // remove entery
  void remove (Key k)
  {
    unsigned key[KEYLENGTH];
      k.fill_key (key);
      remove (key);
  }

  //! add entry to the hash-table
  void add (Key inkey, void * value)
  {
    unsigned key[KEYLENGTH];
    inkey.fill_key (key);
    add (key, value);
  }

  //! get size of the hash-table
  int get_no_of_buckets ()
  {
    return NBUCKETS;
  }

  //! get first entry
  HashEntryPtr * getbucketptr ()
  {
    return bucket;
  }

  //! get i-th bucket
  HashEntryPtr getBucket (int entry)
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
class HTIterator
{
private:
  HashTable * table;
  HashEntryPtr current;
  int index;
  int size;

public:
   ~HTIterator ()
  {
  };

  HTIterator (HashTable * ht);
  HashEntryPtr getNextBucket ();
  void *next ();
  void reset ()
  {
	  current = table->getBucket (0);
	  index = 0;
  };
};

#endif /* HASHTAB_H_ */

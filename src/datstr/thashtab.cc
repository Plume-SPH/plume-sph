/*
 * thashtab.cc
 *
 *  Created on: May 5, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cassert>
#include <climits>
using namespace std;

#include "thashtab.h"

#define THASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

THashTable::THashTable (int size, int prime, double minR[], double maxR[])
{
  int i;
  NBUCKETS = size;
  PRIME = prime;
  SIZE02 = NBUCKETS / 10;
  SIZE01 = NBUCKETS - SIZE02;
  umax = (double) IScale;
//  MAX_ADD = 1000;

  // allocate table-size
  bucket = new THashEntryPtr [NBUCKETS];
  for (i = 0; i < NBUCKETS; i++)
    bucket[i] = NULL;

  for (i = 0; i < DIMENSION; i++)
  {
    minDom[i] = minR[i];
    maxDom[i] = maxR[i];
  }
}

THashTable::~THashTable ()        //evacuate the table
{
  for (int i = 0; i < NBUCKETS; i++)
  {
    THashEntryPtr p = bucket[i];

    while (p)
    {
      THashEntryPtr p_next = p->next;
      delete p;

      p = p_next;
    }
  }
  delete[]bucket;
}

THashEntryPtr THashTable::searchBucket (THashEntryPtr p, unsigned *keyi)
{
  int
    i;

  while (p)
  {
    for (i = 0; i < TKEYLENGTH; i++)
    {
      if (p->key[i] != *(keyi + i))
      {
        p = p->next;
        break;                  //not found, check next element
      }
      else if (i == TKEYLENGTH - 1)
        return p;               //found, return element pointer
    }
  }
  return p;
}

THashEntryPtr THashTable::addElement (int entry, unsigned key[])
{
  THashEntryPtr
    p = new THashEntry (key);

  if ((bucket[entry]))          //this place is already occupied
  {
    THashEntryPtr
      currentPtr = bucket[entry];

    while (currentPtr != 0 && (key[0] > currentPtr->key[0]))
    {
      p->pre = currentPtr;
      currentPtr = currentPtr->next;
    }

    if (currentPtr != 0 && key[0] == currentPtr->key[0])
    {
      while (currentPtr != 0 && (key[1] > currentPtr->key[1]))
      {
        p->pre = currentPtr;
        currentPtr = currentPtr->next;
      }
    }

    if (currentPtr != 0 && key[0] == currentPtr->key[0] && key[1] == currentPtr->key[1])
    {
      while (currentPtr != 0 && (key[2] > currentPtr->key[2]))
      {
        p->pre = currentPtr;
        currentPtr = currentPtr->next;
      }
    }

    if (currentPtr)
      currentPtr->pre = p;
    p->next = currentPtr;
    currentPtr = p->pre;
    if (currentPtr)
      currentPtr->next = p;
    else
      bucket[entry] = p;//if previous previous does not exist, that means this should be the first element in the hash table
  }

  //  p->next = *(bucket+entry);        //add the bucket to the head
  else
    bucket[entry] = p;          //else eliminate it
  return p;
}

void *
THashTable::lookup (unsigned * key)
{
  int entry = hash (key);
  THashEntryPtr p = searchBucket (bucket[entry], key);// not understand...

  if (!p)
    return NULL;                //if not found, return 0
  return p->value;              //if found return a pointer
}

// lookup using key structure
void *
THashTable::lookup (TKey kstr)
{
  unsigned key[TKEYLENGTH];

  kstr.fill_key (key);
  int entry = hash (key);
  THashEntryPtr p = searchBucket (bucket[entry], key);

  if (!p)
    return NULL;
  return p->value;
}

void
THashTable::add (unsigned *key, void *value)
{
  int entry = hash (key);
  THashEntryPtr p = searchBucket (bucket[entry], key);

  if (p == NULL)    //make sure that the same HashEntryPtr does not exist in the hashtable
  {
    p = addElement (entry, key);
    p->value = value;
  }
  return;
}

void
THashTable::remove (unsigned *key)
{
  int entry = hash (key);
  THashEntryPtr p = searchBucket (bucket[entry], key);

  if (p == NULL)
    return;
  if (bucket[entry] == p)
  {
    bucket[entry] = p->next;
    delete p;
  }
  else
  {
    if (!(p->next))
      delete p;

    else
    {
      (p->pre)->next = p->next;
      (p->next)->pre = p->pre;
      delete p;
    }
  }
}

#if CODE_DIMENSION==1
int THashTable::hash (unsigned * key)  //A mask hash table
{
  return (key[0]);
}
#else
int THashTable::hash (unsigned * key)
{
//  int S03 = MAX_ADD;
//  int S02 = (NBUCKETS-S03)/10;
//  int S01 = NBUCKETS-S03-S02;
  int S02 = NBUCKETS/10;
  int S01 = NBUCKETS-S02;
  int i1 = (int) ((double) key[0] / umax * S01);
  int i2 = (int) ((double) key[1] / umax * S02);
//  int i3 = key[2];

//  return (i1 + i2 +i3);
  return (i1 + i2);
}
#endif


/***************************************************
 *          Hashtable iterator
 ***************************************************/
THTIterator::THTIterator (THashTable * ht)
{
  table = ht;
  current = table->getBucket (0);
  size = table->get_no_of_buckets ();
  index = 0;
}

THashEntryPtr THTIterator::getNextBucket ()
{
  while (table->getBucket (++index) == NULL && index < size);
  if (index == size)
    return NULL;
  return table->getBucket (index);
}

void *
THTIterator::next ()
{
  void *value;

  if (current)                  // if current pointer has more links
  {
    value = current->value;
    current = current->next;    // move to next link, after getting the value
    return value;
  }
  current = getNextBucket ();   // get next valid Hashtable entry
  if (current == NULL)
    return NULL;
  value = current->value;
  current = current->next;      // move to next link, after getting the value
  return value;
}

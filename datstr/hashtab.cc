/*
 * hashtab.cc
 *
 *  Created on: Mar 5, 2015
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

#include "hashtab.h"

#define HASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

HashTable::HashTable (int size, int prime, double minR[], double maxR[])
{
  int i;
  NBUCKETS = size;
  PRIME = prime;
  SIZE02 = NBUCKETS / 10;
  SIZE01 = NBUCKETS - SIZE02;
  umax = (double) IScale;


  // allocate table-size
  bucket = new HashEntryPtr [NBUCKETS];
  for (i = 0; i < NBUCKETS; i++)
    bucket[i] = NULL;

  for (i = 0; i < DIMENSION; i++)
  {
    minDom[i] = minR[i];
    maxDom[i] = maxR[i];
  }
}

HashTable::~HashTable ()        //evacuate the table
{
  for (int i = 0; i < NBUCKETS; i++)
  {
    HashEntryPtr p = bucket[i];

    while (p)
    {
      HashEntryPtr p_next = p->next;
      delete p;

      p = p_next;
    }
  }
  delete[]bucket;
}

HashEntryPtr HashTable::searchBucket (HashEntryPtr p, unsigned *keyi)
{
  int
    i;

  while (p)
  {
    for (i = 0; i < KEYLENGTH; i++)
    {
      if (p->key[i] != *(keyi + i))
      {
        p = p->next;
        break;                  //not found, check next element
      }
      else if (i == KEYLENGTH - 1)
        return p;               //found, return element pointer
    }
  }
  return p;
}

HashEntryPtr HashTable::addElement (int entry, unsigned key[])
{
  HashEntryPtr
    p = new HashEntry (key);

  if ((bucket[entry]))          //this place is already occupied
  {
    HashEntryPtr
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

    if (currentPtr)
      currentPtr->pre = p;
    p->next = currentPtr;
    currentPtr = p->pre;
    if (currentPtr)
      currentPtr->next = p;
    else
      bucket[entry] = p;
  }

  //  p->next = *(bucket+entry);        //add the bucket to the head
  else
    bucket[entry] = p;          //else eliminate it
  return p;
}

void *
HashTable::lookup (unsigned * key)
{
  int entry = hash (key);
  HashEntryPtr p = searchBucket (bucket[entry], key);

  if (!p)
    return NULL;                //if not found, return 0
  return p->value;              //if found return a pointer
}

// lookup using key structure
void *
HashTable::lookup (Key kstr)
{
  unsigned key[KEYLENGTH];

  kstr.fill_key (key);
  int entry = hash (key);
  HashEntryPtr p = searchBucket (bucket[entry], key);

  if (!p)
    return NULL;
  return p->value;
}

void
HashTable::add (unsigned *key, void *value)
{
  int entry = hash (key);
  HashEntryPtr p = searchBucket (bucket[entry], key);

  if (p == NULL)                //was (!p)
  {
    p = addElement (entry, key);
    p->value = value;
  }
  return;
}

void
HashTable::remove (unsigned *key)
{
  int entry = hash (key);
  HashEntryPtr p = searchBucket (bucket[entry], key);

  if (!p)
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

// hash-function
int HashTable::hash (unsigned * key)
{
  int i1 = (int) ((double) key[0] / umax * SIZE01);
  int i2 = (int) ((double) key[1] / umax * SIZE02) ;
  return (i1 + i2);
}

/***************************************************
 *          Hashtable iterator
 ***************************************************/
HTIterator::HTIterator (HashTable * ht)
{
  table = ht;
  current = table->getBucket (0);
  size = table->get_no_of_buckets ();
  index = 0;
}

HashEntryPtr HTIterator::getNextBucket ()
{
  while (table->getBucket (++index) == NULL && index < size);
  if (index == size)
    return NULL;
  return table->getBucket (index);
}

void *
HTIterator::next ()
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

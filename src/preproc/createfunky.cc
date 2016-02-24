/*
 * createfunky.cc
 *
 *  Created on: Mar 6, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_HDF5_H
# include <hdf5.h>
# include <hdf5calls.h>
#endif

#include <iomanip>
#include <sstream>
#include <string>
using namespace std;

#include "buckstr.h"
#include "preprocess.h"

string create_filename (string base, string ext, const int padding, int myid)
{
  ostringstream ss;
  ss << base << setw(padding) << setfill('0') << myid << ext;
  return ss.str();
}

void createfunky(int myid, int nhtvars, double * htvars,
                 vector <BucketStruct> & bg, vector <unsigned> & partition)
{

//  int i, j;
  string basename = "funky";
  string exten    = ".h5";
  const int padding = 4;

  // create filenames
  string fname = create_filename (basename, exten, padding, myid);
  hid_t fp = GH5_fopen_serial (fname.c_str(),'w');

  // write HASH TABLE constants
  int dims1[2]={nhtvars,0}; //nhtvars is the number of elements in htvars
  GH5_WriteS (fp ,"/hashtable_constants", dims1, (void *) htvars, 0, 0, DOUBLETYPE);

  // copy vector to an regular array
  int nbuck = (int) bg.size();
  BucketStruct *myBucks = new BucketStruct [nbuck];
  copy(bg.begin(), bg.end(), myBucks);

    // write Background Grid data
  GH5_write_grid_data(fp,"/Buckets", nbuck, myBucks);

  int numpart = (int) partition.size ();
  unsigned * part_table = new unsigned [numpart];
  copy (partition.begin (), partition.end (), part_table);
  dims1[0] = numpart;
  dims1[1] = 0;
  GH5_WriteS (fp, "/partition_table", dims1, (void *) part_table, 0, 0 , UINTTYPE);

  delete [] myBucks;
  delete [] part_table;

  // close file
  GH5_fclose(fp);

  return;
}

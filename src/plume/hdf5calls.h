/*
 * hdf5calls.h
 *
 *  Created on: Mar 7, 2015
 *      Author: zhixuanc
 */

#ifndef HDF5CALLS_H
#define HDF5CALLS_H

#include <hdf5.h>
#include <buckstr.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef PARALLEL_IO
# define GH5_fopen GH5_fopen_parallel
# define GH5_Write GH5_WriteP
#else
# define GH5_fopen GH5_fopen_serial
# define GH5_Write GH5_WriteS
#endif

// datatype flags
const int INTTYPE=1;
const int UINTTYPE=2;
const int FLOATTYPE=3;
const int DOUBLETYPE=4;
const int CHARTYPE=5;

hid_t GH5_fopen_serial (const char *filename, char mode);
hid_t GH5_fopen_parallel (const char *filename, char mode);

inline herr_t GH5_fclose(hid_t fp) { return(H5Fclose(fp)); }
inline herr_t GH5_gclose(hid_t id) { return(H5Gclose(id)); }

int GH5_getsize (hid_t fp, const char *fullpath, hsize_t *dims);
hid_t GH5_gopen (hid_t fp, const char *name, char mode);
hid_t GH5_createdataset (hid_t gid, hid_t dataspace, const char *dsetname, unsigned type);
hid_t GH5_cellstruct ();

// Overload data read functions
herr_t GH5_readdata(hid_t fp, const char *fullpath, double *);
herr_t GH5_readdata(hid_t fp, const char *fullpath, int *);
herr_t GH5_readdata(hid_t fp, const char *fullpath, unsigned *);

// read grid data
herr_t GH5_read_grid_data(hid_t fp, const char *path, BucketStruct *c);

// wirte grid data
herr_t GH5_write_grid_data (hid_t fp, const char *path, int size, BucketStruct *c);

// wirte serial data
herr_t GH5_WriteS (hid_t , const char *, int *, void *, int, int, const int);

// write parallel data
herr_t GH5_WriteP (hid_t , const char *, int *, void *, int, int, const int);

#endif /* HDF5CALLS_H_ */

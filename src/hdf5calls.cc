/*
 * hdf5calls.cc
 *
 *  Created on: Mar 7, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <hdf5.h>
#include <stdlib.h>
#include <assert.h>

#include "hdf5calls.h"
#include "buckstr.h"

const int HDF5_FAIL = -1;

hid_t
GH5_fopen_serial(const char *filename, char mode)
{
  hid_t fp;

  switch (mode)
  {
  case 'w':
    fp = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (!fp)
    {
      fprintf(stderr, "Unable to open new HDF file.\n");
      exit(0);
    }
    break;
    // if the file already exists
  case 'a':
    fp = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    if (fp < 0)
    {
      fprintf(stderr, "Unable to open %s.\
              Make sure its in the current directory\n", filename);
      exit(0);
    }
    break;
  case 'r':
    fp = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fp < 0)
    {
      fprintf(stderr, "Unable to open %s.\
              Make sure its in the current directory\n", filename);
      exit(0);
    }
    break;
  default:
    fprintf(stderr, "GH5 ERROR: Unknown file access option, EXITING\n");
    exit(1);
  }
  return fp;
}

#ifdef PARALLEL_IO
hid_t
GH5_fopen_parallel(const char *filename, char mode)
{
  hid_t fp;
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  switch (mode)
  {
  case 'w':
    // Create a new file collectively
    fp = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    if (!fp)
    {
      fprintf(stderr, "Unable to open new HDF file.\n");
      exit(0);
    }
    break;
  case 'a':
    fp = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    if (fp < 0)
    {
      fprintf(stderr, "Unable to open %s.\
              Make sure its in the current directory\n", filename);
      exit(0);
    }
    break;
  case 'r':
    fp = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);
    if (fp < 0)
    {
      fprintf(stderr, "Unable to open %s.\
              Make sure its in the current directory\n", filename);
      exit(0);
    }
    break;
  default:
    fprintf(stderr, "GH5 ERROR: Unknown file access option, EXITING\n");
    exit(1);
  }
  return fp;
}
#endif

int
GH5_getsize(hid_t fp, const char *fullpath, hsize_t * dims)
{
  int rank;

  // open the data set for access
#ifdef H5_USE_16_API
  hid_t dataset_id = H5Dopen(fp, fullpath);
#else
  hid_t dataset_id = H5Dopen(fp, fullpath, H5P_DEFAULT);
#endif
  // obtain access to datasapce
  hid_t dspace_id = H5Dget_space(dataset_id);

  // get the dimensions
  rank = H5Sget_simple_extent_dims(dspace_id, dims, NULL);

  // close dataset before leaving
  H5Dclose(dataset_id);
  return rank;
}

// Create/open a group within file
hid_t
GH5_gopen(hid_t fp, const char *name, char mode)
{
  hid_t group_id;

  switch (mode)
  {
  case 'w':
#ifdef  H5_USE_16_API
    group_id = H5Gcreate(fp, name, 0);
#else
    group_id = H5Gcreate(fp, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    if (group_id < 0)
    {
      fprintf(stderr, "GH5 ERROR: Failed to create new HDF5 Group: %s\n", name);
      exit(1);
    }
    break;
  case 'r':
#ifdef H5_USE_16_API
    group_id = H5Gopen(fp, name);
#else
    group_id = H5Gopen(fp, name, H5P_DEFAULT);
#endif
    if (group_id < 0)
    {
      fprintf(stderr, "GH5 ERROR: Failed to open HDF5 Group: %s\n", name);
      exit(1);
    }
    break;
    // this should not happen if life is running usual
  }
  return group_id;
}

// New HDF5 data type for buckets
hid_t
GH5_bucketstruct()
{
  hid_t newtype;

  // declare a Struct prototype
  BucketStruct buck;

  // create hdf5 compound dataype
  newtype = H5Tcreate(H5T_COMPOUND, sizeof(buck));

  // insert all members one-by-one
  H5Tinsert(newtype, "Key", HOFFSET(BucketStruct, key), H5T_NATIVE_UINT);
  H5Tinsert(newtype, "neighs", HOFFSET(BucketStruct, neighs), H5T_NATIVE_UINT);
  H5Tinsert(newtype, "buckettype", HOFFSET(BucketStruct, buckettype),
            H5T_NATIVE_INT);
  H5Tinsert(newtype, "myproc", HOFFSET(BucketStruct, myproc), H5T_NATIVE_INT);
  H5Tinsert(newtype, "neigh_proc", HOFFSET(BucketStruct, neigh_proc),
            H5T_NATIVE_INT);
  H5Tinsert(newtype, "xcoord", HOFFSET(BucketStruct, xcoord),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(newtype, "ycoord", HOFFSET(BucketStruct, ycoord),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(newtype, "zcoord", HOFFSET(BucketStruct, zcoord),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(newtype, "elev", HOFFSET(BucketStruct, elev), H5T_NATIVE_DOUBLE);

  // return identifier
  return newtype;
}

herr_t
GH5_readdata(hid_t fp, const char *fullpath, double *buf)
{
  herr_t status;

  // open dataset for reading
#ifdef H5_USE_16_API
  hid_t data_id = H5Dopen(fp, fullpath);
#else
  hid_t data_id = H5Dopen(fp, fullpath, H5P_DEFAULT);
#endif

  // read the data
  status =
    H5Dread(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

  // close dataset
  H5Dclose(data_id);
  return status;
}

herr_t
GH5_readdata(hid_t fp, const char *fullpath, int *buf)
{
  herr_t status;

  // open dataset for reading
#ifdef H5_USE_16_API
  hid_t data_id = H5Dopen(fp, fullpath);
#else
  hid_t data_id = H5Dopen(fp, fullpath, H5P_DEFAULT);
#endif
  // read the data
  status = H5Dread(data_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

  // close dataset
  H5Dclose(data_id);
  return status;
}

herr_t
GH5_readdata(hid_t fp, const char *fullpath, unsigned *buf)
{
  herr_t status;

  // open dataset for reading
#ifdef H5_USE_16_API
  hid_t data_id = H5Dopen(fp, fullpath);
#else
  hid_t data_id = H5Dopen(fp, fullpath, H5P_DEFAULT);
#endif

  // read the data
  status =
    H5Dread(data_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

  // close dataset
  H5Dclose(data_id);
  return status;
}

herr_t
GH5_read_grid_data(hid_t fp, const char *path, BucketStruct * buf)
{
  hid_t data_id, type_id;
  herr_t status;

  type_id = GH5_bucketstruct();
#ifdef H5_USE_16_API
  data_id = H5Dopen(fp, path);
#else
  data_id = H5Dopen(fp, path, H5P_DEFAULT);
#endif
  status = H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

  H5Tclose(type_id);
  H5Dclose(data_id);

  return status;
}

herr_t
GH5_WriteS(hid_t fp, const char *path, int dims[], void *buf,
           int unused1, int unused2, const int DATAYPE)
{
  hid_t data_id;
  hid_t space_id;
  herr_t status;
  hsize_t size[2] = { dims[0], dims[1] };

  //create data spaces
  if ( dims[1] > 1)
    space_id = H5Screate_simple(2, size, 0);
  else
    space_id = H5Screate_simple(1, size, 0);

  //create dataset for the buffer
  hid_t data_type;

  switch (DATAYPE)
  {
  case INTTYPE:
    data_type = H5Tcopy(H5T_NATIVE_INT);
    break;
  case UINTTYPE:
    data_type = H5Tcopy(H5T_NATIVE_UINT);
    break;
  case FLOATTYPE:
    data_type = H5Tcopy(H5T_NATIVE_FLOAT);
    break;
  case DOUBLETYPE:
    data_type = H5Tcopy(H5T_NATIVE_DOUBLE);
    break;
  case CHARTYPE:
    data_type = H5Tcopy(H5T_NATIVE_CHAR);
    break;
  default:
    fprintf(stderr, "ERROR: unknown data type in %s.\n", __FILE__);
    exit(1);
  }

#ifdef H5_USE_16_API
  data_id = H5Dcreate(fp, path, data_type, space_id, H5P_DEFAULT);
#else
  data_id =
    H5Dcreate(fp, path, data_type, space_id, H5P_DEFAULT, H5P_DEFAULT,
              H5P_DEFAULT);
#endif
  //write data
  status = H5Dwrite(data_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

  //close dataset and dataspcae
  status = H5Dclose(data_id);
  if (status > 0)
  {
    return status;
  }
  status = H5Sclose(space_id);

  return status;
}

#ifdef PARALLEL_IO
herr_t
GH5_WriteP (hid_t fp, const char *path, int dims[], void *buf,
            int start_in, int count_in, const int DATAYPE)
{
  hid_t data_id;
  hid_t space_id;
  hsize_t size[2] = { dims[0], dims[1] };
  hsize_t start[2] = { start_in, 0 };
  hsize_t count[2] = { count_in, 0 };
  herr_t status;

  //create data spaces
  if ( dims[1] > 1 )
    space_id = H5Screate_simple(2, size, NULL);
  else
    space_id = H5Screate_simple(1, size, NULL);

  //create dataset for the buffer
  hid_t data_type;

  switch (DATAYPE)
  {
  case INTTYPE:
    data_type = H5Tcopy(H5T_NATIVE_INT);
    break;
  case UINTTYPE:
    data_type = H5Tcopy(H5T_NATIVE_UINT);
    break;
  case FLOATTYPE:
    data_type = H5Tcopy(H5T_NATIVE_FLOAT);
    break;
  case DOUBLETYPE:
    data_type = H5Tcopy(H5T_NATIVE_DOUBLE);
    break;
  case CHARTYPE:
    data_type = H5Tcopy(H5T_NATIVE_CHAR);
    break;
  default:
    fprintf(stderr, "ERROR: Unknown dataype in %s\n", __FILE__);
    exit(13);
  }
#  ifdef H5_USE_16_API
  data_id = H5Dcreate(fp, path, data_type, space_id, H5P_DEFAULT);
#  else
  data_id =
    H5Dcreate(fp, path, data_type, space_id, H5P_DEFAULT, H5P_DEFAULT,
              H5P_DEFAULT);
#  endif

  // create file dataspace independently
  hid_t file_space = H5Dget_space(data_id);

  assert(file_space != HDF5_FAIL);

  // set offsets for the file
  status =
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
  assert(status != HDF5_FAIL);

  // create property list for write
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

  // create memory dataspace independently
  hid_t mem_space = H5Screate_simple(1, count, NULL);

  // write data, data size is non-zero
  if (count_in > 0)
    status = H5Dwrite(data_id, data_type, mem_space, file_space, plist_id, buf);

  //close dataset and dataspcae
  H5Sclose(mem_space);
  H5Sclose(file_space);

  // close properties id
  H5Pclose(plist_id);

  status = H5Dclose(data_id);
  if (status != 0)
  {
    return status;
  }
  status = H5Sclose(space_id);
  return status;
}

#endif

herr_t
GH5_write_grid_data(hid_t fp, const char *path, int size, BucketStruct * buf)
{
  int rank = 1;
  hid_t structid, dataset;
  hsize_t dims[] = { size };

  // create new dataspace
  hid_t space = H5Screate_simple(rank, dims, NULL);

  // new dataype id
  structid = GH5_bucketstruct();

  // create data type
#ifdef H5_USE_16_API
  dataset = H5Dcreate(fp, path, structid, space, H5P_DEFAULT);
#else
  dataset =
    H5Dcreate(fp, path, structid, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  // write data to file
  herr_t status =
    H5Dwrite(dataset, structid, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

  // close space, type and dataset identifiers
  H5Tclose(structid);
  H5Sclose(space);
  H5Dclose(dataset);

  // return status
  return status;
}

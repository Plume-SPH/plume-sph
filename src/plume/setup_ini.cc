/*
 * setup_ini.cc
 *
 *  Created on: Mar 24, 2015
 *      Author: zhixuanc
 */

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <mpi.h>

#include "particler.h"
#include "constant.h"
#include "parameters.h"
#include "sph_header.h"
#include "IndMap.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

bool whether_most_out(double *mincrd, double *maxcrd,  double *mindom, double *maxdom)
{
	bool del = false;
	int i;

	/*In z direction, the lower part of the domain was not extended to 2.5 bucket-size
	 *So do not check mincrd[2]
	 * */
	for (i=0; i<DIMENSION-1; i++)
		if (mincrd[i] == mindom[i])
		{
			del = true;
			return del;
		}

	for (i=0; i<DIMENSION; i++)
		if (maxcrd[i] == maxdom[i])
		{
			del = true;
			return del;
		}

	return del;
}

// Sort Container by int hash function
bool sortByHash(const IndMap &lhs, const IndMap &rhs) { return lhs.hash < rhs.hash; }

// Sort Container by index function
bool sortByIndex(const IndMap &lhs, const IndMap &rhs) { return lhs.index < rhs.index; }

// hash equal
struct HashIs {
    HashIs( int s ) : toFind(s) { }
    bool operator() (const IndMap &n)
        { return n.get_hash() == toFind; }
    int toFind;
};
// index equal
struct IndexIs {
    IndexIs( int s ) : toFind(s) { }
    bool operator() (const IndMap &n)
        { return n.get_index() == toFind; }
    int toFind;
};
//myid equal
struct MyidIs{
    MyidIs( int s ) : toFind(s) { }
    bool operator() (const Initial_Index_Map &n)
        { return n.get_iidex() == toFind; }
    int toFind;
};

const double sqrt2 = 1.41421356237310;

int
setup_ini(int myid, THashTable * P_table, HashTable * BG_mesh,
                TimeProps * timeprops, SimProps * simprops, int numprocs, int* my_comm)
{
  int i, j;
//  int rid, cid;

  vector < TKey > pneighs;

//  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION], sj[DIMENSION];
//  double vel[DIMENSION] = {0., 0., 0.};  //initial velocity need to be set to zero
//  double prss, erg, dens, mss;
//  double range[6];
//  double sml2;
//  double pressj, rhoj;
//  unsigned keyr[TKEYLENGTH], keyc[TKEYLENGTH];
  Particle *pi=NULL;
//  Bucket   *bi=NULL;

#ifdef DEBUG
   bool do_search = false;
   unsigned keycheck[TKEYLENGTH] = {81794040, 83755327, 0};
   unsigned keytemp[TKEYLENGTH] ;
#endif

  vector <IndMap> idset;
  std::vector<IndMap>::iterator vit;
  vit = idset.begin();
  //create a hash table iterator instance
  THTIterator *itr = new THTIterator(P_table);
  IndMap initializer_IM;
  // go through particle table
  j=0;
  while ((pi = (Particle *) itr->next ()))
   if (!pi->is_guest())//count number of non-guest particles.
   {

#ifdef DEBUG
	   if (do_search)
	   {
		   for (i = 0; i < TKEYLENGTH; i++)
			    keytemp[i] = pi->getKey ().key[i];

		    if (find_particle (keytemp, keycheck))
			    cout << "The particle found!" << endl;
	   }
#endif

//	  for (i = 0; i < DIMENSION; i++)
//		  xi[i] = *(pi->get_coords() + i);
//
//	  sml2 = 0.5*(pi->get_smlen ());
//	  range[0]=xi[0]-sml2;
//	  range[1]=xi[0]+sml2;
//	  range[2]=xi[1]-sml2;
//	  range[3]=xi[1]+sml2;
//	  range[4]=xi[2]-sml2;
//	  range[5]=xi[2]+sml2;
//
//	  air_prop_hydro (xi, range, &erg, &prss, &dens, &mss);
//
//      //put data back into particle:
//      pi->put_density(dens);
//      pi->put_energy(erg);
//      pi->put_pressure(prss);
//      pi->put_velocity(vel);
//      pi->put_mass(mss);
//
//      //the second variable need to be updated.
//      pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P);
      //Replace the old data with an lib function
	  initial_air(pi, simprops);
      j++;
   }//end of go through all particles

//  int nrp=j; //number of non-guest particles
//  int nrp_all; //Total number of non-guest particles
//
//  //reset size
//  idset.resize(nrp);
//  //Go through all particles again to create local InMap
//  j = 0;
//  itr->reset();
//  while ((pi = (Particle *) itr->next ()))
//  {
//      //create local IndMap
//      if (!(pi->is_guest()))
//      {
//         for (i = 0; i < TKEYLENGTH; i++)
//             keyr[i] = pi->getKey().key[i];
//
//         for (i = 0; i < TKEYLENGTH; i++)
//    	     idset[j].put_key(keyr[i], i);
//
//         rid = P_table->hash (keyr);//R u sure this will return numbering of particle.
//         idset[j].put_hash(rid);
//         idset[j].put_pid(myid);
//         j++;
//      }
//  }//end of go through all particles
//
//  //Establish Initial index map
//   int IIdMap [2][numprocs];
//   MPI_Gather(&nrp,1,MPI_INT,&(IIdMap[0]),1,MPI_INT,0,MPI_COMM_WORLD);
//
//   if (myid == 0)
//	   IIdMap[0][myid]=nrp;
//
//   MPI_Bcast( &IIdMap, 2*numprocs, MPI_INT, 0, MPI_COMM_WORLD );
//   IIdMap[1][0]=0;
//   for (j=1; j<numprocs; j++)
// 	  IIdMap[1][j]= IIdMap[1][j-1]+IIdMap[0][j-1];
//
//   //Need put the matrix into a vector for the convenience of search
//   vector<Initial_Index_Map> IIdMap_v;
//   IIdMap_v.resize(numprocs);
//   vector<Initial_Index_Map> :: iterator iidit = IIdMap_v.begin();
//
//   for (i=0; i<numprocs; i++)
//   {
//	   iidit->put_id(i);
//	   iidit->put_iidex(IIdMap[1][i]);
//	   iidit++;
//   }
//
//  //sort the local hash-index map
//  sort(idset.begin(), idset.end(), sortByHash);
//  //Assign index to the vector
//  for (i=0; i<nrp; i++) idset[i].put_index(i);
//
//  IndMap * array_local = new IndMap[nrp];
//  IndMap * array;
//  vit = idset.begin();
//  i=0;//local index start from 0;
//  for (vit=idset.begin(); vit!=idset.end(); vit++)
//  {
//	  array_local[i].put_key(vit->get_key());
//	  array_local[i].put_hash(vit->get_hash());
//	  array_local[i].put_pid(vit->get_pid());
//	  array_local[i].put_index(vit->get_index());
//	  i++;
//  }
//
//  //Create 5th datatype for IndMap
//  create_indmpi_struct();
//
//  //exchange index information: IndMap Among neighbors
//  exchange_indmap ( numprocs, myid, my_comm, array_local, &array, nrp, &nrp_all);
//
//  //put data in array back into a vector container
//  vector <IndMap> idset_all;
//  idset_all.resize(nrp_all);
//  std::vector<IndMap>::iterator vit2 ;
//  vit2 = idset_all.begin() ;
//  for (i=0; i<nrp_all; i++)
//  {
//	  vit2->put_key(array[i].get_key());
//	  vit2->put_hash(array[i].get_hash());
//	  vit2->put_pid(array[i].get_pid());
//	  vit2->put_index(array[i].get_index());
//	  vit2++;
//  }
//
//  //allocate space for A and bitr->reset();
//  double **A = new double*[nrp];
//  for(int i = 0; i <nrp; ++i)
//      A[i] = new double[nrp_all];
//  double *b = new double[nrp];
//
//  //initial value for A and b
//  for (rid = 0; rid <nrp; rid++)
//	  b[rid]=0;
//
//  for (rid = 0; rid <nrp; rid++)
//	  for (cid = 0; cid <nrp; cid++)
//		  A[rid][cid] = 0;
//
//  int hash;
//  int gii, lid; //gii: global initial lid: local index
//  int process_id;
//  itr->reset();
//  double wght;
//  while ((pi = (Particle *) itr->next ()))
//   {
//       //need make sure that the particle is not guest, && the particles are air phase.
//	   if ((!(pi->is_guest())) && pi->which_phase() == 1)
//	   {
//		   process_id=myid;
//		   prss = ( pi-> get_pressure ());
//  	       for (i = 0; i < TKEYLENGTH; i++)
//  	 	       keyr[i] = pi->getKey().key[i];
//
//  	       //The following section is find index from key and myid
//  	       hash=P_table->hash(keyitr->reset();r);
//  	       vit2 = find_if (idset_all.begin(), idset_all.end(), HashIs(hash));
//  	       iidit = find_if (IIdMap_v.begin(), IIdMap_v.end(), MyidIs(process_id));
//  	       gii=iidit->get_id();
//  	       if (vit2 != idset_all.end())
//  	           lid = vit2->get_index();//R u sure this will return numbering of particle.
//           rid = gii + lid;
//
//  	       b[rid] = prss;
//
//	       // expanded smoothing length for Momentum equation
//	        double hi = pi->get_smlen();
//		    double supp = 3 * hi;
////	        wnorm =0;
// 	       // list of neighbors
// 	       vector < TKey > pneighs = pi->get_neighs();
// 	       vector < TKey >::iterator p_itr;
// 	       //go through neighbours
// 	       for (p_itr = pneighs.begin(); p_itr != pneighs.end(); p_itr++)
// 	       {
//
// 	           Particle *pj = (Particle *) P_table->lookup(*p_itr);
// 	           assert (pj);
//
// 	           //need find the processor to which the particle belong to
// 	           process_id = pj->get_my_processor ();
//
// 	           double dist_sq = 0;
// 	           double dist = 0;
//
// 	           for (i = 0; i < DIMENSION; i++)
// 	           {
// 	                dx[i] = xi[i] - *(pj->get_coords() + i);
// 	                si[i] = dx[i] / hi;
// 	           }
//
// 	           // if dx < 3 * sqrt(2) * h
// 	           if (in_support(dx, supp))
// 	           {
// 	    	       for (i = 0; i < TKEYLENGTH; i++)
// 	    	 	   keyc[i] = pj->getKey().key[i];
//
// 	     	       //The following section is to find index from key and myid
// 	     	       hash=P_table->hash(keyc);
// 	     	       vit2 = find_if (idset_all.begin(), idset_all.end(), HashIs(hash));
// 	     	       iidit = find_if (IIdMap_v.begin(), IIdMap_v.end(), MyidIs(process_id));
// 	     	       gii=iidit->get_id();
// 	     	       if (vit2 != idset_all.end())
// 	     	           lid = vit2->get_index();//R u sure this will return numbering of particle.
// 	               cid = gii + lid;
//
// 	                      rhoj = (pj->get_mass());
// 	                      pressj = (pj->get_pressure());
//
// 	               wght = weight(si, hi);
// 	               //wnorm += wght * (pj->get_mass()) / (pj->get_density()); //mj is unknown, what should I do
// 	               A[rid][cid]=pressj*wght/rhoj;
// 	             }
// 	           } // end loop over neighs
//	        }
//         }
//
//   // Solve system of equations
//   double *m = new double[nrp];
//   double mss;itr->reset();
//   solve_mass(A, b, m, nrp, nrp_all, myid);
//
//   //put m back
//   itr->reset();
//   while ((pi = (Particle *) itr->next ()))
//     {
//         //need make sure that the particle is not guest
//  	   if ((!(pi->is_guest())) && pi->which_phase() == 1)
//  	   {
//  		   process_id=myid;
//    	   for (i = 0; i < TKEYLENGTH; i++)
//    	 	    keyr[i] = pi->getKey().key[i];
//
//    	   //The following section is find index from key and myid
//    	    hash=P_table->hash(keyr);
//    	    vit2 = find_if (idset_all.begin(), idset_all.end(), HashIs(hash));
//    	    iidit = find_if (IIdMap_v.begin(), IIdMap_v.end(), MyidIs(process_id));
//    	    gii=iidit->get_id();
//    	    if (vit2 != idset_all.end())
//    	        lid = vit2->get_index();//R u sure this will return numbering of particle.
//             rid = gii + lid;
//
//    	     mss=m[rid];
//    	     pi->put_mass(mss);
//  	        }
//       }
//
//  // iterate over hashtable to update secondary state_variables
//  THTIterator *it2 = new THTIterator(P_table);
//  while ((pi = (Particle *) it2->next()))
//    if (pi->need_neigh())
//       pi->update_state_vars();
//

  /*go through all buckets and delete the most outside layer of bucket
   * What I need do is:
   * 1) inactive the most outside bucket ---> do not do this, do this in a separate process
   * 2) clear up these bucket ---> clear up particles to save memory
   * 3) remove these particles in these inactivated buckets
   *
   * note: for guest bucket, just inactive it and delete all guest particles.
   * as the same action will be taken in guest's home process,
   * this will automatically synchronize all process.
   * */
//  double mincrd[DIMENSION], maxcrd[DIMENSION];
//  double mindom[DIMENSION], maxdom[DIMENSION];
//  bool del;
//  for (i = 0; i < DIMENSION; i++)
//  {
//    mindom[i] = *(P_table->get_minDom() + i);
//    maxdom[i] = *(P_table->get_maxDom() + i);
//  }
//
//  HTIterator * itr2 = new HTIterator (BG_mesh);
//  Bucket * Bnd_buck;
//
//  vector < TKey > plist;
//  vector < TKey >::iterator p_itr;
//  unsigned tempkey[TKEYLENGTH];
//
//  while ((Bnd_buck = (Bucket *) itr2->next ()))
//  {
//	   /*It is OK to clear up guest buckets
//	   *---> in that way, synchronization is not necessary!
//	   *---> it is certain kind of trade off between computing and communication.
//	  */
//	   for (i=0; i<DIMENSION; i++)
//                {
//                    mincrd[i] = *(Bnd_buck->get_mincrd () + i);
//		    maxcrd[i] = *(Bnd_buck->get_maxcrd () + i);
//                }
//
//		del = whether_most_out(mincrd, maxcrd, mindom, maxdom);
//
//		if (del)
//		{
//			plist = Bnd_buck->get_particle_list ();
//			for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
//			{
//	    		    pi = (Particle *) P_table->lookup(*p_itr);
//	    		    assert (pi);
//
//	    		for (i = 0; i<TKEYLENGTH; i++)
//	    		   tempkey[i] = pi->getKey().key[i];
//
//	    	           P_table->remove(tempkey);
//			}
//
//			Bnd_buck->empty_plist ();
////			Bnd_buck->mark_inactive ();
//		}
//  }//end of go through all buckets

//  // clean up
//  delete itr2;
  delete itr;
//  delete it2;
//  // need clear up more variables.
//  delete[] array_local;
//  delete[] array;
//  delete[] A;
//  delete[] b;
//  delete[] m;
//  vit = idset.erase (idset.end());
//  vit2 = idset_all.erase (idset_all.end());
//  iidit = IIdMap_v.erase (IIdMap_v.end());;
//  delete [] array_local;
//  delete [] array;

  return 0;
}

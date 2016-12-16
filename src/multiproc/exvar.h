
/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: exvar.h 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifndef EXVAR_H
#  define EXVAR_H

#  include <mpi.h>
extern MPI_Datatype BRIEF_BUCKET_TYPE;
extern MPI_Datatype BUCKET_ADD_TYPEE;
extern MPI_Datatype BUCKET_TYPE;
extern MPI_Datatype PARTICLE_TYPE;
extern MPI_Datatype BND_IMAGE_TYPE;
extern MPI_Datatype LB_VERT_TYPE;
extern MPI_Datatype INVOLVED_HEADER_TYPE;

#endif // EXVAR__H

#
# automake template for src/repartition directory
#
# $Id: Makefile.am,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $
#

noinst_LIBRARIES = libmultiproc.a

libmultiproc_a_SOURCES = \
   BSFC_create_bins.cc \
   BSFC_refine_partition.cc \
   BSFC_update_and_send_elements.cc \
   exvar.h \
   move_bnd_img.cc \
   move_data.cc \
   mpi_struct.cc \
   multiproc.h \
   pack_data.h pack_data.cc \
   repartition_BSFC.h repartition_BSFC.cc \
   send_reflections.cc

libmultiproc_a_CPPFLAGS = -I$(top_srcdir)/datstr \
                          -I$(top_srcdir)/src \
                          -I$(top_srcdir)/debug

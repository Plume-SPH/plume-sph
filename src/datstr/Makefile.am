#
# automake template for map directory
#
# $Id:$
#

noinst_LIBRARIES = libbgmesh.a

libbgmesh_a_SOURCES = bgmesh.h \
                      bnd_image.h bnd_image.cc \
                      bucket.h bucket.cc \
                      hashtab.h hashtab.cc \
                      thashtab.h thashtab.cc \
                      hilbert.h hilbert.c \
                      mesh_update.cc \
                      particle.h particle.cc \
                      properties.h \
                      Involved_header.h\
                      meteo.h meteo.cc

#if PARALLEL_BUILD
libbgmesh_a_CPPFLAGS = -I$(top_srcdir)/src/plume \
                       -I$(top_srcdir)/src/multiproc \
                       -I$(top_srcdir)/src/debug
#else
#libbgmesh_a_CPPFLAGS = -I$(top_srcdir)/src/plume \
#                       -I$(top_srcdir)/src/debug
#endif                       

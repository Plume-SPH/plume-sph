#
# automake template for map directory
#
# $Id:$
#
# the last 2 in INCLUDES and LDADD only needed when debug, 

bin_PROGRAMS = preprocess

preprocess_SOURCES = buckstr.h \
                     createfunky.cc \
                     preprocess.h preprocess.cc  \
                     gisapi.h inc_plane.cc

INCLUDES = -I$(top_srcdir)/plume \
           -I$(top_srcdir)/gisapi \
           -I$(top_srcdir)/datstr 
           

preprocess_LDADD = $(top_builddir)/plume/libgmfghdf5.a \
                   $(top_builddir)/gisapi/libgisapi.a \
                   $(top_builddir)/datstr/libbgmesh.a 

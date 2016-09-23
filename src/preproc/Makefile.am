#
# automake template for map directory
#
# $Id:$
#
# the last 2 in INCLUDES and LDADD only needed when debug, 

bin_PROGRAMS = preprocess

preprocess_SOURCES = buckstr.h \
                     createfunky.cc \
                     preprocess.h preprocess.cc

INCLUDES = -I$(top_srcdir)/src/plume  \
           -I$(top_srcdir)/src/datstr 
           

preprocess_LDADD = $(top_builddir)/src/plume/libgmfghdf5.a \
                   $(top_builddir)/src/datstr/libbgmesh.a 

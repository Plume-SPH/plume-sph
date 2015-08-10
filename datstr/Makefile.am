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
                      Involved_header.h

libbgmesh_a_CPPFLAGS = -I$(top_srcdir)/src \
                       -I$(top_srcdir)/multiproc \
                       -I$(top_srcdir)/debug


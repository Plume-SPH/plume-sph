#
# automake template for src/repartition directory
#
# $Id: Makefile.am,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $
#

noinst_LIBRARIES = libdebug.a

libdebug_a_SOURCES = \
            debug_header.h \
            particle_deb.cc \
            debug_bucket_lib.cc \
            debug_lib.cc 

libdebug_a_CPPFLAGS = -I$(top_srcdir)/datstr \
                      -I$(top_srcdir)/src \
                      -I$(top_srcdir)/multiproc \
                      -I$(top_srcdir)/preproc

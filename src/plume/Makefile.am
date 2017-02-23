#
# automake template for map directory
#
# $Id:$
#

noinst_LIBRARIES = libgmfgsph.a libgmfghdf5.a

libgmfgsph_a_SOURCES = sph_header.h sph_lib.cc
libgmfgsph_a_CPPFLAGS = -I$(top_srcdir)/src/datstr

libgmfghdf5_a_SOURCES = hdf5calls.cc hdf5calls.h
libgmfghdf5_a_CPPFLAGS = -I$(top_srcdir)/src/preproc \
                         -I$(top_srcdir)/src/datstr

bin_PROGRAMS = particler

particler_SOURCES = add_air.cc \
                    add_pressure_ghost.cc \
                    add_wall_ghost.cc \
                    adapt_domain.cc \
                    bcond.cc \
                     calculate_gradient.cc \
                    constant.h \
                    dataread.cc \
                    density_massfrac.cc  \
                    h5part_writer.cc \
                    IndMap.h  \
                    matlab_writer.cc \
                    mom_engr_update.cc \
                    neighs.cc \
                    outforms.h \
                    particler.h \
                    parameters.h \
                    scan_outside_layer.cc \
                    setup_erupt.cc \
                    setup_influx.cc \
                    setup_ini.cc\
                    smooth_velocity.cc \
                    shift_brief_buck.cc \
                    timestep.cc \
                    update_out_layer.cc \
                    update_pos.cc \
                    write_output.cc \
                    main.cc

INCLUDES = -I$(top_srcdir)/src/datstr \
           -I$(top_srcdir)/src/preproc \
           -I$(top_srcdir)/src/multiproc\
           -I$(top_srcdir)/src/debug

particler_LDADD = $(top_builddir)/src/datstr/libbgmesh.a \
                  $(top_builddir)/src/multiproc/libmultiproc.a \
                  $(top_builddir)/src/debug/libdebug.a \
                  libgmfghdf5.a libgmfgsph.a



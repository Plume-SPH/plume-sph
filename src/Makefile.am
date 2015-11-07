#
# automake template for map directory
#
# $Id:$
#

noinst_LIBRARIES = libgmfgsph.a libgmfghdf5.a

libgmfgsph_a_SOURCES = sph_header.h sph_lib.cc

libgmfghdf5_a_SOURCES = hdf5calls.cc hdf5calls.h
libgmfghdf5_a_CPPFLAGS = -I$(top_srcdir)/preproc \
                         -I$(top_srcdir)/datstr

bin_PROGRAMS = particler

particler_SOURCES = add_air.cc \
                    add_pressure_ghost.cc \
                    add_wall_ghost.cc \
                    adapt_domain.cc \
                    bcond.cc \
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
                    setup_erupt.cc\
                    setup_ini.cc\
                    smooth_velocity.cc \
                    timestep.cc \
                    update_out_layer.cc \
                    update_pos.cc \
                    write_output.cc \
                    main.cc

INCLUDES = -I$(top_srcdir)/datstr \
           -I$(top_srcdir)/preproc \
           -I$(top_srcdir)/multiproc\
           -I$(top_srcdir)/debug

particler_LDADD = $(top_builddir)/datstr/libbgmesh.a \
                  $(top_builddir)/multiproc/libmultiproc.a \
                  $(top_builddir)/debug/libdebug.a \
                  libgmfghdf5.a libgmfgsph.a

#funky2ascii_SOURCES = funky2txt.cc hdf5calls.h
#funky2ascii_LDADD = libgmfghdf5.a


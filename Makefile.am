#
# automake template for the top-most source directory
#
# $Id: $
#

#EXTRA_DIST =

if PARALLEL_BUILD
    SUBDIRS = ./src/gisapi ./src/datstr ./src/multiproc ./src/debug ./src/plume ./src/preproc 
else
    SUBDIRS = ./src/gisapi ./src/datstr ./src/debug ./src/plume ./src/preproc
endif

if BUILD_DOCS
    OPT_DIRS = doc
endif

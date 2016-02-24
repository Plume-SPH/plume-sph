#
# automake template for the top-most source directory
#
# $Id: $
#

#EXTRA_DIST =

if PARALLEL_BUILD
    SUBDIRS = gisapi datstr multiproc debug plume preproc 
else
    SUBDIRS = gisapi datstr debug plume preproc
endif

if BUILD_DOCS
    OPT_DIRS = doc
endif

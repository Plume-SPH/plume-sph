#
# automake template for src/gisapi directory
#
# $Id: Makefile.am 120 2007-06-07 19:21:02Z dkumar $
#

noinst_LIBRARIES = libgisapi.a

libgisapi_a_SOURCES = GisApi.cc GisAscFile.cc GisBinFile.cc GisCats.cc \
	GisRasterHdr.cc GisColors.cc GisLines.cc GisLabels.cc GisSPRFile.cc \
	GisGrid.cc GisTriFile.cc GisApi.h GisAscFile.h GisBinFile.h GisCats.h GisRasterHdr.h \
	GisColors.h GisLines.h GisLabels.h GisSPRFile.h GisGrid.h GisTriFile.h
        
if HAVE_GDAL
  libgisapi_a_SOURCES += gmfg_GdalApi.h GdalApi.cc
endif

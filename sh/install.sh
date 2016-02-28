# This is the old way to build and install:
setenv PETSC_ARCH arch-linux2-c-debug
setenv PETSC_DIR /rohit1/data/users/zhixuanc/Soft/petsc
setenv PATH /rohit1/data/users/zhixuanc/Soft/petsc/arch-linux2-c-debug/bin:$PATH
setenv HDF5 "/var/lib/yum/yumdb/h/daf1a055f356acdc317ca1d8f4cccf22dfb72714-hdf5-1.8.5.patch1-9.el6-x86_64"

aclocal
autoconf
automake
automake --add-missing

./configure --without-gdal --with-hdf5=$HDF5 CC=mpicc CXX=mpicxx FC=mpif90 CXXFLAGS="-DHS_USE_16_API -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include" LDFLAGS="-DHS_USE_16_API -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc"

#enable debug: configure
./configure --without-gdal --with-hdf5 --enable-debug CC=mpicc CXX=mpicxx FC=mpif90 CXXFLAGS="-DHS_USE_16_API -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include" LDFLAGS="-DHS_USE_16_API -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc"

#enable multiprocess
./configure --without-gdal --with-hdf5 --enable-debug --enable-parallel CC=mpicc CXX=mpicxx FC=mpif90 CXXFLAGS="-DHS_USE_16_API -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include " LDFLAGS="-DHS_USE_16_API -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc"

make clean
make distclean


make
make install

#The new way to build and install:
setenv PATH /rohit1/data/users/zhixuanc/Soft/petsc/arch-linux2-c-debug/bin:$PATH
setenv HDF5 "/var/lib/yum/yumdb/h/daf1a055f356acdc317ca1d8f4cccf22dfb72714-hdf5-1.8.5.patch1-9.el6-x86_64"

make clean
make distclean

aclocal
autoconf
automake --add-missing

./configure --without-gdal --with-hdf5 --enable-debug --enable-parallel CC=mpicc CXX=mpicxx FC=mpifort CXXFLAGS="-std=c++0x" | tee myconfig.out

make | tee mymake.out
make install | tee mymakeinstall.out


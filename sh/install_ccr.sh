module load hdf
aclocal
autoconf
automake --add-missing

#enable multiprocess
./configure --without-gdal --with-hdf5 --enable-debug --enable-parallel CC=mpicc CXX=mpicxx FC=mpif90

make clean
make distclean


make
make install

#load module
module load hdf5/1.8.15p1
#clean up
make clean
make distclean

#generate Makefile.in and config.h.in and other necessary files ---> need Makefile.am and configure.ac
aclocal
autoconf
automake --add-missing


#configure and generate Makefile --->options are flexible
./configure --without-gdal --with-hdf5 --enable-debug --enable-parallel CC=mpicc CXX=mpicxx FC=mpifort CXXFLAGS="-pg" 2>&1 | tee myconfig.out

make 2>&1 | tee mymake.out
make install 2>&1| tee mymakeinstall.out



#The new way to build and install:
setenv PATH /keith1/data/users/zhixuanc/Soft/petsc/arch-linux2-c-debug/bin:$PATH
setenv HDF5 "/var/lib/yum/yumdb/h/daf1a055f356acdc317ca1d8f4cccf22dfb72714-hdf5-1.8.5.patch1-9.el6-x86_64"

make clean
make distclean

aclocal
autoconf
automake --add-missing
#remove --enable-debug to turn off debug mode
#remove --enable-parallel to use serial version : This is not recommended
./configure --without-gdal --with-hdf5 --enable-debug --enable-parallel CC=mpicc CXX=mpicxx FC=mpifort CXXFLAGS="-pg" | tee myconfig.out

make | tee mymake.out
make install | tee mymakeinstall.out



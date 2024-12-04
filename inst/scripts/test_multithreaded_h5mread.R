# The workhorse behind multithreaded h5mread() is C function
# read_data_4_5_parallel() implemented in src/h5mread_index.c
#
# 4 problems at the moment that get in the way:
#
#   1. The HDF5 library provided by Rhdf5lib was not compiled with
#      --enable-threadsafe so it crashes read_data_4_5_parallel()
#
#   2. So for now we need to compile and link h5mread against our own
#      custom installation of HDF5 compiled with --enable-threadsafe.
#      Unfortunately this flag is not compatible with the HDF5 High Level
#      library so we also need to compile with --disable-hl.
#      This is a real bummer because h5mread does need the HL library
#      (in src/h5dimscales.c).
#      So for the sake of testing read_data_4_5_parallel(), we use the
#      following hack to install h5mread:
#      - add #ifdef XXXooXXX and #endif at the top and bottom of 
#        src/h5dimscales.c
#      - comment out all the lines of the /* h5dimscales.c */ section
#        in src/R_init_h5mread.c
#      This breaks all the functions in R/h5writeDimnames.R but we don't
#      need them to test read_data_4_5_parallel() as they are unrelated to
#      h5chunkdim().
#      We also need to remove Rhdf5lib from the LinkingTo field and replace
#      PKG_LIBS=$(RHDF5LIB_LIBS) with the following lines in src/Makevars
#      (assuming our custom installation of HDF was configured with
#      ./configure --prefix=/home/hpages/hdf5 --enable-threadsafe --disable-hl):
#
#        ## Add -fopenmp flag if OpenMP is available.
#        PKG_CFLAGS=-I/home/hpages/hdf5/include $(SHLIB_OPENMP_CFLAGS)
#        PKG_LIBS=-L/home/hpages/hdf5/lib -lhdf5 -lcrypto -lcurl -lsz -laec -lz -ldl -lm $(SHLIB_OPENMP_CFLAGS)
#
#      Finally we need to add export LD_LIBRARY_PATH=/home/hpages/hdf5/lib
#      to our .profile.
#      After doing all the above, we were able to install the h5mread package.
#
#   3. The sanity checks below fail :-(
#
#   4. read_data_4_5_parallel() is significantly slower than
#      read_data_4_5_sequential() at the moment :-( :-(

library(h5mread)
m2 <- matrix(1:125e6, ncol=1250)  # 100,000 x 1250
HDF5Array::writeHDF5Array(m2, "test.h5", "m2", chunkdim=c(50, 50))
stopifnot(identical(h5mread("test.h5", "m2"), m2))

## TEST 1: We load 15 full chunks.
starts <- list(1:250, 1:150)
m <- h5mread("test.h5", "m2", starts=starts)

## This sanity check fails when using read_data_4_5_parallel()! 
stopifnot(identical(m, m2[starts[[1L]], starts[[2L]]]))

## TEST 2: We load a single value from each of the 50k chunks!
chunkdim <- h5chunkdim("test.h5", "m2")
starts <- list(seq(from=1, to=nrow(m2), by=chunkdim[[1L]]),
               seq(from=1, to=ncol(m2), by=chunkdim[[2L]]))
system.time(m <- h5mread("test.h5", "m2", starts=starts))
#    user  system elapsed 
#   0.683   0.024   0.707  <-- read_data_4_5_sequential()
#   2.341   0.629   2.010  <-- read_data_4_5_parallel(), using 20 threads

## This sanity check fails when using read_data_4_5_parallel()!
stopifnot(identical(m, m2[starts[[1L]], starts[[2L]]]))


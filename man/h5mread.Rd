\name{h5mread}

\alias{get_h5mread_returned_type}
\alias{h5mread}

\title{An alternative to \code{rhdf5::h5read}}

\description{
  An efficient and flexible alternative to \code{rhdf5::h5read()}.
}

\usage{
h5mread(filepath, name, starts=NULL, counts=NULL, noreduce=FALSE,
        as.vector=NA, as.integer=FALSE, as.sparse=FALSE,
        method=0L, use.H5Dread_chunk=FALSE)

get_h5mread_returned_type(filepath, name, as.integer=FALSE)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the HDF5 file where the dataset
    to read from is located, or an \link{H5File} object.

    Note that you must create and use an \link{H5File} object if the HDF5
    file to access is stored in an Amazon S3 bucket. See \code{?\link{H5File}}
    for how to do this.

    Also please note that \link{H5File} objects must NOT be used in the
    context of parallel evaluation at the moment.
  }
  \item{name}{
    The name of the dataset in the HDF5 file.
  }
  \item{starts, counts}{
    \code{starts} and \code{counts} are used to specify the \emph{array
    selection}. Each argument can be either \code{NULL} or a list with
    one list element per dimension in the dataset.

    If \code{starts} and \code{counts} are both \code{NULL}, then the entire
    dataset is read.

    If \code{starts} is a list, each list element in it must be a vector of
    valid positive indices along the corresponding dimension in the dataset.
    An empty vector (\code{integer(0)}) is accepted and indicates an empty
    selection along that dimension. A \code{NULL} is accepted and indicates
    a \emph{full} selection along the dimension so has the same meaning
    as a missing subscript when subsetting an array-like object with \code{[}.
    (Note that for \code{[} a \code{NULL} subscript indicates an empty
    selection.)

    Each list element in \code{counts} must be \code{NULL} or a vector
    of non-negative integers of the same length as the corresponding
    list element in \code{starts}. Each value in the vector indicates how
    many positions to select starting from the associated start value.
    A \code{NULL} indicates that a single position is selected for each
    value along the corresponding dimension.

    If \code{counts} is \code{NULL}, then each index in each \code{starts}
    list element indicates a single position selection along the corresponding
    dimension. Note that in this case the \code{starts} argument is
    equivalent to the \code{index} argument of \code{\link[rhdf5]{h5read}}
    and \code{\link[S4Arrays]{extract_array}} (with the caveat that
    \code{\link[rhdf5]{h5read}} doesn't accept empty selections).

    Finally note that when \code{counts} is not \code{NULL} then the
    selection described by \code{starts} and \code{counts} must be
    \emph{strictly ascending} along each dimension.
  }
  \item{as.vector}{
    Should the data be returned in a vector instead of an array?
    By default (i.e. when set to \code{NA}), the data is returned in an
    ordinary array when reading from a multidimensional dataset, and in
    an ordinary vector when reading from a 1D dataset.
    You can override this by setting \code{as.vector} to \code{TRUE}
    or \code{FALSE}.
  }
  \item{as.integer}{
    If set to \code{TRUE} then the data is loaded in an ordinary array
    (or vector) of \code{type()} \code{"integer"}. This will typically
    reduce the memory footprint of the returned array or vector by half
    if the values in the HDF5 dataset are floating point values.

    Note that, when \code{as.integer=TRUE}, the values loaded from the
    HDF5 dataset get coerced to integers at the C level as early as possible
    so this transformation is very efficient.
  }
  \item{as.sparse}{
    By default \code{h5mread()} returns the data in an ordinary
    array or vector. Use \code{as.sparse=TRUE} to return it in a
    \link[SparseArray]{SparseArray} derivative from the \pkg{SparseArray}
    package. This will significantly reduce the memory footprint
    of the returned object if the HDF5 dataset contains mostly zeros.
  }
  \item{noreduce, method, use.H5Dread_chunk}{
    For testing and advanced usage only. Do not use.
  }
}

\value{
  \code{h5mread()} returns an ordinary array or vector if \code{as.sparse}
  is \code{FALSE} (the default), and a \link[SparseArray]{COO_SparseArray}
  object if \code{as.sparse} is \code{TRUE}.

  \code{get_h5mread_returned_type()} returns the type of the array or
  vector that will be returned by \code{h5mread()}.
  Equivalent to (but more efficient than):
  \preformatted{  typeof(h5mread(filepath, name, rep(list(integer(0)), ndim)))
  } where \code{ndim} is the number of dimensions (a.k.a. \emph{rank} in
  HDF5 jargon) of the dataset.
}

\seealso{
  \itemize{
    \item \link{H5File} objects.

    \item \code{\link[rhdf5]{h5read}} in the \pkg{rhdf5} package.

    \item \code{\link[S4Arrays]{extract_array}} in the \pkg{S4Arrays}
          package.

    \item \link[SparseArray]{COO_SparseArray} objects in the \pkg{SparseArray}
          package.

    \item \code{\link{h5mread_from_reshaped}} to read data from a virtually
          reshaped HDF5 dataset.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------
test_h5 <- system.file("extdata", "test.h5", package="h5mread")
h5ls(test_h5)

m1 <- h5mread(test_h5, "m1")  # 12 x 5 integer matrix

m <- h5mread(test_h5, "m1", starts=list(c(8, 12:7), NULL))
m

## Sanity check:
stopifnot(identical(m1[c(8, 12:7), ], m))

m <- h5mread(test_h5, "m1", starts=list(c(8, 12:7), integer(0)))
m

## Sanity check:
stopifnot(identical(m1[c(8, 12:7), NULL], m))

m2 <- h5mread(test_h5, "m2")  # 4000 x 90 double matrix

m2a <- h5mread(test_h5, "m2", starts=list(31, 1), counts=list(10, 8))
m2a

## Sanity check:
stopifnot(identical(m2[31:40, 1:8], m2a))

m2b <- h5mread(test_h5, "m2", starts=list(31, 1), counts=list(10, 8),
               as.integer=TRUE)
m2b

## Sanity check:
storage.mode(m2a) <- "integer"
stopifnot(identical(m2a, m2b))

a3 <- h5mread(test_h5, "a3")  # 180 x 75 x 4 integer array

starts <- list(c(21, 101), NULL, 3:4)
counts <- list(c( 5,  22), NULL, NULL)
a <- h5mread(test_h5, "a3", starts=starts, counts=counts)
dim(a)
a[1:10, 1:12, ]

## Sanity check:
stopifnot(identical(a3[c(21:25, 101:122), , 3:4, drop=FALSE], a))

## ---------------------------------------------------------------------
## RETURNING THE DATA AS A SPARSE ARRAY
## ---------------------------------------------------------------------

starts <- list(c(21:25, 101:122), NULL, 3:4)
coo <- h5mread(test_h5, "a3", starts=starts, as.sparse=TRUE)
coo

class(coo)  # COO_SparseArray object (see ?COO_SparseArray)
dim(coo)

## Sanity check:
stopifnot(is(coo, "COO_SparseArray"), identical(a, as.array(coo)))

## ---------------------------------------------------------------------
## PERFORMANCE
## ---------------------------------------------------------------------
library(ExperimentHub)
hub <- ExperimentHub()

## With the "sparse" TENxBrainData dataset
## ---------------------------------------
fname0 <- hub[["EH1039"]]
h5ls(fname0)  # all datasets are 1D datasets

index <- list(77 * sample(34088679, 5000, replace=TRUE))
## h5mread() is about 4x faster than h5read():
system.time(a <- h5mread::h5mread(fname0, "mm10/data", index))
system.time(b <- h5read(fname0, "mm10/data", index=index))
stopifnot(identical(a, as.vector(b)))

index <- list(sample(1306127, 7500, replace=TRUE))
## h5mread() is about 14x faster than h5read():
system.time(a <- h5mread::h5mread(fname0, "mm10/barcodes", index))
system.time(b <- h5read(fname0, "mm10/barcodes", index=index))
stopifnot(identical(a, as.vector(b)))

## With the "dense" TENxBrainData dataset
## --------------------------------------
fname1 <- hub[["EH1040"]]
h5ls(fname1)  # "counts" is a 2D dataset

set.seed(33)
index <- list(sample(27998, 300), sample(1306127, 450))
## h5mread() is about 2x faster than h5read():
system.time(a <- h5mread::h5mread(fname1, "counts", index))
system.time(b <- h5read(fname1, "counts", index=index))
stopifnot(identical(a, b))

## Alternatively 'as.sparse=TRUE' can be used to reduce memory usage:
system.time(coo <- h5mread::h5mread(fname1, "counts", index, as.sparse=TRUE))
stopifnot(identical(a, as.array(coo)))

## The bigger the selection, the greater the speedup between
## h5read() and h5mread():
\donttest{
  index <- list(sample(27998, 1000), sample(1306127, 1000))
  ## h5mread() about 4x faster than h5read() (12s vs 48s):
  system.time(a <- h5mread::h5mread(fname1, "counts", index))
  system.time(b <- h5read(fname1, "counts", index=index))
  stopifnot(identical(a, b))

  ## With 'as.sparse=TRUE' (about the same speed as with 'as.sparse=FALSE'):
  system.time(coo <- h5mread::h5mread(fname1, "counts", index, as.sparse=TRUE))
  stopifnot(identical(a, as.array(coo)))
}
}
\keyword{utilities}

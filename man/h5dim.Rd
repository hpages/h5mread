\name{h5dim}

\alias{h5exists}
\alias{h5dim}
\alias{h5chunkdim}

\title{Get the dimensions of an HDF5 dataset}

\description{
  Two convenience functions to obtain the dimensions of an HDF5
  dataset as well as the dimensions of its chunks.
}

\usage{
h5dim(filepath, name, as.integer=TRUE)

h5chunkdim(filepath, name, adjust=FALSE)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to an HDF5 file.
  }
  \item{name}{
    The name (as a single string) of a dataset in the HDF5 file.
  }
  \item{as.integer}{
    By default \code{h5dim()} returns the dimensions of the dataset in
    an integer vector and will raise an error if any dimension is greater
    than \code{.Machine$integer.max} (= 2^31 - 1).

    Use \code{as.integer=FALSE} to support datasets with dimensions greater
    than \code{.Machine$integer.max}. In this case the dimensions are
    returned in a numeric vector.
  }
  \item{adjust}{
    By default \code{h5chunkdim()} returns the dimensions of the chunks
    as reported by \code{H5Pget_chunk} from the C HDF5 API. Note that the
    HDF5 specs allow some or all the dimensions of the chunks to be greater
    than the dimensions of the dataset.

    You can use \code{adjust=TRUE} to request \code{h5chunkdim()} to
    return \emph{adjusted chunk dimensions}, that is, dimensions that do
    not exceed the dimensions of the dataset.
    The \emph{adjusted chunk dimensions} are simply obtained by replacing
    those dimensions in the vector of chunk dimensions that are greater
    than the corresponding dataset dimension with the latter.
  }
}

\value{
  An integer (or numeric) vector of length the number of dimensions
  of the HDF5 dataset.
}

\examples{
test_h5 <- system.file("extdata", "test.h5", package="h5mread")
h5ls(test_h5)

h5dim(test_h5, "m2")
h5chunkdim(test_h5, "m2")

h5dim(test_h5, "a3")
h5chunkdim(test_h5, "a3")
}
\keyword{utilities}

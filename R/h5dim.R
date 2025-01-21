### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .H5Fopen()
### .H5Dopen()
### .H5Gopen()
###

### An undocumented feature of rhdf5::H5Fopen(), rhdf5::H5Dopen(), and
### rhdf5::H5Gopen() is that they won't necessarily throw an error when
### they fail to open the file, dataset, or group, but they can actually
### return a FALSE (with a message).
### The three thin wrappers below detect this situation and throw an error.

.H5Fopen <- function(name, flags=h5default("H5F_ACC_RD"), fapl=NULL)
{
    fid <- suppressMessages(rhdf5::H5Fopen(name, flags=flags, fapl=fapl))
    if (!is(fid, "H5IdComponent"))
        stop(wmsg("failed to open HDF5 file '", name, "'"))
    fid
}

.H5Dopen <- function(h5loc, name, dapl=NULL)
{
    did <- suppressMessages(rhdf5::H5Dopen(h5loc, name, dapl=dapl))
    if (!is(did, "H5IdComponent"))
        stop(wmsg("failed to open HDF5 dataset '", name, "'"))
    did
}

.H5Gopen <- function(h5loc, name)
{
    gid <- suppressMessages(rhdf5::H5Gopen(h5loc, name))
    if (!is(gid, "H5IdComponent"))
        stop(wmsg("failed to open HDF5 group '", name, "'"))
    gid
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5exists()
###

h5exists <- function(filepath, name)
{
    fid <- .H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    H5Lexists(fid, name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5isdataset() and h5isgroup()
###

h5isdataset <- function(filepath, name)
{
    fid <- .H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    did <- try(.H5Dopen(fid, name), silent=TRUE)
    ans <- !inherits(did, "try-error")
    if (ans)
        H5Dclose(did)
    ans
}

h5isgroup <- function(filepath, name)
{
    fid <- .H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    gid <- try(.H5Gopen(fid, name), silent=TRUE)
    ans <- !inherits(gid, "try-error")
    if (ans)
        H5Gclose(gid)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5dim() and h5chunkdim()
### 
    
### Return an object of class H5IdComponent representing an H5 dataset ID.
.get_h5dataset <- function(filepath, name)
{
    if (substr(name, 1L, 1L) != "/")
        name <- paste0("/", name)
    group <- gsub("(.*/)[^/]*$", "\\1", name)
    name <- gsub(".*/([^/]*)$", "\\1", name)
    if (is(filepath, "H5File")) {
        fid <- as(filepath, "H5IdComponent")
    } else {
        fid <- .H5Fopen(filepath, flags="H5F_ACC_RDONLY")
        on.exit(H5Fclose(fid))
    }
    gid <- .H5Gopen(fid, group)
    on.exit(H5Gclose(gid), add=TRUE)
    .H5Dopen(gid, name)
}

### NOT exported but used in the HDF5Array package!
dim_as_integer <- function(dim, filepath, name, what="HDF5 dataset")
{   
    if (is.integer(dim))
        return(dim)
    if (any(dim > .Machine$integer.max)) {
        dim_in1string <- paste0(dim, collapse=" x ")
        if (is(filepath, "H5File"))
            filepath <- path(filepath)
        stop(wmsg("Dimensions of ", what, " are too big: ", dim_in1string),
             "\n\n  ",
             wmsg("(This error is about HDF5 dataset '", name, "' ",
                  "from file '", filepath, "'.)"),
             "\n\n  ",
             wmsg("Please note that the HDF5Array package only ",
                  "supports datasets where each dimension is ",
                  "<= '.Machine$integer.max' (= 2**31 - 1)."))
    }
    as.integer(dim)
}

### The TENxMatrixSeed() constructor from the HDF5Array package calls h5dim()
### with 'as.integer=FALSE' in order to get the dimension of a 1D array of
### length >= 2^31.
h5dim <- function(filepath, name, as.integer=TRUE)
{
    did <- .get_h5dataset(filepath, name)
    on.exit(H5Dclose(did), add=TRUE)
    sid <- H5Dget_space(did)
    on.exit(H5Sclose(sid), add=TRUE)
    dim <- H5Sget_simple_extent_dims(sid)$size
    if (as.integer)
        dim <- dim_as_integer(dim, filepath, name)
    dim
}

### Return NULL or an integer vector parallel to 'h5dim(filepath, name)'.
h5chunkdim <- function(filepath, name, adjust=FALSE)
{
    did <- .get_h5dataset(filepath, name)
    on.exit(H5Dclose(did), add=TRUE)
    pid <- H5Dget_create_plist(did)
    on.exit(H5Pclose(pid), add=TRUE)
    if (H5Pget_layout(pid) != "H5D_CHUNKED")
        return(NULL)
    ## We use rev() to invert the order of the dimensions returned by
    ## H5Pget_chunk(). It seems that H5Pget_chunk() should take care of
    ## this though, for consistency with how rhdf5 handles the order of the
    ## dimensions everywhere else (e.g. see ?H5Sget_simple_extent_dims).
    chunkdim <- rev(H5Pget_chunk(pid))
    chunkdim <- dim_as_integer(chunkdim, filepath, name,
                               what="HDF5 dataset chunks")
    if (adjust) {
        dim <- h5dim(filepath, name, as.integer=FALSE)
        ## A sanity check that should never fail.
        stopifnot(length(chunkdim) == length(dim))
        chunkdim <- as.integer(pmin(dim, chunkdim))
    }
    chunkdim
}


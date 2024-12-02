### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


add_prefix_to_basename <- function(name, prefix=".")
{
    stopifnot(isSingleString(name), isSingleString(prefix))
    slash_idx <- which(safeExplode(name) == "/")
    if (length(slash_idx) == 0L) {
        dname <- ""
        bname <- name
    } else {
        last_slash_idx <- max(slash_idx)
        dname <- substr(name, start=1L, stop=last_slash_idx)
        bname <- substr(name, start=last_slash_idx+1L, stop=nchar(name))
    }
    paste0(dname, prefix, bname)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normarg_h5_filepath() and normarg_h5_name()
###

normarg_h5_filepath <- function(path, what1="'filepath'", what2="the dataset")
{
    if (!isSingleString(path))
        stop(wmsg(what1, " must be a single string specifying the path ",
                  "to the HDF5 file where ", what2, " is located"))
    file_path_as_absolute(path)  # return absolute path in canonical form
}

normarg_h5_name <- function(name, what1="'name'",
                                  what2="the name of a dataset",
                                  what3="")
{
    if (!isSingleString(name))
        stop(wmsg(what1, " must be a single string specifying ",
                  what2, " in the HDF5 file", what3))
    if (name == "")
        stop(wmsg(what1, " cannot be the empty string"))
    if (substr(name, start=1L, stop=1L) == "/") {
        name <- sub("^/*", "/", name)  # only keep first leading slash
    } else {
        name <- paste0("/", name)
    }
    name
}


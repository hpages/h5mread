### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


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


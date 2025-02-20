### =========================================================================
### Summarization of an HDF5 dataset
### -------------------------------------------------------------------------


### Not exported yet, needs more work...
h5summarize <- function(filepath, name, index=NULL, as.integer=FALSE,
                        op, na.rm=FALSE, verbose=FALSE)
{
    if (!is(filepath, "H5File")) {
        filepath <- H5File(filepath, .no_rhdf5_h5id=TRUE)
        on.exit(close(filepath))
    }
    name <- normarg_h5_name(name)

    .Call2("C_h5summarize", filepath, name, index, as.integer,
                            op, na.rm, verbose,
                            PACKAGE="h5mread")
}


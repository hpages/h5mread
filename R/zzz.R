.onLoad <- function(libname, pkgname)
{

}

.onUnload <- function(libpath)
{
    library.dynam.unload("h5mread", libpath)
}


.First.lib <- function(lib,pkg)
{
   library.dynam("adsmooth",pkg,lib)
   cat("adsmooth 0.1-1 loaded\n")
}

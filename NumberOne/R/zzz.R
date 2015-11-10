.PkgEnv <- new.env()
.onLoad <- function(libname, pkgname)
{
  dir.create(paste0(tempdir(),'/NumberOne'))
  path = paste0(tempdir(),'/NumberOne')
  assign("tmpdir",path, envir = .PkgEnv)
  cat('\n Welcome in NumberOne package,
      \n you can easily use it with the user interface by typing NumberOneGUI() in the console. \n\n')
}
.onUnload <- function(libpath)
{
  unlink(get("tmpdir",envir=.PkgEnv), force = T, recursive = T)
}

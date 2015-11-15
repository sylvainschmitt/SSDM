.PkgEnv <- new.env()
.onAttach <- function(libname, pkgname){
  packageStartupMessage('Welcome in SSDM package, you can easily use it with the user interface by typing gui() in the console.')
}
.onLoad <- function(libname, pkgname){
  dir.create(paste0(tempdir(),'/NumberOne'))
  path = paste0(tempdir(),'/NumberOne')
  assign("tmpdir",path, envir = .PkgEnv)
}
.onUnload <- function(libpath)
{
  unlink(get("tmpdir",envir=.PkgEnv), force = T, recursive = T)
}

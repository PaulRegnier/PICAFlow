# Automatically check if a new version is available when PICAFlow is loaded

.onLoad = function(libname, pkgname)
{
  checkNewVersion()
}

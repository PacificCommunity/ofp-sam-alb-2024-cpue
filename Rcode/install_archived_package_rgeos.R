url <- "https://cran.r-project.org/src/contrib/Archive/rgeos/rgeos_0.6-4.tar.gz"
pkgFile <- "rgeos_0.4-1.tar.gz"
download.file(url = url, destfile = pkgFile)

# Install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package tarball
unlink(pkgFile)

require(devtools)
require(roxygen2)

#setwd("~/Dropbox/Projects/055_Karakoc_rpackage/")
#create("EDMhelper")

setwd("~/Dropbox/Projects/055_Karakoc_rpackage/EDMhelper/")
document()
check()

build()
load_all()


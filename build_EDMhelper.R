rm(list=ls())

require(devtools)
require(roxygen2)

#setwd("~/Dropbox/Projects/055_Karakoc_rpackage/")
#create("EDMhelper")
#setwd("~/Dropbox/Projects/055_Karakoc_rpackage/EDMhelper/")

setwd("C:/Users/karakoc/Documents/GitHub/r_package_EDMhelper/EDMhelper")
document()

check()

build()
load_all()




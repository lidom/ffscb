## ################################
##
## Make Package Codes
##
## ################################

getwd()

## Remove pkg 
remove.packages("ffscb")

## Create/update documentation and (re-)write NAMESPACE
devtools::document("ffscb")      

## CRAN-check pkg
# devtools::check("ffscb")       # check the p-package

## Install
devtools::install_local("ffscb")
# library(ffscb)
## #################################


# spnbmd <- read.csv(file = "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/R/fregion_pkg/data/spnbmd.csv")
# View(spnbmd)
# save(spnbmd, file = "ffscb/data/spnbmd.RData")

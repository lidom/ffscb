## ################################
##
## Make Package Codes
##
## ################################

## Remove pkg 
remove.packages("ffscb")

## Create/update documentation and (re-)write NAMESPACE
devtools::document("ffscb")      

## CRAN-check pkg
# devtools::check("ffscb")       # check the package

## Install
devtools::install_local("ffscb", force = TRUE)
##
library("ffscb")
help("ffscb")
## #################################


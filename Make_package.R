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
# devtools::check("ffscb")       # check the p-package

## Install
devtools::install_local("ffscb")
##
library("ffscb")
help("ffscb")
## #################################


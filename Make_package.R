## ################################


## Remove pkg 
remove.packages("ffscb")
## Create/update documentation
devtools::document("ffscb_pkg")      
## CRAN-check pkg
# devtools::check("R/fregion_pkg")       # check the p-package
## Install
devtools::install_local("ffscb_pkg")
## #################################

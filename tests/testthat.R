reticulate::py_config()
reticulate::use_python("/usr/bin/python3")

library(testthat)
library(bigsimr)

test_check("bigsimr")

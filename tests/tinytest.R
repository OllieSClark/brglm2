
if ( requireNamespace("tinytest", quietly=TRUE) ){
  tinytest::test_package("brglm2")
  ## Away from CRAN, run all tests again with the trust-region algorithm
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    op <- options(brglm2.algorithm = "trust_region")
    tinytest::test_package("brglm2")
    options(op)
  }
}


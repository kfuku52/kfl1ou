library(testthat)
library(kfl1ou)

run_kfl1ou_tests <- function() {
  test_check("kfl1ou", stop_on_warning = TRUE)
}

if (.Platform$OS.type == "windows") {
  # Forking is unavailable on Windows. Tests that exercise the parallel branch
  # replace the worker implementation, so make that capability explicit for
  # the full test run; the dedicated fallback test locally restores FALSE and
  # asserts the documented warning.
  with_mocked_bindings(
    run_kfl1ou_tests(),
    l1ou_supports_multicore = function() TRUE,
    .package = "kfl1ou"
  )
} else {
  run_kfl1ou_tests()
}

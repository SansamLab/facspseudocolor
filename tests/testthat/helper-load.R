# Characterization tests exercise internal functions until the public API is
# introduced. Load them from the installed namespace rather than sourcing R/.
test_namespace <- as.list(asNamespace("facspseudocolor"), all.names = TRUE)
list2env(test_namespace, envir = environment())

write_test_csv <- function(path, dna, target) {
  dat <- data.frame(dna = dna, target = target)
  names(dat) <- c("DNA", "Target")
  utils::write.csv(dat, path, row.names = FALSE)
}

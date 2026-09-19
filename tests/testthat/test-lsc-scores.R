source(file.path("..", "..", "R", "lsc_scores.R"))

SIG_PATH <- file.path("..", "..", "data", "signatures")

# A minimal series matrix object, so the tests never touch the network.
fake_sm <- function(expr, gsm = colnames(expr), gpl = "GPL999", gse = "GSE1") {
  list(expr = expr, gpl = gpl, gse = gse,
       meta = list(Sample_geo_accession = list(gsm),
                   Sample_title = list(paste0("s", seq_along(gsm)))))
}

test_that("published pLSC6 coefficients are the ones in the file", {
  sig <- load_signature("pLSC6", sig_dir = SIG_PATH)
  # Elsayed et al. Leukemia 2020, doi:10.1038/s41375-019-0604-8:
  # pLSC6 = DNMT3B*0.189 + GPR56*0.054 + CD34*0.0171 + SOCS2*0.141
  #         + SPINK2*0.109 + FAM30A*0.0516
  published <- c(DNMT3B = 0.189, GPR56 = 0.054, CD34 = 0.0171,
                 SOCS2 = 0.141, SPINK2 = 0.109, FAM30A = 0.0516)
  got <- setNames(sig$coefficient, sig$gene)
  expect_equal(sort(names(got)), sort(names(published)))
  expect_equal(got[names(published)], published)
  expect_true(all(sig$verified))
})

test_that("LSC17 is present but still flagged unverified", {
  sig <- load_signature("LSC17", sig_dir = SIG_PATH)
  expect_equal(nrow(sig), 17)
  expect_false(all(sig$verified))
})

test_that("aliases resolve to the primary gene name", {
  sig <- load_signature("pLSC6", sig_dir = SIG_PATH)
  gpr <- sig[sig$gene == "GPR56", ]
  expect_true("ADGRG1" %in% symbols_for(gpr))
  fam <- sig[sig$gene == "FAM30A", ]
  expect_true("KIAA0125" %in% symbols_for(fam))
})

test_that("scoring an unverified signature is refused unless allowed", {
  expr <- matrix(1, nrow = 17, ncol = 2,
                 dimnames = list(paste0("p", 1:17), c("GSM1", "GSM2")))
  sig <- load_signature("LSC17", sig_dir = SIG_PATH)
  map <- data.frame(probe = paste0("p", 1:17), symbol = sig$gene,
                    stringsAsFactors = FALSE)
  sm <- fake_sm(expr)
  expect_error(score_lsc(sm, "LSC17", probe_map = map, sig_dir = SIG_PATH), "unverified")
})

test_that("the score is the published weighted sum", {
  sig <- load_signature("pLSC6", sig_dir = SIG_PATH)
  expr <- matrix(c(rep(1, 6), rep(2, 6)), nrow = 6,
                 dimnames = list(paste0("p", 1:6), c("GSM1", "GSM2")))
  map <- data.frame(probe = paste0("p", 1:6), symbol = sig$gene,
                    stringsAsFactors = FALSE)
  res <- score_lsc(fake_sm(expr), "pLSC6", probe_map = map, sig_dir = SIG_PATH)
  expect_equal(res$Score[1], sum(sig$coefficient))
  expect_equal(res$Score[2], 2 * sum(sig$coefficient))
  expect_equal(res$GenesFound[1], 6L)
  expect_true(res$CoefficientsVerified[1])
})

test_that("a partial signature is refused rather than scored", {
  sig <- load_signature("pLSC6", sig_dir = SIG_PATH)
  expr <- matrix(1, nrow = 3, ncol = 2,
                 dimnames = list(paste0("p", 1:3), c("GSM1", "GSM2")))
  map <- data.frame(probe = paste0("p", 1:3), symbol = sig$gene[1:3],
                    stringsAsFactors = FALSE)
  expect_error(score_lsc(fake_sm(expr), "pLSC6", probe_map = map, sig_dir = SIG_PATH), "Only 3/6")
})

test_that("duplicate probes collapse on the highest mean, and multi-gene probes drop", {
  expr <- matrix(c(1, 5, 9), nrow = 3, ncol = 2,
                 dimnames = list(c("lo", "hi", "amb"), c("GSM1", "GSM2")))
  map <- data.frame(probe  = c("lo", "hi", "amb"),
                    symbol = c("CD34", "CD34", "CD34///FOO"),
                    stringsAsFactors = FALSE)
  out <- collapse_to_genes(expr, map, "CD34")
  expect_equal(rownames(out$mat), "CD34")
  expect_equal(as.numeric(out$mat[1, 1]), 5)   # the "hi" probe, not "lo" or "amb"
  expect_equal(out$method, "max_mean")
})

test_that("the symbol column is chosen by overlap even when names repeat", {
  # GPL23159 really does publish two columns both called SPOT_ID; selecting by
  # name returns the first, which is the wrong one.
  tab <- data.frame(ID = c("a", "b"), SPOT_ID = c("junk", "junk"),
                    SPOT_ID = c("CD34", "SOCS2"),
                    check.names = FALSE, stringsAsFactors = FALSE)
  j <- pick_symbol_column(tab, prefer_symbols = c("CD34", "SOCS2"))
  expect_equal(j, 3L)
  expect_equal(as.character(tab[[j]]), c("CD34", "SOCS2"))
})

test_that("an un-logged matrix is flagged rather than scored silently", {
  sig <- load_signature("pLSC6", sig_dir = SIG_PATH)
  expr <- matrix(5000, nrow = 6, ncol = 2,
                 dimnames = list(paste0("p", 1:6), c("GSM1", "GSM2")))
  map <- data.frame(probe = paste0("p", 1:6), symbol = sig$gene,
                    stringsAsFactors = FALSE)
  expect_warning(res <- score_lsc(fake_sm(expr), "pLSC6", probe_map = map, sig_dir = SIG_PATH),
                 "un-logged")
  expect_equal(res$InputScale[1], "NOT log-scale")
})

test_that("a series matrix round-trips through the reader", {
  fp <- tempfile(fileext = ".txt")
  writeLines(c(
    '!Series_geo_accession\t"GSE1"',
    '!Series_platform_id\t"GPL999"',
    '!Sample_geo_accession\t"GSM1"\t"GSM2"',
    '!Sample_characteristics_ch1\t"cell line: HL60"\t"cell line: HL60"',
    '!Sample_characteristics_ch1\t"drug treatment: vehicle"\t"drug treatment: AZA"',
    "!series_matrix_table_begin",
    '"ID_REF"\t"GSM1"\t"GSM2"',
    '"p1"\t1.5\t2.5',
    "!series_matrix_table_end"
  ), fp)
  sm <- read_series_matrix(fp)
  expect_equal(sm$gse, "GSE1")
  expect_equal(sm$gpl, "GPL999")
  expect_equal(dim(sm$expr), c(1L, 2L))
  expect_equal(as.numeric(sm$expr["p1", ]), c(1.5, 2.5))
  st <- sample_table(sm)
  expect_equal(st$GSM, c("GSM1", "GSM2"))
  expect_match(st$Characteristics[2], "AZA")
})

# Treatment harmonisation for GEO sample characteristics.
# Shared by R/rank_drugs.R and R/scan_corpus.R so the two cannot diverge.

# ---- harmonisation ---------------------------------------------------------

# Keys that carry a treatment in the seven studies seen so far. GEO has no standard
# for this: the same concept appears as "drug treatment", "agent", "treatment" and
# "cell type treatment" across four studies.
TREAT_KEYS <- "^(drug treatment|agent|treatment|cell type treatment|compound|perturbation|drug)$"
LINE_KEYS  <- "^(cell line|cell type|cell_line)$"
# Some studies record the biology that matters in a genotype field instead of a cell
# line, and some put the line only in the sample title. Both are worth recovering: a
# menin-MLL inhibitor is expected to act in KMT2A-rearranged cells and not in wild type,
# and collapsing that distinction hides the result.
GENO_KEYS  <- "^(genotype|mutation|subtype|karyotype)$"
KNOWN_LINES <- c("THP-?1","HL-?60","KG-?1a?","MOLM-?1[346]","MV4-?11","U-?937","NB4",
                 "OCI-?AML-?\\d+","Kasumi-?\\d+","K-?562","SKM-?1","ME-?1","TF-?1",
                 "HEL","SET-?2","NOMO-?1","EOL-?1","HNT-?34","AML-?193","PL-?21")

# A control arm, however the submitter phrased it.
CONTROL_RE <- paste0(
  "vehicle|dmso|untreated|^\\s*none\\s*$|^\\s*control\\b|negative control|",
  "^\\s*0(\\.0+)?\\s*(u|n|m)m\\b|mock"
)

# Dose, duration and descriptive noise that sit around the agent name.
strip_noise <- function(x) {
  x <- gsub("\\([^)]*\\)", " ", x)                                  # (1 uM), (DMSO 0.01%)
  x <- gsub("\\b\\d+(\\.\\d+)?\\s*(days?|hours?|hrs?|h)\\b", " ", x, ignore.case = TRUE)
  x <- gsub("\\b\\d+(\\.\\d+)?\\s*(u|n|m|p)m\\b", " ", x, ignore.case = TRUE)
  x <- gsub("\\b\\d+(\\.\\d+)?\\s*%", " ", x)
  x <- gsub("\\b(treated|treatment|with|for|of|the|a|an|prodrug|inhibitor)\\b", " ",
            x, ignore.case = TRUE)
  trimws(gsub("\\s+", " ", x))
}

# Returns one row per sample: the raw label, whether it is a control, and a best-effort
# agent name. Agent extraction from free text is a heuristic, so confidence is recorded
# rather than assumed.
harmonise <- function(chars) {
  pairs <- lapply(chars, function(s) {
    bits <- trimws(strsplit(s, "\\|")[[1]])
    kv <- do.call(rbind, lapply(bits, function(b) {
      if (!grepl(":", b)) return(NULL)
      p <- strsplit(b, ":", fixed = TRUE)[[1]]
      c(tolower(trimws(p[1])), trimws(paste(p[-1], collapse = ":")))
    }))
    if (is.null(kv)) matrix(character(0), ncol = 2) else kv
  })

  pick <- function(kv, pattern) {
    if (!nrow(kv)) return(NA_character_)
    hit <- which(grepl(pattern, kv[, 1]))
    if (!length(hit)) NA_character_ else kv[hit[1], 2]
  }

  raw  <- vapply(pairs, pick, character(1), TREAT_KEYS)
  line <- vapply(pairs, pick, character(1), LINE_KEYS)
  geno <- vapply(pairs, pick, character(1), GENO_KEYS)

  is_ctrl <- grepl(CONTROL_RE, raw, ignore.case = TRUE)
  agent   <- ifelse(is_ctrl, "control", strip_noise(raw))
  agent[!is.na(agent) & agent == ""] <- NA_character_

  # A name that is still a sentence is not a drug name.
  conf <- ifelse(is.na(agent), "none",
          ifelse(agent == "control", "control",
          ifelse(nchar(agent) <= 24 & lengths(strsplit(agent, " ")) <= 3,
                 "high", "low")))

  data.frame(TreatmentRaw = raw, CellLine = line, Genotype = geno, IsControl = is_ctrl,
             Agent = agent, AgentConfidence = conf, stringsAsFactors = FALSE)
}


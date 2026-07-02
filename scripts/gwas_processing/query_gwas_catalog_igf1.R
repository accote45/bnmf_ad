#!/usr/bin/env Rscript
# Query GWAS Catalog for IGF-1 summary statistics
# Filters: prioritize key authors, specific ancestries, N >= 5000,
#          exclude burden/CNV/interaction/bivariate GWAS

library(httr)
library(jsonlite)
library(data.table)
library(glue)
library(curl)

base_url <- "https://www.ebi.ac.uk/gwas/rest/api"

# --- Helper: paginate through GWAS Catalog API ---
gwas_paginate <- function(url, max_pages = 50, page_size = 100) {
  all_results <- list()
  page <- 0

  repeat {
    resp <- GET(url,
                add_headers("Accept" = "application/json"),
                query = list(page = page, size = page_size))

    if (status_code(resp) != 200) {
      message(glue("API error: {status_code(resp)} on page {page}"))
      break
    }

    data <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                     simplifyVector = FALSE)
    items <- data[["_embedded"]]

    if (is.null(items) || length(items[[1]]) == 0) break

    all_results <- c(all_results, items[[1]])
    page_info <- data[["page"]]

    total_pages <- page_info[["totalPages"]]
    total_elements <- page_info[["totalElements"]]

    if (page == 0) {
      message(glue("Total elements: {total_elements}, pages: {total_pages}"))
    }

    if (page >= (total_pages - 1) || page >= max_pages) break
    page <- page + 1
    Sys.sleep(0.3)
  }

  return(all_results)
}

# --- Step 1: Search for IGF-1 related studies ---
message("=== Searching GWAS Catalog for IGF-1 studies ===\n")

# Search by disease trait text
search_terms <- c("IGF-1", "IGF1", "insulin-like growth factor 1",
                   "insulin-like growth factor I", "insulin like growth factor")

all_studies <- list()

for (term in search_terms) {
  message(glue("Searching for: '{term}'"))

  url <- glue("{base_url}/studies/search/findByDiseaseTrait")
  resp <- GET(url,
              add_headers("Accept" = "application/json"),
              query = list(diseaseTrait = term, size = 200))

  if (status_code(resp) == 200) {
    data <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                     simplifyVector = FALSE)
    studies <- data[["_embedded"]][["studies"]]
    if (!is.null(studies)) {
      message(glue("  Found {length(studies)} studies"))
      all_studies <- c(all_studies, studies)
    } else {
      message("  Found 0 studies")
    }
  }
  Sys.sleep(0.3)
}

# Also search via EFO trait search
message("\nSearching EFO traits for IGF-1...")
efo_url <- glue("{base_url}/efoTraits/search/findByTraitIgnoreCase")
for (term in c("IGF-1", "IGF1", "insulin-like growth factor")) {
  resp <- GET(efo_url,
              add_headers("Accept" = "application/json"),
              query = list(trait = term, size = 50))

  if (status_code(resp) == 200) {
    data <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                     simplifyVector = FALSE)
    traits <- data[["_embedded"]][["efoTraits"]]
    if (!is.null(traits)) {
      for (tr in traits) {
        message(glue("  EFO trait found: {tr$trait} ({tr$shortForm})"))

        # Get studies for this EFO trait
        trait_studies_url <- glue("{base_url}/efoTraits/{tr$shortForm}/associations")
        # Get the studies link instead
        study_url <- tr[["_links"]][["studies"]][["href"]]
        if (!is.null(study_url)) {
          resp2 <- GET(study_url,
                       add_headers("Accept" = "application/json"),
                       query = list(size = 200))
          if (status_code(resp2) == 200) {
            data2 <- fromJSON(content(resp2, as = "text", encoding = "UTF-8"),
                              simplifyVector = FALSE)
            efo_studies <- data2[["_embedded"]][["studies"]]
            if (!is.null(efo_studies)) {
              message(glue("    {length(efo_studies)} studies for this EFO trait"))
              all_studies <- c(all_studies, efo_studies)
            }
          }
        }
      }
    }
  }
  Sys.sleep(0.3)
}

message(glue("\nTotal raw studies collected: {length(all_studies)}"))

# --- Step 2: Parse study metadata ---
message("\n=== Parsing study metadata ===\n")

parse_study <- function(s) {
  accession <- s[["accessionId"]] %||% NA_character_

  # Extract trait names
  trait_names <- tryCatch({
    traits <- s[["_embedded"]][["efoTraits"]]
    if (!is.null(traits)) {
      paste(sapply(traits, function(t) t[["trait"]] %||% ""), collapse = "; ")
    } else {
      s[["diseaseTrait"]][["trait"]] %||% NA_character_
    }
  }, error = function(e) NA_character_)

  disease_trait <- tryCatch({
    s[["diseaseTrait"]][["trait"]] %||% NA_character_
  }, error = function(e) NA_character_)

  # Author and publication
  author <- s[["publicationInfo"]][["author"]][["fullname"]] %||%
    s[["author"]] %||% NA_character_
  pubdate <- s[["publicationInfo"]][["publicationDate"]] %||%
    s[["publicationDate"]] %||% NA_character_
  pubmedid <- s[["publicationInfo"]][["pubmedId"]] %||% NA_character_
  title <- s[["publicationInfo"]][["title"]] %||% NA_character_

  # Sample size - parse from initialSampleSize
  initial_sample <- s[["initialSampleSize"]] %||% NA_character_

  # Ancestry info from ancestries link
  ancestry_info <- NA_character_

  # Check if full summary stats available
  has_sumstats <- s[["fullPvalueSet"]] %||% FALSE

  data.table(
    accession = accession,
    author = author,
    pubdate = pubdate,
    pubmedid = pubmedid,
    disease_trait = disease_trait,
    efo_trait = trait_names,
    initial_sample = initial_sample,
    has_sumstats = has_sumstats,
    title = title
  )
}

studies_dt <- rbindlist(lapply(all_studies, parse_study), fill = TRUE)

# Deduplicate by accession
studies_dt <- unique(studies_dt, by = "accession")
message(glue("Unique studies: {nrow(studies_dt)}"))

# --- Step 3: Filter studies ---
message("\n=== Applying filters ===\n")

# Show all unique traits found
message("Unique disease traits found:")
for (tr in unique(studies_dt$disease_trait)) {
  message(glue("  - {tr}"))
}

# Extract year from pubdate
studies_dt[, year := as.integer(substr(pubdate, 1, 4))]

# Extract first author last name
studies_dt[, first_author := sub(" .*", "", author)]

# Parse sample sizes from initialSampleSize string
# This field contains text like "5,000 European ancestry individuals"
parse_sample_size <- function(txt) {
  if (is.na(txt)) return(NA_integer_)
  # Extract all numbers, remove commas
  nums <- regmatches(txt, gregexpr("[0-9,]+", txt))[[1]]
  nums <- as.numeric(gsub(",", "", nums))
  if (length(nums) == 0) return(NA_integer_)
  return(sum(nums, na.rm = TRUE))  # total across groups
}

studies_dt[, total_n := sapply(initial_sample, parse_sample_size)]

# Filter: N >= 5000
message(glue("Before N filter: {nrow(studies_dt)} studies"))
studies_dt_filtered <- studies_dt[is.na(total_n) | total_n >= 5000]
message(glue("After N >= 5000 filter: {nrow(studies_dt_filtered)} studies"))

# Filter: exclude burden, CNV, tandem repeat, interaction, bivariate, multivariate
exclude_patterns <- paste0(
  "(?i)(burden|gene-based|cnv|copy number|tandem repeat|",
  "str expansion|interaction|GxE|gene.?environ|",
  "bivariate|multivariate|multi-trait|adjusted for|",
  "conditional|conditional on|mtag)"
)

# Apply exclusion to disease_trait and title
studies_dt_filtered[, exclude := grepl(exclude_patterns, disease_trait, perl = TRUE) |
                      grepl(exclude_patterns, title, perl = TRUE)]
message(glue("Excluding {sum(studies_dt_filtered$exclude, na.rm=TRUE)} studies (burden/CNV/interaction/etc)"))
studies_dt_filtered <- studies_dt_filtered[exclude == FALSE | is.na(exclude)]

# Filter: only keep studies with full summary statistics
message(glue("Studies with summary statistics: {sum(studies_dt_filtered$has_sumstats, na.rm=TRUE)}"))
studies_with_ss <- studies_dt_filtered[has_sumstats == TRUE]

# --- Step 4: Build FTP links and get ancestry details ---
message("\n=== Building FTP download links ===\n")

# GWAS Catalog FTP structure for summary statistics:
# ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCSTXXXXXX-GCSTXXXXXX/GCSTXXXXXX/
# Or via HTTPS:
# https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/

build_ftp_link <- function(accession) {
  # Determine the directory range (grouped by 1000s)
  num <- as.integer(gsub("GCST", "", accession))
  range_start <- (floor((num - 1) / 1000) * 1000) + 1
  range_end <- range_start + 999
  range_dir <- glue("GCST{sprintf('%06d', range_start)}-GCST{sprintf('%06d', range_end)}")

  glue("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{range_dir}/{accession}/")
}

studies_with_ss[, ftp_link := sapply(accession, build_ftp_link)]

# --- Step 5: Get ancestry info for each study ---
message("Fetching ancestry details for studies with summary stats...")

get_ancestry <- function(accession) {
  url <- glue("{base_url}/studies/{accession}/ancestries")
  resp <- GET(url, add_headers("Accept" = "application/json"))

  if (status_code(resp) != 200) return(NA_character_)

  data <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                   simplifyVector = FALSE)
  ancestries <- data[["_embedded"]][["ancestries"]]

  if (is.null(ancestries)) return(NA_character_)

  # Extract ancestry categories and sample sizes for discovery stage
  ancestry_info <- sapply(ancestries, function(a) {
    stage <- a[["type"]] %||% ""
    n <- a[["numberOfIndividuals"]] %||% NA

    # Get ancestry group
    groups <- a[["ancestralGroups"]]
    if (!is.null(groups) && length(groups) > 0) {
      group_names <- sapply(groups, function(g) g[["ancestralGroup"]] %||% "")
      group_str <- paste(group_names, collapse = "/")
    } else {
      group_str <- "Unknown"
    }

    # Get country of recruitment
    country <- tryCatch({
      countries <- a[["countryOfRecruitment"]]
      if (!is.null(countries) && length(countries) > 0) {
        paste(sapply(countries, function(c) c[["countryName"]] %||% ""), collapse = ",")
      } else { "" }
    }, error = function(e) "")

    glue("{stage}: {group_str} (N={n})")
  })

  paste(ancestry_info, collapse = " | ")
}

ancestry_results <- character(nrow(studies_with_ss))
for (i in seq_len(nrow(studies_with_ss))) {
  ancestry_results[i] <- get_ancestry(studies_with_ss$accession[i])
  Sys.sleep(0.3)
}
studies_with_ss[, ancestry := ancestry_results]

# --- Step 6: Priority scoring ---
message("\n=== Scoring and ranking ===\n")

# Priority authors
priority_authors <- c("Karczewski", "Verma", "Chen", "Sakaue", "Gurdasani")

studies_with_ss[, is_priority_author := any(sapply(priority_authors, function(a) {
  grepl(a, first_author, ignore.case = TRUE) | grepl(a, author, ignore.case = TRUE)
})), by = accession]

# Priority ancestry keywords
priority_ancestry <- c("European", "African", "East Asian", "South Asian",
                        "Hispanic", "Latino", "Latin American",
                        "multi-ancestry", "trans-ethnic", "Multi-ancestry")

studies_with_ss[, has_priority_ancestry := any(sapply(priority_ancestry, function(a) {
  grepl(a, ancestry, ignore.case = TRUE) | grepl(a, initial_sample, ignore.case = TRUE)
})), by = accession]

# Score
studies_with_ss[, priority_score := 0L]
studies_with_ss[is_priority_author == TRUE, priority_score := priority_score + 10L]
studies_with_ss[has_priority_ancestry == TRUE, priority_score := priority_score + 5L]
studies_with_ss[, priority_score := priority_score + pmin(as.integer(total_n / 10000), 10L)]

# Sort by priority score descending
setorder(studies_with_ss, -priority_score, -total_n)

# --- Step 7: Output results ---
message("\n=== RESULTS: IGF-1 GWAS with Summary Statistics ===\n")

# Select output columns
output <- studies_with_ss[, .(
  Accession = accession,
  First_Author = first_author,
  Year = year,
  Disease_Trait = disease_trait,
  Ancestry = ancestry,
  Discovery_N = total_n,
  Priority_Author = is_priority_author,
  FTP_Link = ftp_link
)]

# Print nicely
message(glue("Found {nrow(output)} IGF-1 GWAS with downloadable summary statistics\n"))

cat("\n--- Priority Author Studies ---\n\n")
priority <- output[Priority_Author == TRUE]
if (nrow(priority) > 0) {
  for (i in seq_len(nrow(priority))) {
    cat(glue("  [{i}] {priority$Accession[i]}"), "\n")
    cat(glue("      Author: {priority$First_Author[i]} ({priority$Year[i]})"), "\n")
    cat(glue("      Trait:  {priority$Disease_Trait[i]}"), "\n")
    cat(glue("      N:      {format(priority$Discovery_N[i], big.mark=',')}"), "\n")
    cat(glue("      Ancestry: {priority$Ancestry[i]}"), "\n")
    cat(glue("      FTP:    {priority$FTP_Link[i]}"), "\n\n")
  }
} else {
  cat("  None found.\n\n")
}

cat("\n--- Other Studies ---\n\n")
other <- output[Priority_Author == FALSE]
if (nrow(other) > 0) {
  for (i in seq_len(nrow(other))) {
    cat(glue("  [{i}] {other$Accession[i]}"), "\n")
    cat(glue("      Author: {other$First_Author[i]} ({other$Year[i]})"), "\n")
    cat(glue("      Trait:  {other$Disease_Trait[i]}"), "\n")
    cat(glue("      N:      {format(other$Discovery_N[i], big.mark=',')}"), "\n")
    cat(glue("      Ancestry: {other$Ancestry[i]}"), "\n")
    cat(glue("      FTP:    {other$FTP_Link[i]}"), "\n\n")
  }
} else {
  cat("  None found.\n\n")
}

# Also save to file
outfile <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/sumstats/igf1_gwas_catalog_query.tsv"
fwrite(output, outfile, sep = "\t")
message(glue("\nResults saved to: {outfile}"))

# Print as a compact table too
message("\n=== Compact Table ===")
print(output[, .(Accession, First_Author, Year, Disease_Trait, Discovery_N, Priority_Author)])

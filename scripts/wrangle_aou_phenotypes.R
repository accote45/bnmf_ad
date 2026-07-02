suppressPackageStartupMessages({
  library(tidyverse)
  library(bigrquery)
  library(data.table)
  library(lubridate)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

out_dir <- file.path(Sys.getenv("HOME"), "workspaces/duplicateofstatinresponse/multiancestry_polygenic", "phenotypes")
data_dir <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--out-dir") {
    out_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--data-dir") {
    data_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--help") {
    cat("Usage: Rscript wrangle_aou_phenotypes.R [--out-dir /path/to/output] [--data-dir /path/to/cached/exports]\n")
    quit(status = 0)
  } else {
    cat(sprintf("Unknown argument: %s\n", args[i])); quit(status = 1)
  }
}

if (!is.null(data_dir)) {
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("AoU Phenotype Wrangling: T2D & CAD\n")
cat("============================================\n\n")
cat(sprintf("Output directory: %s\n", out_dir))
cat(sprintf("Data directory:   %s\n", ifelse(is.null(data_dir), "(none — all queries via BigQuery)", data_dir)))
cat(sprintf("CDR: %s\n", Sys.getenv("WORKSPACE_CDR")))
cat(sprintf("Billing project: %s\n\n", Sys.getenv("GOOGLE_PROJECT")))

# ---------------------------------------------------------------------------
# Helper: Export BigQuery results and read as data.frame
# ---------------------------------------------------------------------------
read_bq_export_from_workspace_bucket <- function(export_path) {
  col_types <- cols(
    gender = col_character(),
    race = col_character(),
    ethnicity = col_character(),
    sex_at_birth = col_character(),
    standard_concept_code = col_character()
  )
  bind_rows(
    map(
      system2("gsutil", args = c("ls", export_path),
              stdout = TRUE, stderr = TRUE),
      function(csv) {
        message(str_glue("Loading {csv}."))
        chunk <- read_csv(
          pipe(str_glue("gsutil cat {csv}")),
          col_types = col_types,
          show_col_types = FALSE
        )
        chunk
      }
    )
  )
}

run_bq_export <- function(sql, export_name, cache_dir = NULL) {
  cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(export_name, ".csv")) else NULL

  # 1) Local cache: fastest path
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cat(sprintf("  [local cache] Reading %s from %s\n", export_name, cache_file))
    df <- read_csv(cache_file, show_col_types = FALSE)
    cat(sprintf("  Loaded %d rows x %d cols\n", nrow(df), ncol(df)))
    return(df)
  }

  export_path <- file.path(
    Sys.getenv("WORKSPACE_BUCKET"),
    "bq_exports",
    Sys.getenv("OWNER_EMAIL"),
    export_name,
    paste0(export_name, "_*.csv")
  )

  # 2) GCS shards: check if a previous export already exists in the bucket
  gcs_files <- tryCatch(
    system2("gsutil", args = c("ls", export_path), stdout = TRUE, stderr = TRUE),
    warning = function(w) character(0)
  )
  gcs_exists <- length(gcs_files) > 0 && !any(grepl("CommandException", gcs_files))

  if (gcs_exists) {
    cat(sprintf("  [GCS cache] Found existing shards for %s, skipping BigQuery\n", export_name))
  } else {
    # 3) No cache anywhere — run the BigQuery export
    cat(sprintf("  Exporting %s via BigQuery...\n", export_name))
    bq_table_save(
      bq_dataset_query(
        Sys.getenv("WORKSPACE_CDR"),
        sql,
        billing = Sys.getenv("GOOGLE_PROJECT")
      ),
      export_path,
      destination_format = "CSV"
    )
  }

  cat(sprintf("  Reading exported data from GCS...\n"))
  df <- read_bq_export_from_workspace_bucket(export_path)
  cat(sprintf("  Loaded %d rows x %d cols\n", nrow(df), ncol(df)))

  # Save locally for next run
  if (!is.null(cache_file)) {
    write_csv(df, cache_file)
    cat(sprintf("  [local cache] Saved to %s\n", cache_file))
  }

  df
}

# ============================================================================
# 1) DEMOGRAPHICS
# ============================================================================
cat("\n--- Querying demographics ---\n")

demographic_sql <- paste("
    SELECT
        person.person_id,
        person.gender_concept_id,
        p_gender_concept.concept_name as gender,
        person.birth_datetime as date_of_birth,
        person.race_concept_id,
        p_race_concept.concept_name as race,
        person.ethnicity_concept_id,
        p_ethnicity_concept.concept_name as ethnicity,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth 
    FROM
        `person` person 
    LEFT JOIN
        `concept` p_gender_concept 
            ON person.gender_concept_id = p_gender_concept.concept_id 
    LEFT JOIN
        `concept` p_race_concept 
            ON person.race_concept_id = p_race_concept.concept_id 
    LEFT JOIN
        `concept` p_ethnicity_concept 
            ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id 
    LEFT JOIN
        `concept` p_sex_at_birth_concept 
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id", sep="")

demo_df <- run_bq_export(demographic_sql, "demographics", cache_dir = data_dir)

demo_df <- demo_df %>%
  filter(sex_at_birth %in% c("Male", "Female"))

cat(sprintf("  Demographics after filtering: %d individuals\n", nrow(demo_df)))

# ============================================================================
# 2) SELF-REPORT SURVEY (shared across T2D and CAD)
# ============================================================================
cat("\n--- Querying self-report survey (concept 1740639) ---\n")

survey_sql <- paste("
    SELECT
        answer.person_id,
        answer.survey_datetime,
        answer.answer  
    FROM
        `ds_survey` answer   
    WHERE
        (
            question_concept_id IN (SELECT
                DISTINCT concept_id                         
            FROM
                `cb_criteria` c                         
            JOIN
                (SELECT
                    CAST(cr.id as string) AS id                               
                FROM
                    `cb_criteria` cr                               
                WHERE
                    concept_id IN (1740639)                               
                    AND domain_id = 'SURVEY') a 
                    ON (c.path like CONCAT('%', a.id, '.%'))                         
            WHERE
                domain_id = 'SURVEY'                         
                AND type = 'PPI'                         
                AND subtype = 'QUESTION')
        )", sep="")

survey_df <- run_bq_export(survey_sql, "survey_conditions", cache_dir = data_dir)

cat(sprintf("  Survey responses: %d\n", nrow(survey_df)))

# Derive enrollment date per participant (max survey_datetime)
enrollment_df <- survey_df %>%
  mutate(survey_datetime = ymd_hms(survey_datetime)) %>%
  group_by(person_id) %>%
  summarise(enrollment_date = max(survey_datetime, na.rm = TRUE), .groups = "drop")

cat(sprintf("  Enrollment dates derived for %d participants\n", nrow(enrollment_df)))

# ============================================================================
# 3) T2D STATUS
# ============================================================================
cat("\n--- Querying T2D conditions (OMOP concept 201820) ---\n")

diabetes_sql <- paste("
    SELECT
        c_occurrence.person_id,
        c_standard_concept.concept_name as standard_concept_name,
        c_source_concept.concept_name as source_concept_name,
        c_source_concept.concept_code as source_concept_code,
        c_occurrence.condition_start_datetime
    FROM
        ( SELECT
            *
        FROM
            `condition_occurrence` c_occurrence
        WHERE
            (
                condition_concept_id IN (SELECT
                    DISTINCT c.concept_id
                FROM
                    `cb_criteria` c
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id
                    FROM
                        `cb_criteria` cr
                    WHERE
                        concept_id IN (201820)
                        AND full_text LIKE '%_rank1]%'      ) a
                        ON (c.path LIKE CONCAT('%.', a.id, '.%')
                        OR c.path LIKE CONCAT('%.', a.id)
                        OR c.path LIKE CONCAT(a.id, '.%')
                        OR c.path = a.id)
                WHERE
                    is_standard = 1
                    AND is_selectable = 1)
            )) c_occurrence
    LEFT JOIN
        `concept` c_standard_concept
            ON c_occurrence.condition_concept_id = c_standard_concept.concept_id
    LEFT JOIN
        `concept` c_source_concept
            ON c_occurrence.condition_source_concept_id = c_source_concept.concept_id
", sep = "")

diabetes_df <- run_bq_export(diabetes_sql, "diabetes_conditions", cache_dir = data_dir)

# Combine condition + survey sources
cat("  Combining T2D sources (conditions + survey)...\n")

t2d_all_df <- bind_rows(
  # T2D from condition table
  diabetes_df %>%
    transmute(
      person_id,
      t2d = 1L,
      t2d_date = ymd_hms(condition_start_datetime),
      source = "condition"
    ),

  # T2D from self-report survey
  survey_df %>%
    filter(
      str_detect(
        answer,
        regex("type 1 diabetes|type 2 diabetes|other/unknown diabetes",
              ignore_case = TRUE)
      ),
      str_detect(answer, regex("yes", ignore_case = TRUE))
    ) %>%
    transmute(
      person_id,
      t2d = 1L,
      t2d_date = ymd_hms(survey_datetime),
      source = "survey"
    )
)

# Deduplicate: one row per person, earliest date
t2d_df <- t2d_all_df %>%
  group_by(person_id) %>%
  summarise(
    t2d = 1L,
    t2d_date = min(t2d_date, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  T2D cases: %d unique individuals\n", nrow(t2d_df)))

# --- 4a) CAD conditions ---
cat("\n--- Querying CAD conditions (OMOP concept 317576) ---\n")

cad_sql <- paste("
    SELECT
        c_occurrence.person_id,
        c_standard_concept.concept_name as standard_concept_name,
        c_occurrence.condition_start_datetime
    FROM
        ( SELECT
            *
        FROM
            `condition_occurrence` c_occurrence
        WHERE
            (
                condition_concept_id IN (SELECT
                    DISTINCT c.concept_id
                FROM
                    `cb_criteria` c
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id
                    FROM
                        `cb_criteria` cr
                    WHERE
                        concept_id IN (317576)
                        AND full_text LIKE '%_rank1]%'      ) a
                        ON (c.path LIKE CONCAT('%.', a.id, '.%')
                        OR c.path LIKE CONCAT('%.', a.id)
                        OR c.path LIKE CONCAT(a.id, '.%')
                        OR c.path = a.id)
                WHERE
                    is_standard = 1
                    AND is_selectable = 1)
                OR  condition_source_concept_id IN (SELECT
                    DISTINCT c.concept_id
                FROM
                    `cb_criteria` c
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id
                    FROM
                        `cb_criteria` cr
                    WHERE
                        concept_id IN (
                            1326588, 1326589, 1326590, 1326591,
                            1569126, 1569127,
                            35207684, 35207685, 35207686, 35207687,
                            35207688, 35207689, 35207691, 35207692,
                            35207693, 35207694, 35207695, 35207696,
                            35207697, 35207698, 35207699, 35207700,
                            35207701, 35207702,
                            44819697, 44819699, 44819700, 44819702,
                            44820857, 44820858, 44820859, 44820860,
                            44820861, 44820862, 44820863, 44820864,
                            44823111, 44824237,
                            44825428, 44825429, 44825430,
                            44826635, 44826636,
                            44827782, 44827783,
                            44828972, 44828973,
                            44830079,
                            44831236, 44831237, 44831238,
                            44832372, 44832373, 44832374, 44832375, 44832376,
                            44833561,
                            44834718, 44834719, 44834720, 44834721,
                            44834723, 44834724, 44834725,
                            44835926, 44835927, 44835928, 44835929,
                            44837099,
                            45533436, 45557536, 45562340,
                            45572079, 45572080, 45576865,
                            45596194, 45605779, 45605781
                        )
                        AND full_text LIKE '%_rank1]%'      ) a
                        ON (c.path LIKE CONCAT('%.', a.id, '.%')
                        OR c.path LIKE CONCAT('%.', a.id)
                        OR c.path LIKE CONCAT(a.id, '.%')
                        OR c.path = a.id)
                WHERE
                    is_standard = 0
                    AND is_selectable = 1)
            )) c_occurrence
    LEFT JOIN
        `concept` c_standard_concept
            ON c_occurrence.condition_concept_id = c_standard_concept.concept_id
", sep = "")

cad_condition_df <- run_bq_export(cad_sql, "cad_conditions", cache_dir = data_dir)

# --- 4b) CAD procedures ---
cat("\n--- Querying CAD procedures ---\n")

cad_procedure_sql <- paste("
    SELECT
        procedure.person_id,
        p_standard_concept.concept_name as standard_concept_name,
        procedure.procedure_datetime
    FROM
        ( SELECT
            *
        FROM
            `procedure_occurrence` procedure
        WHERE
            (
                procedure_source_concept_id IN (SELECT
                    DISTINCT c.concept_id
                FROM
                    `cb_criteria` c
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id
                    FROM
                        `cb_criteria` cr
                    WHERE
                        concept_id IN (
                            2001500, 2001501, 2001502, 2001504,
                            2001507, 2001509, 2001510, 2001511,
                            2001512, 2001513, 2001514, 2001515, 2001517,
                            2106964,
                            2107216, 2107217, 2107218, 2107219,
                            2107220, 2107221, 2107222, 2107223, 2107224,
                            2107230, 2107231,
                            2107242, 2107243, 2107244,
                            2313796,
                            2313801, 2313802, 2313803, 2313804,
                            2313810, 2313811,
                            43527908, 43527909,
                            43527994, 43527995, 43527996,
                            43527997, 43527998, 43527999,
                            43528000, 43528001, 43528002, 43528003, 43528004
                        )
                        AND full_text LIKE '%_rank1]%'      ) a
                        ON (c.path LIKE CONCAT('%.', a.id, '.%')
                        OR c.path LIKE CONCAT('%.', a.id)
                        OR c.path LIKE CONCAT(a.id, '.%')
                        OR c.path = a.id)
                WHERE
                    is_standard = 0
                    AND is_selectable = 1)
            )) procedure
    LEFT JOIN
        `concept` p_standard_concept
            ON procedure.procedure_concept_id = p_standard_concept.concept_id
", sep = "")

cad_procedure_df <- run_bq_export(cad_procedure_sql, "cad_procedures", cache_dir = data_dir)

# --- 4c) CAD from self-report survey ---
cat("  Extracting CAD self-reports from survey...\n")

cad_survey_df <- survey_df %>%
  filter(
    str_detect(answer, regex("coronary artery|heart attack", ignore_case = TRUE)),
    str_detect(answer, regex("yes", ignore_case = TRUE))
  ) %>%
  transmute(
    person_id,
    cad_date = ymd_hms(survey_datetime),
    source = "survey"
  )

cat(sprintf("  CAD survey reports: %d\n", nrow(cad_survey_df)))

# --- 4d) Combine all CAD sources ---
cat("  Combining all CAD sources...\n")

cad_all_df <- bind_rows(
  cad_condition_df %>%
    transmute(
      person_id,
      cad_date = ymd_hms(condition_start_datetime),
      source = "condition"
    ),

  cad_procedure_df %>%
    transmute(
      person_id,
      cad_date = ymd_hms(procedure_datetime),
      source = "procedure"
    ),

  cad_survey_df
)

# Deduplicate: one row per person, earliest date
cad_df <- cad_all_df %>%
  group_by(person_id) %>%
  summarise(
    cad = 1L,
    cad_date = min(cad_date, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  CAD cases: %d unique individuals\n", nrow(cad_df)))

# ============================================================================
# 5) ASSEMBLE FINAL PHENOTYPE FILE
# ============================================================================
cat("\n--- Assembling phenotype file ---\n")

ancestry_df <- read_tsv(file.path(Sys.getenv("HOME"), "workspaces/duplicateofstatinresponse/genome_useful/ancestry_preds.tsv"))

ancestry_parsed_df <- ancestry_df %>%
  # remove [ and ] from the string
  mutate(
    pca_features_clean = str_remove_all(pca_features, "\\[|\\]")
  ) %>%
  # split into 16 PC columns
  separate(
    pca_features_clean,
    into = paste0("PC", 1:16),
    sep = ",\\s*",
    convert = TRUE
  ) %>%
  # make ID match your main dataset
  rename(person_id = research_id) %>%
  # keep just what we need
  select(person_id, ancestry_pred, PC1:PC16)

pheno_df <- demo_df %>%
  select(person_id, sex_at_birth, race, ethnicity, date_of_birth) %>%
  left_join(enrollment_df, by = "person_id") %>%
  mutate(
    age = as.numeric(interval(date_of_birth, enrollment_date) / years(1))
  ) %>%
  left_join(t2d_df, by = "person_id") %>%
  left_join(cad_df, by = "person_id") %>%
  mutate(
    t2d = replace_na(t2d, 0L),
    cad = replace_na(cad, 0L)
  ) %>%
  left_join(ancestry_parsed_df, by = "person_id")

cat(sprintf("  Total individuals: %d\n", nrow(pheno_df)))
cat(sprintf("  T2D cases: %d (%.1f%%)\n",
  sum(pheno_df$t2d), 100 * mean(pheno_df$t2d)))
cat(sprintf("  CAD cases: %d (%.1f%%)\n",
  sum(pheno_df$cad), 100 * mean(pheno_df$cad)))

# Write output
out_file <- file.path(out_dir, "aou_t2d_cad_phenotypes.csv")
write_csv(pheno_df, out_file)

cat(sprintf("\nOutput written to: %s\n", out_file))
cat(sprintf("  Dimensions: %d rows x %d cols\n", nrow(pheno_df), ncol(pheno_df)))
cat(sprintf("  Columns: %s\n", paste(names(pheno_df), collapse = ", ")))
cat("\nHead:\n")
print(head(pheno_df, 10))

cat("\n============================================\n")
cat("Phenotype wrangling complete!\n")
cat("============================================\n")

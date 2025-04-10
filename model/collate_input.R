#!/usr/bin/env Rscript

# Tung Nguyen
# Takes the inputs of each of the iVar and samtools mpileup into one compressed parquet file

# --- Package Loading ---
# Use pacman for easy package management
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cran.r-project.org")
pacman::p_load(argparse, duckdb, stringr, glue, install = TRUE) # Ensure argparse is loaded

# --- Command Line Argument Parsing ---
parser <- ArgumentParser(description = "Combine depth and variant TSV files into compressed Parquet using DuckDB.")

parser$add_argument("--depth-dir", type = "character", required = TRUE,
                    help = "Directory containing input depth TSV files (e.g., '*_depth.tsv')")
parser$add_argument("--variants-dir", type = "character", required = TRUE,
                    help = "Directory containing input variant TSV files (e.g., '*.variants.tsv')")
parser$add_argument("--depth-out", type = "character", required = TRUE,
                    help = "Output file path for the combined depth data (Parquet format)")
parser$add_argument("--variant-out", type = "character", required = TRUE,
                    help = "Output file path for the combined variant data (Parquet format)")
parser$add_argument("--compression-level", type = "integer", default = 6,
                    help = "Zstd compression level for Parquet output [default: %(default)s]")

# Parse arguments
args <- parser$parse_args()

# --- Input Validation ---
if (!dir.exists(args$depth_dir)) {
  stop(glue::glue("Error: Depth directory not found: {args$depth_dir}"))
}
if (!dir.exists(args$variants_dir)) {
  stop(glue::glue("Error: Variants directory not found: {args$variants_dir}"))
}

# --- DuckDB Connection ---
cat("Initializing DuckDB...\n")
# Establish DuckDB connection (in-memory)
con <- dbConnect(duckdb::duckdb())
# Ensure connection is closed even if errors occur
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
cat("DuckDB connection established.\n")

# --- Main Processing ---

# 1. Process Depth Files using DuckDB SQL UNION ALL
cat("\n--- Processing Depth Files ---\n")
cat(glue::glue("Input directory: {args$depth_dir}\n"))
cat(glue::glue("Output file: {args$depth_out}\n"))

depth_files <- list.files(path = args$depth_dir, pattern = "_depth\\.tsv$", full.names = TRUE, recursive = TRUE)

if (length(depth_files) == 0) {
  warning(glue::glue("No depth files found in '{args$depth_dir}' matching pattern '*_depth.tsv'. Skipping depth processing."))
} else {
  cat(glue::glue("Found {length(depth_files)} depth files.\n"))
  cat("Generating DuckDB SQL query for depth files...\n")

  # Prepare column definitions for read_csv (assuming no header)
  col_defs <- "{'gene':'VARCHAR', 'genepos':'VARCHAR', 'depth':'VARCHAR'}"

  # Build a list of SELECT statements, one for each file
  select_statements <- lapply(depth_files, function(file_path) {
    file_basename <- basename(file_path)
    number_name <- stringr::str_remove(file_basename, "_depth\\.tsv$")
    sql_file_path <- gsub("\\\\", "/", file_path) # Ensure forward slashes for SQL

    glue::glue("
      SELECT
          ROW_NUMBER() OVER () AS pos,      -- Generate row number for this file (best effort for order)
          '{number_name}' AS name,          -- Add extracted name
          depth,                            -- Select depth column
          genepos,                          -- Select genepos column
          gene                              -- Select gene column
      FROM read_csv(
          '{sql_file_path}',                -- Read *this specific* file
          delim = '\\t',
          header = false,                   -- No header in depth files
          columns = {col_defs},             -- Define column names and types
          auto_detect = false
      )
    ")
  })

  # Combine all SELECT statements with UNION ALL
  combined_sql <- paste(select_statements, collapse = "\nUNION ALL\n")

  # Create the final COPY query
  copy_depth_query <- glue::glue("
  COPY (
      {combined_sql}
  ) TO '{args$depth_out}' (
      FORMAT PARQUET,
      CODEC 'ZSTD',
      COMPRESSION_LEVEL {args$compression_level}
  );")

  # Execute the query
  cat("Executing DuckDB query to process depth files...\n")
  tryCatch({
    dbExecute(con, copy_depth_query)
    cat(glue::glue("Depth data successfully written to: {args$depth_out}\n"))
  }, error = function(e) {
      # Stop execution on error, provide context
      stop(glue::glue("DuckDB query for depth files failed: {e$message}\n"))
  })
}

# 2. Process Variant Files using DuckDB read_csv_auto
cat("\n--- Processing Variant Files ---\n")
cat(glue::glue("Input directory: {args$variants_dir}\n"))
cat(glue::glue("Output file: {args$variant_out}\n"))

# Create a pattern suitable for DuckDB's globbing
variants_pattern <- file.path(args$variants_dir, "*.variants.tsv")
variants_pattern <- gsub("\\\\", "/", variants_pattern) # Ensure forward slashes

# Check if files exist using R before passing pattern to DuckDB
# This provides a clearer error if the directory is empty or pattern matches nothing
variant_check_files <- list.files(path = args$variants_dir, pattern = "\\.variants\\.tsv$", recursive = TRUE)

if (length(variant_check_files) == 0) {
    warning(glue::glue("No variant files found in '{args$variants_dir}' matching pattern '*.variants.tsv'. Skipping variant processing."))
} else {
    cat(glue::glue("Found {length(variant_check_files)} variant files (using R check).\n"))
    cat(glue::glue("Using DuckDB pattern: {variants_pattern}\n"))

    # Use DuckDB's read_csv_auto directly within a COPY statement
    copy_variant_query <- glue::glue("
    COPY (
        SELECT
            -- Extract base filename for the NAME column
            regexp_replace(filename, '^.*/', '') AS NAME,
            * EXCLUDE(filename) -- Select all original columns from TSV
        FROM read_csv_auto(
            '{variants_pattern}',
            header = true,      -- Variants files have headers
            delim = '\\t',       -- TSV files
            filename = true     -- Need filename column
        )
    ) TO '{args$variant_out}' (
        FORMAT PARQUET,
        CODEC 'ZSTD',
        COMPRESSION_LEVEL {args$compression_level}
    );")

    cat("Executing DuckDB query to process variant files...\n")
    tryCatch({
      dbExecute(con, copy_variant_query)
      cat(glue::glue("Variant data successfully written to: {args$variant_out}\n"))
    }, error = function(e) {
        # Stop execution on error, provide context
        stop(glue::glue("DuckDB query for variant files failed: {e$message}\n"))
    })
}

cat("\n--- Script Finished ---\n")

#' @title Format LGD trap data for SCRAPI / SCRAPI2
#'
#' @description Reads raw biological data from the LGD trapping database (as
#'   returned by \code{\link{query_lgd_trap}} or a saved CSV) and prepares it
#'   for use as the \code{smoltData} argument of \code{\link{SCRAPI}} or
#'   \code{\link{SCRAPI2}}. The main tasks are standardising 6-letter
#'   \code{GenStock} codes and optionally writing the formatted data to a CSV.
#'
#' @param input a data frame of raw trap records, or a path to a CSV file.
#' @param species one of \code{"chnk"} or \code{"sthd"}. Controls which
#'   life-stage label is expected and whether \code{fwAge} is present.
#' @param exportFile optional file path (without \code{.csv} extension). When
#'   supplied the formatted data frame is written to
#'   \code{<exportFile>.csv}.
#' @param check_codes logical. When \code{TRUE} (default), prints a frequency
#'   table of \code{GenStock} values after correction so codes can be verified
#'   before running SCRAPI.
#'
#' @return A data frame ready for use as \code{smoltData} in
#'   \code{\link{SCRAPI}} or \code{\link{SCRAPI2}}, with columns:
#'   \code{MasterID}, \code{CollectionDate}, \code{WeekNumber}, \code{Rear},
#'   \code{GenStock}, \code{GenSex}, \code{BioScaleFinalAge},
#'   \code{LGDFLmm}, \code{GenRun}, \code{SpawnYear}, \code{MPG},
#'   \code{GenStockProb}, \code{GenParentHatchery}, \code{GenBY},
#'   \code{LGDMarkAD}, \code{GenRear}, and (steelhead only) \code{fwAge}.
#'
#' @details
#' \strong{GenStock code corrections applied automatically:}
#'
#' The LGD database occasionally returns 7-character codes; SCOBI / SCRAPI
#' expects exactly 6 characters. The following substitutions are applied:
#'
#' \tabular{ll}{
#'   Raw code  \tab Corrected code \cr
#'   LOWSALM   \tab LOSALM \cr
#'   LOWCLWR   \tab LOCLWR \cr
#'   LOWGRAN   \tab LOGRAN \cr
#'   LOWSNAK   \tab LOSNAK \cr
#'   MIDUPP    \tab MIDUPP \cr
#' }
#'
#' Any code that is not exactly 6 characters after correction triggers a
#' warning listing the offending values.
#'
#' @examples
#' \dontrun{
#' # From a saved CSV
#' chnk <- lgr2SCRAPI("MY2025.CHNKsmolts_trapBioData_allVariables.csv",
#'                     species    = "chnk",
#'                     exportFile = "MY2025CHNK.trapData_formatted")
#'
#' # Directly from query result
#' sthd_raw <- query_lgd_trap("MY2025", srr_prefix = "3", life_stage = "JV")
#' sthd <- lgr2SCRAPI(sthd_raw, species = "sthd",
#'                     exportFile = "MY2025STHD.trapData_formatted")
#' }
#'
#' @export
lgr2SCRAPI <- function(input,
                        species,
                        exportFile  = NULL,
                        check_codes = TRUE) {

  species <- match.arg(species, c("chnk", "sthd"))

  # ---- read input ----------------------------------------------------------
  if (is.character(input)) {
    dat <- read.csv(input, header = TRUE, stringsAsFactors = FALSE)
  } else {
    dat <- as.data.frame(input, stringsAsFactors = FALSE)
  }

  # ---- standardise GenStock codes ------------------------------------------
  corrections <- c(
    LOWSALM = "LOSALM",
    LOWCLWR = "LOCLWR",
    LOWGRAN = "LOGRAN",
    LOWSNAK = "LOSNAK"
  )

  if ("GenStock" %in% names(dat)) {
    bad <- dat$GenStock %in% names(corrections)
    if (any(bad, na.rm = TRUE))
      dat$GenStock[bad] <- corrections[dat$GenStock[bad]]

    # warn on any remaining non-6-character codes (excluding NA)
    non_na   <- dat$GenStock[!is.na(dat$GenStock) & dat$GenStock != "NA"]
    bad_len  <- unique(non_na[nchar(non_na) != 6])
    if (length(bad_len) > 0)
      warning("GenStock codes with length != 6 after correction: ",
              paste(bad_len, collapse = ", "),
              "\nCheck these before running SCRAPI.")
  }

  # ---- select output columns -----------------------------------------------
  core_cols <- c("MasterID", "CollectionDate", "WeekNumber",
                 "Rear", "GenStock", "GenSex",
                 "BioScaleFinalAge", "LGDFLmm",
                 "GenRun", "SpawnYear", "MPG",
                 "GenStockProb", "GenParentHatchery",
                 "GenBY", "LGDMarkAD", "GenRear", "LGDLifeStage")

  if (species == "sthd" && "fwAge" %in% names(dat))
    core_cols <- c(core_cols, "fwAge")

  keep <- intersect(core_cols, names(dat))
  out  <- dat[, keep, drop = FALSE]

  # ---- print code check ----------------------------------------------------
  if (check_codes && "GenStock" %in% names(out)) {
    cat("\nGenStock frequency table after correction (", species, "):\n", sep = "")
    print(sort(table(out$GenStock), decreasing = TRUE))
  }

  # ---- export --------------------------------------------------------------
  if (!is.null(exportFile)) {
    path <- if (grepl("\\.csv$", exportFile)) exportFile else paste0(exportFile, ".csv")
    write.csv(out, path, row.names = FALSE)
    cat("\nFormatted data written to:", path, "\n")
  }

  invisible(out)
}

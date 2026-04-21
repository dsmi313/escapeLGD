#' @title Prepare stratum-level data for Bayesian GE model
#'
#' @description Wrangle raw PIT tag detection data into stratum-level summaries
#'   for the JAGS guidance efficiency model. Two distinct fish pools are created:
#'   (1) a psi estimation pool of upstream-tagged fish classified by their first
#'   route through LGR (GRS vs. UND), and (2) direct LGR counts of all fish
#'   detected at GRJ or GRS during each stratum period. A spill covariate is
#'   attached by joining \code{spill_weekly} to strata by overlapping date ranges.
#'
#' @param dat_up data frame of PIT tag detections. Required columns: \code{tag}
#'   (character), \code{site} (character site code), \code{det_date} (Date).
#' @param strat_assign data frame mapping each calendar date to a trap stratum.
#'   Required columns: \code{date} (Date), \code{stratum} (character or integer),
#'   \code{stratum_idx} (integer, 1-based, must increase monotonically).
#' @param spill_weekly data frame of weekly spill values. Required columns:
#'   \code{week_start} (Date) and \code{spill.per} (spill in raw units, e.g. kcfs).
#' @param downstream_sites character vector of site codes at downstream detection
#'   sites used to identify fish that have passed through LGR.
#'
#' @return A data frame with one row per stratum and columns:
#'   \code{stratum}, \code{stratum_idx}, \code{n_GRS_pool}, \code{n_UND},
#'   \code{n_pool}, \code{n_GRJ_obs}, \code{n_GRS_obs}, \code{spill_val}.
#'
#' @importFrom dplyr filter arrange group_by slice ungroup transmute left_join
#'   mutate case_when full_join summarise distinct count rename across
#' @importFrom tidyr pivot_wider replace_na crossing
#' @export
prep_ge_data <- function(dat_up,
                         strat_assign,
                         spill_weekly,
                         downstream_sites = c("GOJ","LMJ","MCJ","JDJ",
                                              "B2J","BCC","TWX",
                                              "PD5","PD6","PD7","PD8","PDW")) {

  # --- Pool A: psi estimation pool ---
  # Each upstream-tagged fish is classified by its first LGR route (GRS or UND).
  # GRJ fish are excluded — they never entered the spillway, so they cannot inform
  # psi (P(detected at GRS | passed through spillway)).
  down_first <- dat_up %>%
    filter(site %in% downstream_sites) %>%
    arrange(tag, det_date) %>%
    group_by(tag) %>% slice(1) %>% ungroup() %>%
    transmute(tag, down_date = det_date)

  lgr_first <- dat_up %>%
    filter(site %in% c("GRJ", "GRS")) %>%
    arrange(tag, det_date) %>%
    group_by(tag) %>% slice(1) %>% ungroup() %>%
    transmute(tag, lgr_route = site)

  pool_a <- down_first %>%
    left_join(lgr_first, by = "tag") %>%
    mutate(lgr_class = case_when(
      lgr_route == "GRS" ~ "GRS",
      is.na(lgr_route)   ~ "UND",
      TRUE               ~ "GRJ"
    )) %>%
    filter(lgr_class %in% c("GRS", "UND")) %>%
    left_join(strat_assign, by = c("down_date" = "date")) %>%
    filter(!is.na(stratum))

  psi_pool <- pool_a %>%
    group_by(stratum, stratum_idx) %>%
    summarise(n_GRS_pool = sum(lgr_class == "GRS"),
              n_UND      = sum(lgr_class == "UND"),
              n_pool     = n_GRS_pool + n_UND,
              .groups = "drop")

  # --- Pool B: direct LGR counts (one row per fish per stratum) ---
  # distinct(tag, site, stratum) collapses multiple detections of the same fish
  # at the same site within a stratum to a single count.
  lgr_counts <- dat_up %>%
    filter(site %in% c("GRJ", "GRS")) %>%
    distinct(tag, site, det_date) %>%
    left_join(strat_assign, by = c("det_date" = "date")) %>%
    filter(!is.na(stratum)) %>%
    distinct(tag, site, stratum) %>%
    count(stratum, site) %>%
    pivot_wider(names_from = site, values_from = n, values_fill = 0)
  if (!"GRJ" %in% names(lgr_counts)) lgr_counts$GRJ <- 0L
  if (!"GRS" %in% names(lgr_counts)) lgr_counts$GRS <- 0L
  lgr_counts <- rename(lgr_counts, n_GRJ_obs = GRJ, n_GRS_obs = GRS)

  # --- Spill covariate: mean spill.per for calendar weeks overlapping each stratum ---
  strat_dates <- strat_assign %>%
    group_by(stratum) %>%
    summarise(s_min = min(date), s_max = max(date), .groups = "drop")

  spill_strat <- crossing(
      strat_dates,
      mutate(spill_weekly, week_end = week_start + 6)
    ) %>%
    filter(week_start <= s_max, week_end >= s_min) %>%
    group_by(stratum) %>%
    summarise(spill_val = mean(spill.per, na.rm = TRUE), .groups = "drop") %>%
    mutate(spill_val = replace_na(spill_val, 0))

  # --- Combine all stratum summaries ---
  full_join(psi_pool, lgr_counts, by = "stratum") %>%
    left_join(spill_strat, by = "stratum") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    arrange(stratum_idx)
}

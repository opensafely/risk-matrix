# _______________________________________________________________________________________
# Purpose:
# Report vaccination coverage and disease burden in different population subgroups for each campaign
# _______________________________________________________________________________________


# Preliminaries ----

# Import libraries
library("tidyverse")
library("dtplyr")
library("lubridate")
library("glue")
library("here")
library("arrow")
library("survival")
library("splines")
library("parglm")
library("rms")

# Import custom functions
source(here("analysis", "0-lib", "design.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cohort_id <- "rsv_2023_24"
} else {
  cohort_id <- args[[1]]
}

cohort_info <- cohort_info |>
  filter(cohort_id == .env$cohort_id)
# list2env(campaign_info, globalenv())

# Set width of temporal bins for frequencies over time in Kaplan-Meier plots in days
temporal_resolution_km <- 7L

# overwrite cohort end to match rounded Kaplan-Meier curve
cohort_info$cohort_end_days <- ceiling_any(cohort_info$cohort_end_days, temporal_resolution_km)

# Create output directory ----
output_dir <- here("output", "3-report", glue("report_{cohort_id}"))
fs::dir_create(output_dir)
options(width = 200)

# Import prepared data ----
data_prepared <- read_feather(
  here("output", "2-prepare", glue("prepare_{cohort_id}"), glue("prepare_{cohort_id}.arrow"))
)

capture.output(
  skimr::skim_without_charts(data_prepared),
  file = fs::path(output_dir, "data_prepared_skim.txt"),
  split = FALSE
)

#################################################################################################################
# 1. Baseline characteristics
#################################################################################################################

# ________________________________
# 1.A Baseline characteristics 
# ________________________________

baseline_summary <- function(data, ...) {
  group_names <- c(...)
  
  data |>
    count(
      across(all_of(group_names)),
      cohort_id,
      name = "n_group2"
    ) |>
    group_by(
      across(all_of(group_names)),
      cohort_id
    ) |>
    mutate(
      n_group2 = round_sdc(n_group2)
    ) |>
    ungroup()
}

baseline_table <-
  baseline_combos |>
  mutate(
    summary = map2(
      group1,
      group2,
      \(x, y) {
        
        if (is.na(y)) y <- NULL
        
        lookup <- c(
          group1_value = x,
          group2_value = y
        )
        
        baseline_summary(data_prepared, x, y) |>
          mutate(across(c(all_of(c(x, y))), as.character)) |>
          rename(any_of(lookup))
      }
    )
  ) |>
  unnest(summary) |>
  select(
    cohort_id,
    group1,
    group1_value,
    group2,
    group2_value,
    everything()
  )

# Add group1 totals
baseline_table <- baseline_table |>
  mutate(group2_value = replace_na(group2_value, "All")) |>
  left_join(
    baseline_table |> filter(group2 == "all") |> select(group1, group1_value, n_group1 = n_group2),
    by = c("group1", "group1_value")
  ) |>
  relocate(n_group1, .after = group1_value) |>
  mutate(pct = round(100 * n_group2 / n_group1, 1))

write_csv(
  baseline_table,
  fs::path(output_dir, "baseline_summary.csv")
)


# ________________________________
# 1.B Baseline vaccination history 
# ________________________________

baseline_vax_summary <- function(data, ...) {
  group_names <- c(...)

  data |>
    count(
      across(all_of(group_names)),
      cohort_id,
      baseline_vax_status,
      name = "n_round10"
    ) |>
    group_by(
      across(all_of(group_names)),
      cohort_id
    ) |>
    mutate(
      n_total_round10 = sum(n_round10),
      n_total_round10 = round_sdc(n_total_round10),
      n_round10 = round_sdc(n_round10),
      pct = round(100 * n_round10 / n_total_round10, 1)
    ) |>
    ungroup()
}

baseline_vax_table <-
  level_combos |>
  mutate(
    summary = map2(
      group1,
      group2,
      \(x, y) {

        if (is.na(y)) y <- NULL

        lookup <- c(
          group1_value = x,
          group2_value = y
        )

        baseline_vax_summary(data_prepared, x, y) |>
          mutate(across(c(all_of(c(x, y))), as.character)) |>
          rename(any_of(lookup))
      }
    )
  ) |>
  unnest(summary) |>
  select(
    cohort_id,
    group1,
    group1_value,
    group2,
    group2_value,
    everything()
  )

write_csv(
  baseline_vax_table,
  fs::path(output_dir, "baseline_vax_history.csv")
)

#################################################################################################################
# 2. Follow-up descriptive analyses
#################################################################################################################

# ________________________________
# 2.A Vaccination uptake
# Kaplan-Meier estimates of cumulative vaccination coverage
# ________________________________

# This code borrows heavily from the KM reusable action https://github.com/opensafely-actions/kaplan-meier-function/blob/main/analysis/km.R
# The reason not to use the KM reusable action directly is that we need to reuse it multiple times across many different stratification variables
# This would create a large number of project.yaml actions, so it's easier to do the repeats within a single script

## tests ----

# follow-up times used as the log offset in the models (narrow outcomes exist for RSV only)
modelled_times <- intersect(
  c("vax_time", "primary_care_time", "admitted_time", "admitted_primary_time",
    "disease_death_time", "admitted_narrow_time", "admitted_primary_narrow_time",
    "disease_death_narrow_time"),
  names(data_prepared)
)

walk(
  modelled_times,
  \(time_var) {
    times_count <- table(cut(data_prepared[[time_var]], c(-Inf, 0, 1, Inf), right = FALSE, labels = c("<0", "0", ">0")), useNA = "ifany")
    
    if (!identical(as.integer(times_count), c(0L, 0L, nrow(data_prepared)))) {
      print(times_count)
      stop(glue("all {time_var} values must be strictly positive"))
    }
  }
)

## Function to calculate KM estimates for a given stratification variable ----

# group variables are provided as characters via dots (...)
# resolution argument is the precision used for the time dimension. If zero, then original resolution is used.
km_estimates <- function(data, group_name1, group_name2, event_name, event_time, event_indicator, resolution = 0) {

  group_names <- c(group_name1, group_name2)

  data_outcome <-
    data |>
    select(
      patient_id,
      all_of(group_names),
      event_time = any_of(event_time),
      event_indicator = any_of(event_indicator)
    ) |>
    mutate(event_time = ceiling_any(event_time, resolution))

  data_km <-
    data_outcome |>
    group_by(across(all_of(group_names))) |>
    nest() |>
    mutate(
      surv_obj_tidy = map(
        data, ~ {
          survfit(
            Surv(event_time, event_indicator) ~ 1,
            data = .x,
            conf.type = "log-log"
          ) |>
            broom::tidy() |>
            add_row(
              time = 0, # assumes time origin is zero
              n.risk = 0,
              n.event = 0,
              n.censor = 0,
              estimate = 1,
              conf.low = 1,
              conf.high = 1,
              .before = 1L
            ) |>
            complete(
              time = seq(0L, cohort_info$cohort_end_days, resolution), # fill in 1 row for each period (defined by resolution) of follow up
              fill = list(n.event = 0L, n.censor = 0L) # fill in zero events on those days
            ) |>
            fill(
              n.risk,
              .direction = "up"
            ) |>
            fill(
              estimate, conf.low, conf.high,
              .direction = "down"
            )
        }
      ),
    ) |>
    select(-data) |>
    unnest(surv_obj_tidy) |>
    mutate(

      # disclosure control
      n.risk = ceiling_any(n.risk, sdc_precision),
      estimate = plyr::round_any(estimate, sdc_precision / nth(n.risk, 2)), # use 2nd value as this skips the t=0 row where n.risk=0
      conf.low = plyr::round_any(conf.low, sdc_precision / nth(n.risk, 2)),
      conf.high = plyr::round_any(conf.high, sdc_precision / nth(n.risk, 2)),

      # cumulative incidence
      cmlinc = 1 - estimate,
      cmlinc.low = 1 - conf.high,
      cmlinc.high = 1 - conf.low,
    )

  rm(data_outcome)

  # print(data_km)

  # write tables that capture underlying plotting data
  data_km_nozero <-
    data_km |>
    arrange(across(all_of(group_names))) |>
    ungroup() |>
    filter(time != 0) |>
    select(
      all_of(group_names),
      time,
      cmlinc,
      cmlinc.low,
      cmlinc.high,
    )

  # commented out because outputting the entire data creates a dataset that is too large to output-check on it's own
  # write_csv(data_km_nozero, fs::path(output_dir, glue("km_{event_name}_{paste0(group_names, collapse='_')}.csv")))
  return(data_km_nozero)

}

# for testing function interactively
# km_estimates_vax <- partial(
#   km_estimates, data = data_prepared, event_name = "vax", event_time = "vax_time", event_indicator = "vax_indicator", resolution = temporal_resolution_km
# )
# km_estimates_vax("all", "all")
# km_estimates_vax("all", "ageband6")
# km_estimates_vax("ageband6", "all")
# km_estimates_vax("ageband6", "sex")

## Function to get km estimates over all subgroup combinations ----

get_all_km_estimates <- function(data, event_name, event_time, event_indicator, resolution) {
  # loop over all group1 and group2 variable combinations and combine into one big dataset
  km_estimates_table <-
    level_combos |>
    mutate(
      km_summary = map2(
        group1, group2,
        .f = \(x, y) {

          if (is.na(y)) y <- NULL
          lookup <- c(group1_value = x, group2_value = y)

          km_estimates(data = data, group_name1 = x, group_name2 = y, event_name = event_name, event_time = event_time, event_indicator = event_indicator, resolution = temporal_resolution_km) |>
            mutate(across(c(all_of(c(x, y))), as.character)) |>
            rename(all_of(c(group1_value = x, group2_value = y)))
        }
      )
    ) |>
    unnest(km_summary) |>
    select(group1, group1_value, group2, group2_value, everything()) |>
    mutate(
      cohort_end_milestone = (time == cohort_info$cohort_end_days) * 1L,
    )

  # Write table to CSV files, limited to group2 == all
  km_estimates_table_sub <- km_estimates_table |> filter(group2=="all")
  # iwalk(
  #    split(km_estimates_table_sub, km_estimates_table_sub$group1),
  #    ~ write_csv(.x, fs::path(output_dir, glue("km_estimates_{event_name}_table_{.y}.csv")))
  # )
  write_csv(
    km_estimates_table_sub,
    fs::path(output_dir, glue("km_estimates_{event_name}_table.csv"))
  )
  
  # Write table to a CSV file containing _only_ reporting milestones
  km_estimates_milestones <-
    km_estimates_table |>
    filter(cohort_end_milestone == 1L) |>
    mutate(
      milestone_date = cohort_info$cohort_end_date,
      milestone = "cohort_end"
    ) |>
    rename(days_since_cohort_start = time) |>
    select(-cohort_end_milestone)

  # Code to split outputs if 5,000 rows exceeded
  # iwalk(
  #   split(km_estimates_milestones, km_estimates_milestones$milestone),
  #   ~ write_csv(.x, fs::path(output_dir, glue("km_estimates_{event_name}_milestones_{.y}.csv")))
  # )
  write_csv(
    km_estimates_milestones,
    fs::path(output_dir, glue("km_estimates_{event_name}_milestones_cohort_end.csv"))
  )
  return(km_estimates_table)
}

## _______________________________________________________________________________________
cat("## Report KM cumulative incidence of vaccination in a standardised table", "\n")
## _______________________________________________________________________________________

km_estimates_vax_table <- get_all_km_estimates(data = data_prepared, event_name =  "vax", event_time = "vax_time", event_indicator = "vax_indicator")
#km_estimates_vax_alive_table <- get_all_km_estimates(filter(data_prepared, alive_and_registered), "vax_alive", event_time = "vax_time", event_indicator = "vax_indicator")

# ________________________________
# 2.A.i Kaplan-Meier plots (for server inspection only)
# ________________________________
## function to make KM plots for the data ----

km_plot <- function(km_data, event_name, group1, group2) {

  group_names <- c(group1, group2)

  # plot km curves locally for checking (but probs not for release as these can be reconstructed from released data)
  coverage_plot <-
    bind_rows(
      # this bit adds an extra row of data so that the firsrt line of the KM plot shows
      km_data |> summarise(
        time = 0,
        cmlinc = 0, cmlinc.low = 0, cmlinc.high = 0, lagtime = 0,
        .by = c("group1", "group1_value", "group2", "group2_value")
      ),
      km_data |> mutate(
        lagtime = lag(time, 1L, 0L), # assumes the time-origin is zero

        .by = c("group1", "group1_value", "group2", "group2_value")
      )
    ) |>
    ggplot() +
    geom_step(aes(x = time, y = cmlinc, group = group2_value, colour = group2_value), direction = "vh") +
    geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high, group = group2_value, fill = group2_value), alpha = 0.1, colour = "transparent") +
    facet_grid(rows = "group1_value") +
    scale_color_brewer(type = "qual", palette = "Set1", na.value = "grey") +
    scale_fill_brewer(type = "qual", palette = "Set1", guide = "none", na.value = "grey") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(xlim = c(0, NA)) +
    labs(
      x = "Days since cohort start",
      y = "Cumulative vaccination coverage",
      colour = NULL,
      title = NULL
    ) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(colour = "black"),
      panel.grid.minor.x = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(.05, .95),
      legend.justification = c(0, 1),
    )

  ggsave(fs::path(output_dir, glue("km_{event_name}_{paste0(group_names, collapse='_')}.png")), plot = coverage_plot)
  # print(coverage_plot)
}

## _______________________________________________________________________________________
cat("## Plot KM cumulative incidences of for inspection the server", "\n")
## _______________________________________________________________________________________

# plot km estimates for inspection the server
walk2(
  level_combos$group1, level_combos$group2,
  .f = \(x, y) {
    km_plot(
      km_data = km_estimates_vax_table |> filter(group1 == x, group2 == y),
      event_name = "vax",
      group1 = x,
      group2 = y
    )
  }
)

#################################################################################################################
# 2.B Disease burden
# Incidence of disease-specific outcomes
#################################################################################################################

# Models:
#   1 = unadjusted (primary estimate)
#   2 = age (restricted cubic spline, 3 knots: rcs(age, 3)) + sex
#   3 = model 2 + all OTHER core clinical subgroups (core subgroups only; for all other
#       subgroups model 3 rows are output with NA estimates so every variable has the
#       same set of rows/columns)

# list of reference categories, which overwrites default reference category if needed
contrasts_reference_levels <- list(
  ageband6 = "65-84",
  ageband12 = "65-74"
)

## Build the model formula for a given subgroup and model number
build_model_spec <- function(data, subgroup, model) {
  
  terms <- subgroup
  
  if (model >= 2) {
    
    # add sex unless sex is the subgroup of interest or does not vary
    if (subgroup != "sex" && n_distinct(data$sex) > 1) {
      terms <- c(terms, "sex")
    }
    
    # add age as restricted cubic spline with 3 knots
    # knots at the 10th, 50th and 90th percentiles of age among individuals in this stratum
    # (the rms default for 3 knots but computed here explicitly)
    # dropped if age is the subgroup, or if knots would coincide (too little age variation)
    if (!(subgroup %in% age_subgroups)) {
      age_knots <- round(quantile(data$age, probs = c(0.1, 0.5, 0.9), names = FALSE), 1)
      if (n_distinct(age_knots) == 3) {
        terms <- c(
          terms,
          glue("rcs(age, c({age_knots[1]}, {age_knots[2]}, {age_knots[3]}))")
        )
      }
    }
  }

  if (model == 3) {
    # all other core subgroups, dropping any that don't vary in this stratum
    other_core <- setdiff(core_subgroups, subgroup)
    other_core <- other_core[map_lgl(other_core, \(v) n_distinct(data[[v]]) > 1)]
    terms <- c(terms, other_core)
  }
  
  list(
    formula = as.formula(paste("event_indicator ~", paste(terms, collapse = " + "))),
    covariates = paste(setdiff(terms, subgroup), collapse = "; ")
  )
}

## Function to output IRs and IRRs for disease burden comparing levels of one subgroup ----
adjusted_estimates <- function(data, subgroup, event_time, event_indicator, model = 1) {

  # model 3 is only defined for core clinical subgroups
  # for other subgroups, still output descriptive rows with NA estimates
  model_applicable <- !(model == 3 && !(subgroup %in% core_subgroups))
  
  if (model_applicable) {
    spec <- build_model_spec(data, subgroup, model)
  } else {
    spec <- list(formula = NULL, covariates = NA_character_)
  }
  
  # variables needed for this model
  model_vars <- unique(c(subgroup, "sex", "age", if (model == 3 && model_applicable) core_subgroups))
  
  # aggregate to covariate patterns: people, events, person-time
 
  # prepare dataset
  data_outcome <-
    data |>
    mutate(
      event_time = .data[[event_time]],
      event_indicator = .data[[event_indicator]]
    ) |>
    select(
      all_of(model_vars),
      event_time,
      event_indicator
    ) |>
    # to save memory, reduce dataset size by counting all distinct rows
    # then use these counts as weights when modelling / summarising later
    summarise(
      count = n(),
      .by = all_of(c(model_vars, "event_time", "event_indicator"))
      )
  
  # how many possible values of the group are there
  n_values <- n_distinct(data_outcome[[subgroup]])
  
  # all levels of variable that exist in the data, sorted by factor (or alphanumeric if not)
  # the first of these is the default reference level
  all_levels <- unique(as.character(sort(data_outcome[[subgroup]])))

  # explicitly include reference_level, and define contrasts if unusual
  # for example, if we wanted to use not the default reference level for a factor
  # this object is used inside a glm() call, to be passed to the `contrasts` argument
  if (
    (subgroup %in% names(contrasts_reference_levels)) &
      (n_values > 0) &
      ifelse(is.null(contrasts_reference_levels[[subgroup]]), FALSE, contrasts_reference_levels[[subgroup]] %in% all_levels)
  ) {
    # create the object passed to `contrasts` argument in glm call to use a different reference category
    subgroup_contrasts <- list(contr.treatment(all_levels,  which(all_levels == contrasts_reference_levels[[subgroup]])))
    names(subgroup_contrasts) <- subgroup
    reference_level <- contrasts_reference_levels[[subgroup]]
  } else {
    subgroup_contrasts <- NULL
    reference_level <- all_levels[1]
  }
  
  # -------------------------------
  # Descriptive summary
  # summarise total people, events, and person-time
  # and add contrast label for merging with model output later
  #--------------------------------
  data_summary <-
    data_outcome |>
    mutate(label = .data[[subgroup]]) |>
    arrange(label) |>
    summarise(
      variable = subgroup,
      n_obs = round_sdc(sum(count)),
      n_event = round_sdc(sum(event_indicator * count)),
      exposure = round_sdc(sum(event_time * count)),
      .by = label
    ) |>
    mutate(
      label = as.character(label),
      reference_row = label == reference_level,
      contrast = glue("{subgroup}{label}")
    )

  # -------------------------------
  # Poisson regression
  # Incidence rate ratios (IRR)
  # -------------------------------

  # all-NA estimate frame, with exactly the columns a successful model fit returns.
  # used wherever a model is not fitted, so that every subgroup x model combination
  # returns an identically shaped result
  na_estimates <-
    data_summary |>
    select(variable, label) |>
    mutate(
      estimate = NA_real_,
      std.error = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    )
  
  
  if (n_values > 1 && model_applicable) {

    parglm_control <- parglm.control(maxit = 40, nthreads = 4)

    # fit the model
    # if there is an error, just return an empty dataset rather than fail
    data_poisson0 <-
      tryCatch(
        expr = {
          data_outcome |>
            parglm(
              data = _,
              formula = spec$formula,
              family = poisson,
              offset = log(event_time),
              control = parglm_control,
              weights = count,
              contrasts = subgroup_contrasts
            ) |>
            broom.helpers::tidy_and_attach(tidy_fun = broom.helpers::tidy_parameters, ci_method = "wald") |>
            broom.helpers::tidy_add_term_labels() |>
            filter(variable == subgroup) |>
            select(variable, label, estimate, std.error, conf.low, conf.high)
          # note: the filter above is the same as doing marginaleffects::avg_comparisons(model, type = "link", variables = subgroup, comparison = "difference"),
          # as long as there are no interaction terms between subgroup and anything else
          # we use broom.helpers functions because it gives us the really nice variable and label info formatting for the outputted tidy dataset
          # if we want to use avg_comparisons in future, then attach the nicely formatted meta info onto a broom::tidy(avg_comparisons) object
          # or see "get_estimates_using_marginaleffects.R" script for a clue
        },
        error = function(e) {
          cat("error for subgroup", subgroup, "model", model, ":", conditionMessage(e), "\n")
          na_estimates
          }
      )
  } else {
    data_poisson0 <- na_estimates
  }

  # combine summary and model outputs
  data_joined <-
    full_join(
      data_summary,
      data_poisson0,
      by = c("variable", "label"),
    ) 
  
  # incidence rates per 1,000 person-years, calculated once
  ir_calc <- calc_ir(data_joined$n_event, data_joined$exposure)
  
  data_poisson <- data_joined |>
    transmute(
      variable, label, reference_row,
      n_obs, n_event, exposure,
      ir_1000py = ir_calc$ir,
      ir_1000py_lower = ir_calc$ir_lower_ci,
      ir_1000py_upper = ir_calc$ir_upper_ci,
      irr = exp(estimate),
      irr.low = exp(conf.low),
      irr.high = exp(conf.high),
      irr.ln.std.error = std.error,
      model = as.integer(model),
      covariates = spec$covariates
    )

  rm(data_outcome)
  gc()
  return(data_poisson)
}

# adjusted_estimates(data_prepared, "ageband7", "admitted_time", "admitted_indicator")

#################################################################################################################
# 3. Comparative analyses
#################################################################################################################

# 3.A Poisson regression
# Incidence rate ratios (IRRs)

## function to loop over all group combinations and report IRR of level0 versus levelX ----
# for a given outcome, loop over all groups combinations, obtaining contrasts for each using adjusted_estimates function, and combining into one file
# specifically, compare level2 groups amongst each other, for all people meeting level1 group criteria
get_all_estimates <- function(data, event_name, event_time, event_indicator, models = 1:3) {

  l1ticker <<- ""

  estimates_list <-
    level_combos |>
    mutate(
      estimates = map2(
        group1, group2,
        function(group1, group2) {
          # just to keep track of how far along we are
          if (l1ticker != group1) {
            print(group1)
            l1ticker <<- group1
          }

          summary_data <-
            data |>
            mutate(
              label1 = data[[group1]],
            ) |>
            nest(.by = c(label1), .key = "group1_subset") |>
            mutate(
              estimates = map(group1_subset, \(group1_subset) {
                map(models, \(m) {
                  adjusted_estimates(group1_subset, group2, event_time, event_indicator, model = m)
                }) |>
                  bind_rows()
              })
            ) |>
            select(-group1_subset) |>
            unnest(estimates) |>
            select(-variable) |>
            rename(label2 = label) |>
            mutate(
              across(c(label1, label2), as.character) # to ensure the unnest() works later
            )

          return(summary_data)
        }
      )
    ) |>
    unnest(estimates) |>
    select(group1, label1, group2, label2, model, everything()) # reorder columns

  # reshape to one row per subgroup level, with model-specific columns side by side
  # descriptive columns (counts, person-time, IR) are identical across models, so they
  # identify the row; only the model-specific columns are spread out with an _m1/_m2/_m3 suffix
  model_cols <- c("irr", "irr.low", "irr.high", "irr.ln.std.error", "covariates")
  
  estimates_wide <-
    estimates_list |>
    pivot_wider(
      names_from = model,
      values_from = all_of(model_cols),
      names_glue = "{.value}_m{model}",
      names_vary = "slowest"   # group columns by model
    ) |>
    select(-covariates_m1)     # always empty

  # check the reshape didn't split any row (i.e. descriptive columns really were identical)
  stopifnot(nrow(estimates_wide) * length(models) == nrow(estimates_list))
  
  # redact statistics derived from groups of 5 or fewer, which round_sdc() maps to zero.
  # grouped by contrast, because an IRR discloses about every level it compares
  estimates_wide <-
    estimates_wide |>
    mutate(redacted_obs = n_obs == 0) |>
    mutate(across(c(exposure, starts_with("ir_1000py")), \(x) if_else(redacted_obs, NA, x))) |>
    mutate(redacted = any(n_obs == 0 | n_event == 0), .by = c(group1, label1, group2)) |>
    mutate(across(starts_with("irr"), \(x) if_else(redacted, NA, x)))
  
  write_csv(estimates_wide, fs::path(output_dir, glue("contrasts_{event_name}.csv")))

}

## _______________________________________________________________________________________
cat("## Get all IRR contrasts for vaccination and burden", "\n")
## _______________________________________________________________________________________

# IRR for vaccination, in the usual way
get_all_estimates(data_prepared, "vax", "vax_time", "vax_indicator")

# IRR for vaccination, only looking at those who survived or stayed registered to the end of the season (collider bias! but matches UKHSA reports)
# get_all_estimates(filter(data_prepared, alive_and_registered), "vax_alive", "vax_time", "vax_indicator")

# IRR for burden, in the usual way
get_all_estimates(data_prepared, "primary_care", "primary_care_time", "primary_care_indicator")
get_all_estimates(data_prepared, "admitted", "admitted_time", "admitted_indicator")
get_all_estimates(data_prepared, "admitted_primary", "admitted_primary_time", "admitted_primary_indicator")
get_all_estimates(data_prepared, "disease_death", "disease_death_time", "disease_death_indicator")

# Exclusion of immunosuppressed
get_all_estimates(data_prepared |> filter(immunosuppressed==0), "primary_care_no_IS", "primary_care_time", "primary_care_indicator")
get_all_estimates(data_prepared |> filter(immunosuppressed==0), "admitted_no_IS", "admitted_time", "admitted_indicator")
get_all_estimates(data_prepared |> filter(immunosuppressed==0), "admitted_primary_no_IS", "admitted_primary_time", "admitted_primary_indicator")
get_all_estimates(data_prepared |> filter(immunosuppressed==0), "disease_death_no_IS", "disease_death_time", "disease_death_indicator")

# Exclusion of immunosuppressed_crd
get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "primary_care_no_IS_CRD", "primary_care_time", "primary_care_indicator")
get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "admitted_no_IS_CRD", "admitted_time", "admitted_indicator")
get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "admitted_primary_no_IS_CRD", "admitted_primary_time", "admitted_primary_indicator")
get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "disease_death_no_IS_CRD", "disease_death_time", "disease_death_indicator")

# RSV only: narrow phenotype for severe outcomes
if ("admitted_narrow_indicator" %in% names(data_prepared)) {
  get_all_estimates(data_prepared, "admitted_narrow", "admitted_narrow_time", "admitted_narrow_indicator")
  get_all_estimates(data_prepared, "admitted_primary_narrow", "admitted_primary_narrow_time", "admitted_primary_narrow_indicator")
  get_all_estimates(data_prepared, "disease_death_narrow", "disease_death_narrow_time", "disease_death_narrow_indicator")
  
  get_all_estimates(data_prepared |> filter(immunosuppressed==0), "admitted_narrow_no_IS", "admitted_narrow_time", "admitted_narrow_indicator")
  get_all_estimates(data_prepared |> filter(immunosuppressed==0), "admitted_primary_narrow_no_IS", "admitted_primary_narrow_time", "admitted_primary_narrow_indicator")
  get_all_estimates(data_prepared |> filter(immunosuppressed==0), "disease_death_narrow_no_IS", "disease_death_narrow_time", "disease_death_narrow_indicator")
  
  get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "admitted_narrow_no_IS_CRD", "admitted_narrow_time", "admitted_narrow_indicator")
  get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "admitted_primary_narrow_no_IS_CRD", "admitted_primary_narrow_time", "admitted_primary_narrow_indicator")
  get_all_estimates(data_prepared |> filter(immunosuppressed_crd==0), "disease_death_narrow_no_IS_CRD", "disease_death_narrow_time", "disease_death_narrow_indicator")
}

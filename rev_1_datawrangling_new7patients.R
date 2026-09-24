# analysis Revision 1 - new 7 patients addition
# LIBRARIES ----
pacman::p_load(
  ggplot2,
  readxl,
  Hmisc,
  gtsummary,
  summarytools,
  kableExtra,
  expss,
  rmarkdown,
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable,     # converting tables to pretty images
  ggstatsplot,
  palmerpenguins,
  survival,
  DiagrammeR,
  glue,
  igraph,
  ggraph,
  progressr,
  zoo,
  changepoint,
  ggpubr,
  MASS,
  mice,
  nnet,
  DescTools,
  officer,
  dplyr,
  pROC,
  patchwork
)
# LIBRARIES ENDE----


# call usual functions
db.source.file.path <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "c:/sci/cdk.stat/", "C:/app/r.apps/cdk.stat/")
source(paste0(db.source.file.path, "sp.utilities.r"))
source(paste0(db.source.file.path, "nirs_helper.r"))

db.language <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "Meine Ablage", "My Drive")
db.path_local <- paste0("C:/sci/rcode/nirs/")

# DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDB
# LAST STABLE DB
#db_clin <- read.csv2(paste0(db.path_local, "clinical_new_clean_v2.csv"))
db_clin <- read.csv2(paste0(db.path_local, "clinical_nirs_24092026_1.csv"))
# load main database
#db_msr <- read.csv2(paste0(db.path_local, "nirs_measurements_clean_v2.csv"))
db_msr <- read.csv2(paste0(db.path_local, "nirs_24092026.csv"))

# load new 7 patients
db_new <- readxl::read_xlsx("G:/Meine Ablage/sci/dbs/graz/nirs/nirs_addition_23092026.xlsx")

# save new version

write.csv2(db_clin, paste0(db.path_local, "clinical_nirs_24092026_4.csv"))
write.csv2(db_msr, paste0(db.path_local, "nirs_24092026_1.csv"))


# Reshape ------------------------------------------------------------
# 0. Reshape Left/Right measurements to long format
# 1. Combine date + time into datetime
# 2. Calculate sample interval and elapsed seconds

db_new_long <- db_new |>
  
  mutate(
    patient_id = as.integer(patient_id),
    
    # combine real calendar date with clock time
    datetime = as.POSIXct(
      paste(
        format(date, "%Y-%m-%d"),
        format(time, "%H:%M:%S")
      ),
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Europe/Vienna"
    )
  ) |>
  
  pivot_longer(
    cols = c(value_L, value_R),
    names_to = "hemisphere_side",
    values_to = "nirs_value_raw"
  ) |>
  
  mutate(
    hemisphere_side = recode(
      hemisphere_side,
      value_L = "Left",
      value_R = "Right"
    )
  ) |>
  
  arrange(
    patient_id,
    hemisphere_side,
    datetime
  ) |>
  
  group_by(
    patient_id,
    hemisphere_side
  ) |>
  
  mutate(
    sample_interval_seconds = as.numeric(
      difftime(
        datetime,
        lag(datetime),
        units = "secs"
      )
    ),
    
    elapsed_seconds = as.numeric(
      difftime(
        datetime,
        first(datetime),
        units = "secs"
      )
    )
  ) |>
  
  ungroup()


n_distinct(db_new_long$patient_id)
# should be 7

nrow(db_new_long)
# should be 2 * 6298 = 12596

range(db_new_long$datetime, na.rm = TRUE)

table(db_new_long$hemisphere_side)

timing_check <- db_new_long |>
  group_by(patient_id, hemisphere_side) |>
  summarise(
    n = n(),
    median_interval = median(sample_interval_seconds, na.rm = TRUE),
    q1_interval = quantile(sample_interval_seconds, 0.25, na.rm = TRUE),
    q3_interval = quantile(sample_interval_seconds, 0.75, na.rm = TRUE),
    min_interval = min(sample_interval_seconds, na.rm = TRUE),
    max_interval = max(sample_interval_seconds, na.rm = TRUE),
    n_zero_interval = sum(sample_interval_seconds == 0, na.rm = TRUE),
    n_negative_interval = sum(sample_interval_seconds < 0, na.rm = TRUE),
    .groups = "drop"
  )

view(timing_check)

zero_report <- db_new_long |>
  group_by(patient_id, hemisphere_side) |>
  summarise(
    n_measurements = n(),
    
    n_zero = sum(nirs_value_raw == 0, na.rm = TRUE),
    pct_zero = 100 * n_zero / n_measurements,
    
    n_missing = sum(is.na(nirs_value_raw)),
    pct_missing = 100 * n_missing / n_measurements,
    
    n_zero_or_missing = sum(
      is.na(nirs_value_raw) | nirs_value_raw == 0
    ),
    
    pct_zero_or_missing =
      100 * n_zero_or_missing / n_measurements,
    
    .groups = "drop"
  )

zero_report

zero_report_patient <- db_new_long |>
  group_by(patient_id) |>
  summarise(
    n_measurements = n(),
    n_zero = sum(nirs_value_raw == 0, na.rm = TRUE),
    n_missing = sum(is.na(nirs_value_raw)),
    n_zero_or_missing =
      sum(is.na(nirs_value_raw) | nirs_value_raw == 0),
    pct_zero_or_missing =
      100 * n_zero_or_missing / n_measurements,
    .groups = "drop"
  )

zero_report_patient

max_interp_gap_seconds <- 30

interpolate_short_gaps <- function(
    time,
    value,
    max_gap_seconds = 30
) {
  
  out <- value
  
  missing_idx <- which(is.na(value))
  valid_idx   <- which(!is.na(value))
  
  # cannot interpolate with fewer than 2 valid observations
  if (length(valid_idx) < 2 || length(missing_idx) == 0) {
    return(out)
  }
  
  for (i in missing_idx) {
    
    left_candidates  <- valid_idx[valid_idx < i]
    right_candidates <- valid_idx[valid_idx > i]
    
    # no extrapolation at beginning/end
    if (
      length(left_candidates) == 0 ||
      length(right_candidates) == 0
    ) {
      next
    }
    
    left_i  <- max(left_candidates)
    right_i <- min(right_candidates)
    
    gap_seconds <- time[right_i] - time[left_i]
    
    # interpolate only if gap is short enough
    if (
      is.finite(gap_seconds) &&
      gap_seconds <= max_gap_seconds
    ) {
      
      out[i] <-
        value[left_i] +
        (
          value[right_i] - value[left_i]
        ) *
        (
          (time[i] - time[left_i]) /
            (time[right_i] - time[left_i])
        )
    }
  }
  
  out
}

db_new_long <- db_new_long |>
  mutate(
    # preserve original raw measurement
    nirs_value_clean = if_else(
      is.na(nirs_value_raw) | nirs_value_raw == 0,
      NA_real_,
      as.numeric(nirs_value_raw)
    )
  ) |>
  
  group_by(patient_id, hemisphere_side) |>
  
  mutate(
    nirs_value_interpolated =
      interpolate_short_gaps(
        time = elapsed_seconds,
        value = nirs_value_clean,
        max_gap_seconds = max_interp_gap_seconds
      ),
    
    was_interpolated =
      is.na(nirs_value_clean) &
      !is.na(nirs_value_interpolated)
  ) |>
  
  ungroup()

interpolation_report <- db_new_long |>
  group_by(patient_id, hemisphere_side) |>
  summarise(
    n_total = n(),
    n_original_valid = sum(!is.na(nirs_value_clean)),
    n_invalid_or_missing = sum(is.na(nirs_value_clean)),
    n_interpolated = sum(was_interpolated),
    n_remaining_missing =
      sum(is.na(nirs_value_interpolated)),
    
    pct_interpolated =
      100 * n_interpolated / n_total,
    
    pct_remaining_missing =
      100 * n_remaining_missing / n_total,
    
    .groups = "drop"
  )

interpolation_report

db_new_ready <- db_new_long |>
  transmute(
    patient_id = as.integer(patient_id),
    datetime = datetime,
    
    # not known until db_clin is merged
    hemisphere_status = NA_character_,
    
    # preserve anatomical side for later mapping
    hemisphere_side = hemisphere_side,
    
    nirs_value_raw = as.numeric(nirs_value_raw),
    nirs_value_interpolated =
      as.numeric(nirs_value_interpolated),
    
    sample_interval_seconds =
      as.numeric(sample_interval_seconds),
    
    elapsed_seconds =
      as.numeric(elapsed_seconds),
    
    angiography_start_point = FALSE,
    recanalisation_point = FALSE,
    first_try_point = FALSE,
    
    source = "new_7"
  )

db_msr_clean <- db_msr |>
  transmute(
    patient_id = as.integer(patient_id),
    
    datetime = as.POSIXct(
      datetime,
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Europe/Vienna"
    ),
    
    hemisphere_status = hemisphere_status,
    
    # old data no longer contain anatomical L/R information
    hemisphere_side = NA_character_,
    
    nirs_value_raw =
      as.numeric(nirs_value_raw),
    
    nirs_value_interpolated =
      as.numeric(nirs_value_interpolated),
    
    sample_interval_seconds =
      as.numeric(sample_interval_seconds),
    
    elapsed_seconds =
      as.numeric(elapsed_seconds),
    
    angiography_start_point =
      angiography_start_point,
    
    recanalisation_point =
      recanalisation_point,
    
    first_try_point =
      first_try_point,
    
    source = "original"
  )

db_msr_combined <- bind_rows(
  db_msr_clean,
  db_new_ready
) |>
  arrange(
    patient_id,
    datetime,
    hemisphere_side,
    hemisphere_status
  ) |>
  mutate(
    row_id = row_number()
  )

nrow(db_msr_clean)
nrow(db_new_ready)
nrow(db_msr_combined)

n_distinct(db_msr_clean$patient_id)
n_distinct(db_new_ready$patient_id)
n_distinct(db_msr_combined$patient_id)

# now we have 177 patients
# try to reach to db_clin database
# load main database
db.language <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "Meine Ablage", "My Drive")
db.path <- paste0("G:/", db.language , "/sci/dbs/graz/nirs/")

orig_db <- read.delim(
  paste0(db.path, "db_total_cleaned_11022026.csv2"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) |> 
  as_tibble()

orig_db <- orig_db |>
  mutate(
    hemikraniektomie = case_when(
      hemikraniektomie == "Ja, Hemikraniektomie" ~ TRUE,
      hemikraniektomie == "Nein" ~ FALSE,
      TRUE ~ FALSE
    ),
    a_thrombolyse = case_when(
      a_thrombolyse == "Ja" ~ TRUE,
      TRUE ~ FALSE
    ),
    ipsilat_extrakran_hochgr_aci_stenose_verschluss = case_when(
      ipsilat_extrakran_hochgr_aci_stenose_verschluss == "Ja" ~ TRUE,
      TRUE ~ FALSE
    )
  )


# this 3 patients are not in the orig database!
#> intersect(orig_db$patient_id, db_new$patient_id)
#[1] 1001745103 1001906278 1001895711 1000639656
#> setdiff(db_new$patient_id, intersect(orig_db$patient_id, db_new$patient_id))
#[1] 1001887933 1001874454 1000299807


# transfering 4 patients from orig_db to db_clin
library(dplyr)

# ------------------------------------------------------------
# Patient IDs to transfer
# ------------------------------------------------------------

ids_to_add <- c(
  1001745103,
  1001906278,
  1001895711,
  1000639656
)


# ------------------------------------------------------------
# 1. Check that all requested patients exist in orig_db
# ------------------------------------------------------------

missing_in_orig <- setdiff(
  ids_to_add,
  unique(orig_db$patient_id)
)

missing_in_orig
# should return numeric(0)


# ------------------------------------------------------------
# 2. Check that they are not already present in db_clin
# ------------------------------------------------------------

already_in_db_clin <- intersect(
  ids_to_add,
  unique(db_clin$patient_id)
)

already_in_db_clin

db_clin$patient_id <- as.integer(db_clin$patient_id)
db_clin$alter <- as.integer(db_clin$alter)
db_clin$andere_embolektomiesysteme <- as.logical(db_clin$andere_embolektomiesysteme)
db_clin$m1_verschluss_angio <- as.logical(db_clin$m1_verschluss_angio)
db_clin$m2_verschluss_angio <- as.logical(db_clin$m2_verschluss_angio)
db_clin$t_gabelverschluss_angio <- as.logical(db_clin$t_gabelverschluss_angio)
db_clin$distaler_aci_verschluss_angio <- as.logical(db_clin$distaler_aci_verschluss_angio)
db_clin$anderer_verschluss_angio <- as.logical(db_clin$anderer_verschluss_angio)
db_clin$nihss_bei_aufnahme_nur_summe <- as.integer(db_clin$nihss_bei_aufnahme_nur_summe)
db_clin$m1_verschluss <- as.logical(db_clin$m1_verschluss)
db_clin$t_gabelverschluss <- as.logical(db_clin$t_gabelverschluss)
db_clin$hemikraniektomie <- as.logical(db_clin$hemikraniektomie)

# ideally numeric(0)
orig_db$thrombozytenzahl_bei_aufnahme <- as.integer(orig_db$thrombozytenzahl_bei_aufnahme)
orig_db$serumglucose_bei_aufnahme <- as.integer(orig_db$serumglucose_bei_aufnahme)
orig_db$crp <- as.integer(orig_db$crp)
orig_db$haematokrit <- as.integer(orig_db$haematokrit)
orig_db$iv_thrombolyse_mit_rtpa <- as.logical(orig_db$iv_thrombolyse_mit_rtpa)
orig_db$anzahl_verwendeter_devices <- as.integer(orig_db$anzahl_verwendeter_devices)
orig_db$penumbra_system <- as.integer(orig_db$penumbra_system)
orig_db$persistierender_gefaessverschluss <- as.integer(orig_db$persistierender_gefaessverschluss)
orig_db$nur_partielle_rekanalisation <- as.logical(orig_db$nur_partielle_rekanalisation)
orig_db$trevo_stent <- as.character(orig_db$trevo_stent)
orig_db$andere_embolektomiesysteme <- as.logical(orig_db$andere_embolektomiesysteme)
orig_db$hemikraniektomie <- as.logical(orig_db$hemikraniektomie)
str(orig_db$andere_embolektomiesysteme)

db_clin <- db_clin |> dplyr::select(-rr_syst_angiographie_beginn)
orig_db <- orig_db |> dplyr::select(-rr_syst_angiographie_beginn)

db_clin <- db_clin |> dplyr::select(-rr_syst_angiographie_ende)
orig_db <- orig_db |> dplyr::select(-rr_syst_angiographie_ende)

db_clin <- db_clin |> dplyr::select(-rr_diast_angiographie_ende)
orig_db <- orig_db |> dplyr::select(-rr_diast_angiographie_ende)

db_clin <- db_clin |> dplyr::select(-maximaler_rr_syst_praeprozedural)
orig_db <- orig_db |> dplyr::select(-maximaler_rr_syst_praeprozedural)

db_clin <- db_clin |> dplyr::select(-maximaler_rr_diast_praeprozedural)
orig_db <- orig_db |> dplyr::select(-maximaler_rr_diast_praeprozedural)

db_clin <- db_clin |> dplyr::select(-maximaler_rr_diast_intraprozedural)
orig_db <- orig_db |> dplyr::select(-maximaler_rr_diast_intraprozedural)

db_clin <- db_clin |> dplyr::select(-maximaler_rr_syst_intraprozedural)
orig_db <- orig_db |> dplyr::select(-maximaler_rr_syst_intraprozedural)

db_clin <- db_clin |> dplyr::select(-rechts2)
orig_db <- orig_db |> dplyr::select(-rechts2)

db_clin <- db_clin |> dplyr::select(-links2)
orig_db <- orig_db |> dplyr::select(-links2)

db_clin <- db_clin |> dplyr::select(-nihss_aufn_max)
orig_db <- orig_db |> dplyr::select(-nihss_aufn_max)

db_clin <- db_clin |> dplyr::select(-dauer_der_maschinellen_beatmung)
orig_db <- orig_db |> dplyr::select(-dauer_der_maschinellen_beatmung)

db_clin <- db_clin |> dplyr::select(-mrs_bei_entlassung)
orig_db <- orig_db |> dplyr::select(-mrs_bei_entlassung)

db_clin <- db_clin |> dplyr::select(-vorbehandlung_mit_thrombozytenfh)
orig_db <- orig_db |> dplyr::select(-vorbehandlung_mit_thrombozytenfh)

db_clin <- db_clin |> dplyr::select(-vorbehandlung_orale_antikoagulanzien)
orig_db <- orig_db |> dplyr::select(-vorbehandlung_orale_antikoagulanzien)

db_clin <- db_clin |> dplyr::select(-clopidogrel2)
orig_db <- orig_db |> dplyr::select(-clopidogrel2)

db_clin <- db_clin |> dplyr::select(-gpiib_iiia_antagonisten)
orig_db <- orig_db |> dplyr::select(-gpiib_iiia_antagonisten)

db_clin <- db_clin |> dplyr::select(-solitaire_stent)
orig_db <- orig_db |> dplyr::select(-solitaire_stent)

db_clin <- db_clin |> dplyr::select(-trevo_stent)
orig_db <- orig_db |> dplyr::select(-trevo_stent)

patients_to_add <- orig_db |>
  filter(patient_id %in% ids_to_add) |>
  mutate(
    .patient_order = match(patient_id, ids_to_add)
  ) |>
  arrange(.patient_order) |>
  dplyr::select(-.patient_order)

patients_to_add |>
  count(patient_id)

patients_to_add_clin <- patients_to_add |>
  dplyr::select(any_of(names(db_clin)))

db_clin <- bind_rows(
  db_clin,
  patients_to_add_clin
)



# futher with new loaded db_msr and db_clin
# calculate distances and write angigography first try and opening

# add all 7 in the mix
ids_to_add <- c(ids_to_add, c("1001887933", "1001874454", "1000299807"))
library(dplyr)

library(dplyr)
library(stringr)

db_clin <- db_clin |>
  mutate(
    onset_time_dt = as.POSIXct(
      onset_time,
      format = "%d.%m.%Y %H:%M",
      tz = "Europe/Vienna"
    ),
    
    thrombectomy_time_dt = as.POSIXct(
      thrombectomy_time,
      format = "%d.%m.%Y %H:%M",
      tz = "Europe/Vienna"
    ),
    
    time_diff_onset_thrombectomy_time = case_when(
      patient_id %in% ids_to_add &
        !is.na(onset_time_dt) &
        !is.na(thrombectomy_time_dt) ~
        as.numeric(
          difftime(
            thrombectomy_time_dt,
            onset_time_dt,
            units = "mins"
          )
        ),
      TRUE ~ time_diff_onset_thrombectomy_time
    ),
    
    # keep lyse date as plain text
    lyse_date_tmp = as.character(datum_lysebeginn),
    
    # extract only HH:MM or HH:MM:SS from the time field
    lyse_clock_tmp = str_extract(
      as.character(uhrzeit_lysebeginn),
      "\\d{1,2}:\\d{2}(?::\\d{2})?$"
    ),
    
    # add seconds if only HH:MM is present
    lyse_clock_tmp = if_else(
      !is.na(lyse_clock_tmp) &
        str_count(lyse_clock_tmp, ":") == 1,
      paste0(lyse_clock_tmp, ":00"),
      lyse_clock_tmp
    ),
    
    # combine directly: no intermediate timezone/date conversion
    lyse_time = as.POSIXct(
      paste(
        lyse_date_tmp,
        lyse_clock_tmp
      ),
      format = "%d.%m.%Y %H:%M:%S",
      tz = "Europe/Vienna"
    ),
    
    time_diff_onset_lyse = case_when(
      patient_id %in% ids_to_add &
        !is.na(onset_time_dt) &
        !is.na(lyse_time) ~
        as.numeric(
          difftime(
            lyse_time,
            onset_time_dt,
            units = "mins"
          )
        ),
      TRUE ~ time_diff_onset_lyse
    ),
    time_diff_lyse_angio = case_when(
      patient_id %in% ids_to_add &
        !is.na(lyse_time) &
        !is.na(thrombectomy_time_dt) ~
        as.numeric(
          difftime(
            thrombectomy_time_dt,
            lyse_time,
            units = "mins"
          )
        ),
      TRUE ~ time_diff_onset_lyse
    ),
  ) |>
  dplyr::select(
    -onset_time_dt,
    -thrombectomy_time_dt,
    -lyse_date_tmp,
    -lyse_clock_tmp
  )


# display
db_clin |> filter(patient_id %in% ids_to_add) |> 
  dplyr::select(patient_id, stroke_typ, onset_time, thrombectomy_time, 
                lyse_time, time_diff_lyse_angio, first_try_datetime, opening_datetime)

db_clin |> filter(patient_id %in% ids_to_add) |> 
  dplyr::select(patient_id, stroke_typ, onset_time, datum_lysebeginn, 
                uhrzeit_lysebeginn, lyse_time, time_diff_onset_lyse)

db_clin |> filter(patient_id %in% ids_to_add) |> 
  dplyr::select(patient_id, tici_grp, stroke_typ, mt_duration)

# mt duration calculation
mt_duration <- db_new_ready |> group_by(patient_id) |>  
  filter(patient_id %in% ids_to_add) |> 
  dplyr::select(patient_id, elapsed_seconds) |> 
  summarise(mt_dura <- max(elapsed_seconds)/60)
mt_duration

# put mt_duration in the database
idx <- match(db_clin$patient_id, mt_duration$patient_id)

db_clin$mt_duration[!is.na(idx)] <-
  mt_duration[[2]][idx[!is.na(idx)]]

# so to do 24.09. is
# extract angio data from the console for all 7 patients
# put this data into database db_msr
# calculate differences
# run the statistics again

# AFFECTED HEMISPHERE DESIGNATION!!! for new 7

ids_to_fix <- c(
  1001887933,
  1001874454,
  1000299807,
  1001745103,
  1001906278,
  1001895711,
  1000639656
)

# ------------------------------------------------------------
# 1. Correct datetime for two patients
# 2. Assign Affected / Unaffected from db_clin$side
# ------------------------------------------------------------

# first get stroke side for the 7 patients
side_lookup <- db_clin |>
  filter(patient_id %in% ids_to_fix) |>
  dplyr::select(patient_id, side) |>
  distinct(patient_id, .keep_all = TRUE)

db_msr_combined <- db_msr_combined |>
  
  # temporarily parse datetime
  mutate(
    datetime_tmp = as.POSIXct(
      datetime,
      format = "%d.%m.%Y %H:%M:%S",
      tz = "Europe/Vienna"
    ),
    
    # shift exactly +/- 1 hour
    datetime_tmp = case_when(
      patient_id == 1001906278 ~ datetime_tmp + 3600,
      patient_id == 1000299807 ~ datetime_tmp - 3600,
      TRUE ~ datetime_tmp
    )
  ) |>
  
  # add clinical stroke side
  left_join(
    side_lookup,
    by = "patient_id"
  ) |>
  
  mutate(
    hemisphere_status = case_when(
      
      # only modify the 7 new patients
      patient_id %in% ids_to_fix &
        !is.na(side) &
        !is.na(hemisphere_side) &
        tolower(side) == tolower(hemisphere_side) ~
        "Affected",
      
      patient_id %in% ids_to_fix &
        !is.na(side) &
        !is.na(hemisphere_side) &
        tolower(side) != tolower(hemisphere_side) ~
        "Unaffected",
      
      # keep existing values for everyone else
      TRUE ~ hemisphere_status
    ),
    
    # return datetime to original character format
    datetime = format(
      datetime_tmp,
      "%d.%m.%Y %H:%M:%S"
    )
  ) |>
  
  dplyr::select(
    -datetime_tmp,
    -side
  )

library(dplyr)
library(lubridate)

ids_to_fix <- c(
  1001887933,
  1001874454,
  1000299807,
  1001745103,
  1001906278,
  1001895711,
  1000639656
)

# ------------------------------------------------------------
# Clinical event times
# ------------------------------------------------------------

event_lookup <- db_clin |>
  filter(patient_id %in% ids_to_fix) |>
  dplyr::select(
    patient_id,
    first_try_datetime,
    opening_datetime
  ) |>
  distinct(patient_id, .keep_all = TRUE) |>
  mutate(
    first_try_minute = floor_date(
      as.POSIXct(
        first_try_datetime,
        format = "%d.%m.%Y %H:%M",
        tz = "Europe/Vienna"
      ),
      unit = "minute"
    ),
    
    opening_minute = floor_date(
      as.POSIXct(
        opening_datetime,
        format = "%d.%m.%Y %H:%M",
        tz = "Europe/Vienna"
      ),
      unit = "minute"
    )
  ) |>
  dplyr::select(
    patient_id,
    first_try_minute,
    opening_minute
  )

db_msr_combined <- db_msr_combined |>
  mutate(
    datetime_dt = as.POSIXct(
      datetime,
      format = "%d.%m.%Y %H:%M:%S",
      tz = "Europe/Vienna"
    ),
    
    datetime_minute = floor_date(
      datetime_dt,
      unit = "minute"
    )
  ) |>
  
  left_join(
    event_lookup,
    by = "patient_id"
  ) |>
  
  group_by(patient_id) |>
  
  mutate(
    # --------------------------------------------------------
    # Find first actual NIRS timestamp within first-try minute
    # --------------------------------------------------------
    first_try_first_time = if (
      any(
        !is.na(first_try_minute) &
        datetime_minute == first_try_minute
      )
    ) {
      min(
        datetime_dt[
          !is.na(first_try_minute) &
            datetime_minute == first_try_minute
        ],
        na.rm = TRUE
      )
    } else {
      as.POSIXct(NA)
    },
    
    # --------------------------------------------------------
    # Find first actual NIRS timestamp within opening minute
    # --------------------------------------------------------
    opening_first_time = if (
      any(
        !is.na(opening_minute) &
        datetime_minute == opening_minute
      )
    ) {
      min(
        datetime_dt[
          !is.na(opening_minute) &
            datetime_minute == opening_minute
        ],
        na.rm = TRUE
      )
    } else {
      as.POSIXct(NA)
    },
    
    # --------------------------------------------------------
    # Flag BOTH hemispheres at that timestamp
    # --------------------------------------------------------
    first_try_point = case_when(
      patient_id %in% ids_to_fix &
        !is.na(first_try_first_time) &
        datetime_dt == first_try_first_time ~ TRUE,
      
      patient_id %in% ids_to_fix ~ FALSE,
      
      TRUE ~ first_try_point
    ),
    
    recanalisation_point = case_when(
      patient_id %in% ids_to_fix &
        !is.na(opening_first_time) &
        datetime_dt == opening_first_time ~ TRUE,
      
      patient_id %in% ids_to_fix ~ FALSE,
      
      TRUE ~ recanalisation_point
    )
  ) |>
  
  ungroup() |>
  
  dplyr::select(
    -datetime_dt,
    -datetime_minute,
    -first_try_minute,
    -opening_minute,
    -first_try_first_time,
    -opening_first_time
  )

db_msr_combined |>
  filter(
    patient_id %in% ids_to_fix,
    first_try_point | recanalisation_point
  ) |>
  dplyr::select(
    patient_id,
    datetime,
    hemisphere_side,
    hemisphere_status,
    first_try_point,
    recanalisation_point
  ) |>
  arrange(
    patient_id,
    datetime,
    hemisphere_side
  )

ids_to_fix <- c(
  1001887933,
  1001874454,
  1000299807,
  1001745103,
  1001906278,
  1001895711,
  1000639656
)

db_msr <- db_msr |>
  dplyr::mutate(
    first_try_point = dplyr::case_when(
      patient_id %in% ids_to_fix &
        hemisphere_status == "Unaffected" ~ FALSE,
      TRUE ~ first_try_point
    ),
    
    recanalisation_point = dplyr::case_when(
      patient_id %in% ids_to_fix &
        hemisphere_status == "Unaffected" ~ FALSE,
      TRUE ~ recanalisation_point
    )
  )

# adding missing variables
db_clin_old <- read.csv2(paste0(db.path_local, "clinical_170_clean_v1.csv"))
db_clin_old <- read.csv2(paste0(db.path_local, "db_total_cleaned_11022026.csv2"))

db_clin_old <- read.csv2(
"db_total_cleaned_11022026.csv2",
header = TRUE,
sep = "\t",
dec = ",",
stringsAsFactors = FALSE,
check.names = FALSE
)

vars_to_add <- c(
  "vorbehandlung_mit_thrombozytenfh",
  "vorbehandlung_orale_antikoagulanzien",
  "nihss_entl_max"
)

lookup_old <- db_clin_old |>
  dplyr::select(
    patient_id,
    dplyr::all_of(vars_to_add)
  ) |>
  dplyr::distinct(patient_id, .keep_all = TRUE)

db_clin <- db_clin |>
  dplyr::left_join(
    lookup_old,
    by = "patient_id"
  )

db_clin <- db_clin |>
  dplyr::mutate(
    
    # --------------------------------------------------------
    # Time of admission
    # --------------------------------------------------------
    zeitpunkt = dplyr::case_when(
      zeitpunkt == "Standard duty hours (Moâ€“Fr 7:00â€“15:00)" ~
        "Standard duty hours (Mo–Fr 7:00–15:00)",
      TRUE ~ zeitpunkt
    ),
    
    # --------------------------------------------------------
    # Final infarct
    # --------------------------------------------------------
    final_infarct = dplyr::case_when(
      final_infarct %in% c("0â€“33%", "0–33%", "0-33%") ~ "0–33%",
      final_infarct %in% c("34â€“66%", "34–66%", "34-66%") ~ "34–66%",
      final_infarct %in% c("67â€“100%", "67–100%", "67-100%") ~ "67–100%",
      TRUE ~ final_infarct
    ),
    
    # --------------------------------------------------------
    # Stenosis in treated vessel
    # --------------------------------------------------------
    stenose_im_therapierten_gefaess = dplyr::case_when(
      stenose_im_therapierten_gefaess ==
        "Mildâ€“moderate stenosis" ~
        "Mild–moderate stenosis",
      TRUE ~ stenose_im_therapierten_gefaess
    ),
    
    # --------------------------------------------------------
    # 3-month mRS
    # --------------------------------------------------------
    mRS_3mo_grp = dplyr::case_when(
      mRS_3mo_grp %in% c("0â€“2", "0–2", "0-2") ~ "0–2",
      mRS_3mo_grp %in% c("3â€“6", "3–6", "3-6") ~ "3–6",
      TRUE ~ mRS_3mo_grp
    )
  )

db_clin |>
  dplyr::select(
    zeitpunkt,
    final_infarct,
    stenose_im_therapierten_gefaess,
    mRS_3mo_grp
  ) |>
  lapply(unique)

db_clin <- db_clin |>
  dplyr::mutate(
    stenose_im_therapierten_gefaess = dplyr::case_when(
      stenose_im_therapierten_gefaess == "0 = Keine Stenose" ~ "No stenosis",
      TRUE ~ stenose_im_therapierten_gefaess
    )
  )

db_clin$vessel_type_preangio[,db$clin$patient_id==1001887933 | patient_id == 1000299807] <- "M1"

idx <- match(db_clin$patient_id, c("1001887933", "1000299807"))

db_clin$vessel_type_preangio[!is.na(idx)] <- "M1"
  
idx <- match(db_clin$patient_id, c("1001887933"))
db_clin$final_infarct[!is.na(idx)] <- "0–33%"

db_clin <- db_clin |> mutate(
  geschlecht = case_when(
    geschlecht == "männlich" ~ "men",
    geschlecht == "weiblich" ~ "women",
    TRUE~geschlecht
  )
)

db_clin <- db_clin |> mutate(
  infarktfruehzeichen = case_when(
    infarktfruehzeichen == "Ja, aber  < 1/3 ACM"~ "Less than 1/3 of MCA",
    infarktfruehzeichen == "Nein" ~ "No",
    TRUE ~ infarktfruehzeichen
  )  
)

table(db_clin$stenose_im_therapierten_gefaess)

db_clin <- db_clin |> mutate(
  stenose_im_therapierten_gefaess = case_when(
    stenose_im_therapierten_gefaess=="" ~ "Unknown",
    TRUE ~ stenose_im_therapierten_gefaess
  )
)

db_clin <- db_clin |>
  dplyr::mutate(
    schlaganfall_aetiologie_toast = dplyr::case_when(
      
      schlaganfall_aetiologie_toast %in% c(
        "CE - Cardioembolic",
        "Kardioembolisch"
      ) ~ "CE - Cardioembolic",
      
      schlaganfall_aetiologie_toast ==
        "LAA - Atherothrombotic" ~ "LAA - Atherothrombotic",
      
      schlaganfall_aetiologie_toast ==
        "Other causes" ~ "Other causes",
      
      schlaganfall_aetiologie_toast %in% c(
        "Unbek. Ursache (Diagnostik nicht komplett)",
        "Unclassified",
        "unknown",
        "Unknown etiology"
      ) ~ "Unknown etiology",
      
      TRUE ~ schlaganfall_aetiologie_toast
    )
  )

table(db_t$schlaganfall_aetiologie_toast)

table(db_t$stenose_im_therapierten_gefaess)

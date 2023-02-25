# Name: ADVS
#
# Label: Vital Signs Analysis Dataset
#
# Input: adsl, vs
library(admiral)# Contains example datasets from the CDISC pilot project
library(tidyverse)
library(lubridate)
library(stringr)
library(purrr)
library(xportr)
library(metacore)
library(metatools)
# Load source datasets ----

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
adsl <- haven::read_xpt("adam/adsl.xpt")
# source('programs/adsl.R')
# files_ls <- ls()[which(!str_detect(ls(),'adsl'))]
# rm(list = files_ls)
# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

sdtms <- c('dm','ds','ex','ae','lb','suppdm','suppds','sv','vs','qs','sc','suppae')

read_sdtm <- \(x){
  y <- haven::read_xpt(paste0("sdtm/",x,".xpt")) %>% convert_blanks_to_na()
  assign(x,y,envir = .GlobalEnv)
}

walk(sdtms, ~read_sdtm(.x))

# Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~VSTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "SYSBP", "SYSBP", "Systolic Blood Pressure (mmHg)", 1,
  "DIABP", "DIABP", "Diastolic Blood Pressure (mmHg)", 2,
  "PULSE", "PULSE", "Pulse Rate (beats/min)", 3,
  "WEIGHT", "WEIGHT", "Weight (kg)", 4,
  "HEIGHT", "HEIGHT", "Height (cm)", 5,
  "TEMP", "TEMP", "Temperature (C)", 6
)
attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"


# Assign ANRLO/HI, A1LO/HI
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "SYSBP", 90, 130, 70, 140,
  "DIABP", 60, 80, 40, 90,
  "PULSE", 60, 100, 40, 110,
  "TEMP", 36.5, 37.5, 35, 38
)
# ASSIGN AVALCAT1
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "HEIGHT", 1, ">100 cm",
  "HEIGHT", 2, "<= 100 cm"
)

# User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate()`.
format_avalcat1n <- function(param, aval) {
  case_when(
    param == "HEIGHT" & aval > 140 ~ 1,
    param == "HEIGHT" & aval <= 140 ~ 2
  )
}

format_racen <- function(x) {
  case_when(
    x == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 6,
    x == 'ASIAN' ~ 3,
    x == 'BLACK OR AFRICAN AMERICAN' ~ 2,
    x == 'WHITE' ~ 1,
    TRUE ~ NA_real_
  )
}

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P,TRT01AN, TRT01PN,SITEID,AGE, AGEGR1N, RACE, RACEN, SEX, SAFFL, AGEGR1)

advs1 <- vs %>%
  # Join ADSL with VS (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
derive_vars_dt(
  new_vars_prefix = "A",
  dtc = VSDTC
) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

advs2 <- advs1 %>%
  ## Add PARAMCD only - add PARAM etc later ----
derive_vars_merged_lookup(
  dataset_add = param_lookup,
  new_vars = vars(PARAMCD),
  by_vars = vars(VSTESTCD)
) %>%
  ## Calculate AVAL and AVALC ----
mutate(
  AVAL = VSSTRESN,
  AVALC = VSSTRESC
)


## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
advs3 <- advs2 %>%
  # Derive Timing
  mutate(
    ATPTN = VSTPTNUM,
    ATPT = VSTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(case_when(
      VISIT == "BASELINE" ~ "0",
      str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", "")),
      TRUE ~ NA_character_
    ))
  )

## Derive a new record as a summary record (e.g. mean of the triplicates at each time point) ----
# advs <- advs %>%
#   derive_summary_records(
#     by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT, ADY),
#     filter = !is.na(AVAL),
#     analysis_var = AVAL,
#     summary_fun = mean,
#     set_values_to = vars(DTYPE = "AVERAGE")
#   )

advs4 <- advs3 %>%
  ## Calculate ONTRTFL ----
derive_var_ontrtfl(
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  filter_pre_timepoint = AVISIT == "Baseline"
)

## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
# Also accommodates the ranges A1LO, A1HI
advs5 <- advs4 %>%
  derive_vars_merged(dataset_add = range_lookup, by_vars = vars(PARAMCD)) %>%
  # Calculate ANRIND
  derive_var_anrind()

## Derive baseline flags ----
advs6 <- advs5 %>% mutate(BASETYPE=ATPT) %>%
  # Calculate BASETYPE
  # derive_var_basetype(
  #   basetypes = rlang::exprs(
  #     "AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
  #     "AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
  #     "AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
  #     "LAST" = is.na(ATPTN) & !is.na(AVISIT)
  #   )
  # ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPT),
      order = vars(ADT, VISITNUM, VSSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) &
                #ADT <= TRTSDT & !is.na(BASETYPE) & is.na(AVISIT))
                ADT <= TRTSDT & !is.na(AVISIT))
  )

## Derive baseline information ----
advs7 <- advs6 %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()


## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
advs8 <- advs7 %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = ANL01FL,
      by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT),
      order = vars(ADT, AVAL),
      mode = "last"
    ),
    filter = !is.na(AVISITN) # & AVISITN>0 & ONTRTFL == "Y"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
advs9 <- advs8 %>%
  # Assign TRTA, TRTP
  # Create End of Treatment Record
  # restrict_derivation(
  #   derivation = derive_var_extreme_flag,
  #   args = params(
  #     by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
  #     order = vars(ADT),
  #     new_var = EOTFL,
  #     mode = "last"
  #   ),
  #   filter = (4 < VISITNUM &
#     VISITNUM <= 13 & ANL01FL == "Y")
# ) %>%
# filter(EOTFL == "Y") %>%
# mutate(
#   AVISIT = "End of Treatment",
#   AVISITN = 99
# ) %>%
# union_all(advs) %>%
# select(-EOTFL) %>%
mutate(
  TRTP = TRT01P,
  TRTA = TRT01A
)

## Get ASEQ and AVALCATx and add PARAM/PARAMN ----
advs10 <- advs9 %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN),
    check_type = "error"
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  derive_vars_merged(dataset_add = avalcat_lookup, by_vars = vars(PARAMCD, AVALCA1N)) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -VSTESTCD), by_vars = vars(PARAMCD))


# Add all ADSL variables
advs11 <- advs10 %>%
  derive_vars_merged(
    dataset_add = select(adsl, STUDYID, USUBJID, TRT01AN, TRT01PN),
    by_vars = vars(STUDYID, USUBJID, TRT01AN, TRT01PN)
  ) %>%
  mutate(TRTPN=TRT01PN,
         TRTAN=TRT01AN
         #PCHG= janitor::round_half_up(PCHG,1),
         # CHG=janitor::round_half_up(CHG,1),
         # BASETYPE=ifelse(is.na(ATPT), NA_character_, BASETYPE)
         )

advs <- advs11 %>%
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = vars(ADT, AVISITN, ATPTN, AVAL),
    mode = "last",
    filter = (is.na(ABLFL) & between(AVISITN,4,26) & !is.na(AVISIT)),
    set_values_to = vars(
      AVISIT = "End of Treatment",
      AVISITN = 99
    )
  ) %>%
  mutate(RACEN=format_racen(RACE))

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

advs_spec <- metacore %>% select_dataset('ADVS')

advs <- advs %>% drop_unspec_vars(advs_spec)

advs <- advs %>%
  check_variables(advs_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(advs_spec) %>% # Orders the columns according to the spec
  sort_by_key(advs_spec) %>%
  xportr_label(advs_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(advs_spec) %>%
  xportr_format(advs_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADVS")

# Save output ----
xportr_write(advs, "adam/advs.xpt")

# dir <- tempdir() # Change to whichever directory you want to save the dataset in
# saveRDS(advs, file = file.path(dir, "advs.rds"), compress = "bzip2")


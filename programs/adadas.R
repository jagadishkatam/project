# Name: ADADAS
#
# Label: ADAS-Cog Analysis Dataset
#
# Description: Based on CDISC Pilot data, create ADADAS analysis dataset
#
# Input: adsl, eg
library(admiral)# Contains example datasets from the CDISC pilot project
library(tidyverse)
library(lubridate)
library(stringr)
library(purrr)
library(xportr)
library(metacore)
library(metatools)

# Load source datasets ----
adsl <- haven::read_xpt("adam/adsl.xpt")
# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

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
param_lookup <- qs %>% filter(QSCAT=="ALZHEIMER'S DISEASE ASSESSMENT SCALE") %>%
  select(QSTESTCD, PARAM=QSTEST) %>% unique() %>%
  mutate(PARAMCD=QSTESTCD, PARAMN=row_number())

format_AWRANGE <- \(x){
  case_when(
    x=='BASELINE' ~ '<=1',
    x=='WEEK 8' ~ '2-84',
    x=='WEEK 16' ~ '85-140',
    x=='WEEK 24' ~ '>140',
    TRUE ~ NA_character_
  )
}

format_AWTARGET <- \(x){
  case_when(
    x=='BASELINE' ~ 1,
    x=='WEEK 8' ~ 56,
    x=='WEEK 16' ~ 112,
    x=='WEEK 24' ~ 168,
    TRUE ~ NA_real_
  )
}


format_AWLO <- \(x){
  case_when(
    # x=='BASELINE' ~ NA,
    x=='WEEK 8' ~ 2,
    x=='WEEK 16' ~ 85,
    x=='WEEK 24' ~ 140,
    TRUE ~ NA_real_
  )
}


format_AWHI <- \(x){
  case_when(
    x=='BASELINE' ~ 1,
    x=='WEEK 8' ~ 84,
    x=='WEEK 16' ~ 140,
    # x=='WEEK 24' ~ NA,
    TRUE ~ NA_real_
  )
}

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P)

adas <- qs %>% filter(QSCAT=="ALZHEIMER'S DISEASE ASSESSMENT SCALE") %>%
  # Join ADSL & EG (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADTM, ADY ----
derive_vars_dtm(
  new_vars_prefix = "A",
  dtc = QSDTC,
) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADTM))

adas <- adas %>%
  ## Add PARAMCD only (add PARAM, etc later) ----
derive_vars_merged_lookup(
  dataset_add = param_lookup,
  new_vars = vars(PARAMCD, PARAM, PARAMN),
  by_vars = vars(QSTESTCD)
) %>%
  ## Calculate AVAL and AVALC ----
mutate(
  AVAL = QSSTRESN,
  AVALC = QSSTRESC
)

## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adas <- adas %>%
  # Derive Timing
  mutate(
    ADT = date(ADTM),
    AVISIT = case_when(
      ADY <=1 ~ 'BASELINE',
      between(ADY,2,84) ~ 'WEEK 8',
      between(ADY,85,140) ~ 'WEEK 16',
      ADY > 140 ~ 'WEEK 24',
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(
      case_when(
        ADY <=1 ~ 0,
        between(ADY,2,84) ~ 8,
        between(ADY,85,140) ~ 16,
        ADY > 140 ~ 24,
        TRUE ~ NA_real_
      )
    ),
  )

vis <- adas %>% select(AVISIT, AVISITN) %>% unique()

visall <- merge(param_lookup %>% filter(QSTESTCD=='ACTOT'),vis, all = T)

adas <- derive_locf_records(
  data = adas,
  dataset_expected_obs = visall,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, PARAM, PARAMN),
  order = vars(AVISITN, AVISIT)
)

## Derive baseline flags ----
adas <- adas %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, PARAMCD),
      order = vars(ADT, VISITNUM, QSSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) | !is.na(AVALC)) & (QSBLFL=='Y')
  )

## Derive baseline information ----
adas <- adas %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()

## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
adas <- adas %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, PARAMCD, AVISIT),
      order = vars(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN)
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
adas <- adas %>%
  # Assign TRTA, TRTP
  mutate(TRTP = TRT01P, TRTA = TRT01A)

## Get ASEQ and AVALCAT1/CHGCAT1 and add PARAM/PARAMN ----
adas <- adas %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM),
    check_type = "error"
  )

# Add all ADSL variables
adadas <- adas %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  ) %>% mutate(AWRANGE=format_AWRANGE(AVISIT),
               AWTARGET=format_AWTARGET(AVISIT),
               AWLO=format_AWLO(AVISIT),
               AWHI=format_AWHI(AVISIT),
               AWTDIFF=ADY-AWTARGET,
               TRTPN=TRT01PN,
               AWU='DAYS')



metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adas_spec <- metacore %>% select_dataset('ADADAS')

adadas <- adadas %>% drop_unspec_vars(adas_spec)

adadas <- adadas %>%
  check_variables(adas_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adas_spec) %>% # Orders the columns according to the spec
  sort_by_key(adas_spec) %>%
  # xportr_type(adae_spec) %>%
  xportr_label(adas_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(adas_spec) %>%
  xportr_format(adas_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADADAS")

# Save output ----
xportr_write(adadas, "adam/adadas.xpt")



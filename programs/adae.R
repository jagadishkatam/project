# Name: ADAE
#
# Label: Adverse Event Analysis Dataset
#
# Input: ae, adsl, ex_single
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
# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

sdtms <- c('dm','ds','ex','ae','lb','suppdm','suppds','sv','vs','qs','sc','suppae')

read_sdtm <- \(x){
  y <- haven::read_xpt(paste0("sdtm/",x,".xpt")) %>% convert_blanks_to_na()
  assign(x,y,envir = .GlobalEnv)
}

walk(sdtms, ~read_sdtm(.x))

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, RACEN)

adae <- ae %>%
  # join adsl to ae
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = vars(STUDYID, USUBJID)
  ) %>%
  ## Derive analysis start time ----
  derive_vars_dt(
    new_vars_prefix = "AST",
    dtc = AESTDTC,
    highest_imputation = "D"
  ) %>%
  ## Derive analysis end time ----
derive_vars_dt(
  new_vars_prefix = "AEN",
  dtc = AEENDTC
 ) %>%
  ## Derive analysis end/start date ----
  #derive_vars_dtm_to_dt(vars(ASTDTM, AENDTM)) %>%
  ## Derive analysis start relative day and  analysis end relative day ----
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ASTDT, AENDT)
  ) %>%
  ## Derive analysis duration (value and unit) ----
  derive_vars_duration(
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT,
    in_unit = "days",
    out_unit = "days",
    add_one = TRUE,
    trunc_out = FALSE
  )

ex_ext <- derive_vars_dtm(
  ex,
  dtc = EXSTDTC,
  new_vars_prefix = "EXST",
  flag_imputation = "none"
)

adae <- adae %>%
  ## Derive last dose date/time ----
  # derive_var_last_dose_date(
  #   ex_ext,
  #   filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
  #     !is.na(EXSTDTM),
  #   dose_date = EXSTDTM,
  #   analysis_date = ASTDT,
  #   new_var = LDOSEDTM,
  #   single_dose_condition = (EXSTDTC == EXENDTC),
  #   output_datetime = TRUE
  # ) %>%
  ## Derive severity / causality / ... ----
  mutate(
    ASEV = AESEV,
    AREL = AEREL
  ) %>%
  mutate(TRTEMFL=ifelse(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT>=TRTSDT, 'Y', NA_character_)) %>%
  ## Derive treatment emergent flag ----
  # derive_var_trtemfl(
  #   start_date = ASTDT,
  #   end_date = AENDT,
  #   trt_start_date = TRTSDT,
  #   trt_end_date = TRTEDT,
  #   end_window =
  # ) %>%
  ## Derive occurrence flags: first occurrence of most severe AE ----
  # create numeric value ASEVN for severity
  mutate(
    ASEVN = as.integer(factor(ASEV, levels = c("MILD", "MODERATE", "SEVERE", "DEATH THREATENING")))
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID),
      order = vars(desc(ASEVN), ASTDT, AESEQ),
      new_var = AOCCIFL,
      mode = "first"
    ),
    filter = TRTEMFL == "Y"
  )

cq1 <- paste(c('APPLICATION', 'DERMATITIS', 'ERYTHEMA', 'BLISTER'),collapse = '|')
cq2 <- paste(c('COLD SWEAT', 'HYPERHIDROSIS', 'ALOPECIA'), collapse = '|')

# Join all ADSL with AE
adae <- adae %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  mutate(CQ01NAM=ifelse(str_detect(AEDECOD,cq1) |
                          (AEBODSYS=='SKIN AND SUBCUTANEOUS TISSUE DISORDERS' & !(str_detect(AEDECOD,cq2))),
                        'DERMATOLOGIC EVENTS',
                        NA_character_))

AOCC01FL <- adae %>% filter(!is.na(CQ01NAM) & TRTEMFL=='Y') %>% arrange(USUBJID, ASTDT, AESEQ) %>% group_by(USUBJID) %>%
  slice_head(n=1) %>%
  mutate(AOCC01FL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCC01FL) %>% ungroup()

AOCC02FL <- adae %>% filter(AESER=='Y' & TRTEMFL=='Y') %>% arrange(USUBJID, ASTDT, AESEQ) %>% group_by(USUBJID) %>%
  slice_head(n=1) %>%
  mutate(AOCC02FL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCC02FL) %>% ungroup()

AOCC03FL <- adae %>% filter(AESER=='Y' & TRTEMFL=='Y') %>% arrange(USUBJID, AEBODSYS, ASTDT, AESEQ) %>% group_by(USUBJID,AEBODSYS) %>%
  slice_head(n=1) %>% ungroup() %>%
  mutate(AOCC03FL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCC03FL)

AOCC04FL <- adae %>% filter(AESER=='Y' & TRTEMFL=='Y') %>% arrange(USUBJID, AEBODSYS,AEDECOD, ASTDT, AESEQ) %>% group_by(USUBJID,AEBODSYS, AEDECOD) %>%
  slice_head(n=1) %>% ungroup() %>%
  mutate(AOCC04FL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCC04FL)


AOCCFL <- adae %>% filter(TRTEMFL=='Y') %>% arrange(USUBJID, ASTDT, AESEQ) %>% group_by(USUBJID) %>%
  slice_head(n=1) %>% ungroup() %>%
  mutate(AOCCFL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCCFL)

AOCCSFL <- adae %>% filter(TRTEMFL=='Y') %>% arrange(USUBJID, AEBODSYS, ASTDT, AESEQ) %>%
  group_by(USUBJID,AEBODSYS) %>%
  slice_head(n=1) %>% ungroup() %>%
  mutate(AOCCSFL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCCSFL)

AOCCPFL <- adae %>% filter(TRTEMFL=='Y') %>% arrange(USUBJID, AEBODSYS,AEDECOD, ASTDT, AESEQ) %>% group_by(USUBJID,AEBODSYS, AEDECOD) %>%
  slice_head(n=1) %>% ungroup() %>%
  mutate(AOCCPFL='Y') %>% select(USUBJID, ASTDT, AESEQ, AOCCPFL)

adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCC01FL, USUBJID, ASTDT, AESEQ, AOCC01FL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
)

adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCC02FL, USUBJID, ASTDT, AESEQ, AOCC02FL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
)

adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCC03FL, USUBJID, ASTDT, AESEQ, AOCC03FL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
)

adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCC04FL, USUBJID, ASTDT, AESEQ, AOCC04FL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
)

adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCCFL, USUBJID, ASTDT, AESEQ, AOCCFL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
)


adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCCSFL, USUBJID, ASTDT, AESEQ, AOCCSFL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
)

format_racen <- function(x) {
  case_when(
    x == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 6,
    x == 'ASIAN' ~ 3,
    x == 'BLACK OR AFRICAN AMERICAN' ~ 2,
    x == 'WHITE' ~ 1,
    TRUE ~ NA_real_
  )
}

adae <- derive_vars_merged(
  adae,
  dataset_add = select(AOCCPFL, USUBJID, ASTDT, AESEQ, AOCCPFL),
  by_vars = vars( USUBJID, ASTDT, AESEQ)
) %>%
  mutate(TRTA=TRT01A, TRTAN=TRT01AN, RACEN=format_racen(RACE),
         ADURU=ifelse(!is.na(ADURN) & is.na(ASTDTF), 'DAY', NA_character_),
         ADURN=ifelse(!is.na(ADURN) & is.na(ASTDTF), ADURN, NA_real_)
         )

metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adae_spec <- metacore %>% select_dataset('ADAE')

adae <- adae %>% drop_unspec_vars(adae_spec)

adae <- adae %>%
  check_variables(adae_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adae_spec) %>% # Orders the columns according to the spec
  sort_by_key(adae_spec) %>%
  # xportr_type(adae_spec) %>%
  xportr_label(adae_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(adae_spec) %>%
  xportr_format(adae_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADAE")

# Save output ----
xportr_write(adae, "adam/adae.xpt")

# dir <- tempdir() # Change to whichever directory you want to save the dataset in
# saveRDS(adae, file = file.path(dir, "adae.rds"), compress = "bzip2")


# Name: ADLBH
#
# Label: Analysis Dataset Lab Hematology
#
# Input: adsl, lb
library(admiral)# Contains example datasets from the CDISC pilot project
library(tidyverse)
library(lubridate)
library(stringr)
library(purrr)
library(xportr)
library(metacore)
library(metatools)
library(data.table)
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

sdtms <- c('dm','ds','ex','ae','lb','suppdm','suppds','sv','vs','qs','sc','supplb','suppae')

read_sdtm <- \(x){
  y <- haven::read_xpt(paste0("sdtm/",x,".xpt")) %>% convert_blanks_to_na()
  assign(x,y,envir = .GlobalEnv)
}

walk(sdtms, ~read_sdtm(.x))

# Lookup tables ----
metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adlbh_spec <- metacore %>% select_dataset('ADLBH')

supplb <- supplb %>% mutate(IDVARVAL=trimws(IDVARVAL))
lb1 <- metatools::combine_supp(lb, supplb) %>% as_tibble()

# lb1 <- create_var_from_codelist(lb, adlbh_spec, LBTESTCD, Term, decode_to_code = F)

paramcd <- adlbh_spec$codelist %>% as_tibble() %>% unnest(everything()) %>% filter(code_id=='PARAMCD_ADLBH') %>%
  transmute(LBTESTCD=code, PARAMCD=code, PARAM=decode)

paramn <- adlbh_spec$codelist %>% as_tibble() %>% unnest(everything()) %>% filter(code_id=='PARAMN_ADLBH') %>%
  transmute(PARAMN=as.numeric(code), PARAM=decode)

param_lookup <- paramcd %>%
  derive_vars_merged(
    dataset_add = paramn,
    by_vars = vars(PARAM)
  )


format_aind <- function(x,y,z) {
  case_when(
    x < y ~ 'L',
    between(x, y, z) ~ 'N',
    x > z ~ 'H',
    TRUE ~ 'N'
  )
}

format_bind <- function(x,y,z) {
  case_when(
    x < y ~ 'L',
    between(x, y, z) ~ 'N',
    x > z ~ 'H'
  )
}

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(SUBJID, TRTSDT, TRTEDT, TRT01A, TRT01P,TRT01AN, TRT01PN, SITEID,AGE, AGEGR1N, RACE, RACEN, SEX, SAFFL, AGEGR1, COMP24FL, DSRAEFL)

adlbh1 <- lb1 %>% filter(LBCAT=='HEMATOLOGY' & !(VISITNUM %in% c(6,201))) %>%
  # Join ADSL with VS (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
derive_vars_dt(
  new_vars_prefix = "A",
  dtc = LBDTC
) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlbh2 <- adlbh1 %>%
  ## Add PARAMCD only - add PARAM etc later ----
derive_vars_merged_lookup(
  dataset_add = param_lookup,
  new_vars = vars(PARAMCD, PARAM, PARAMN),
  by_vars = vars(LBTESTCD)
) %>%
  ## Calculate AVAL and AVALC ----
mutate(
  AVAL = LBSTRESN,
  AVALC = LBSTRESC,
  A1HI = LBSTNRHI,
  A1LO = LBSTNRLO,
  ANRIND = LBNRIND,
  ANRHI = (1.5*A1HI),
  ANRLO = (0.5*A1LO),
  PARCAT1 = 'HEM',
  ANRIND = format_aind(AVAL, ANRLO, ANRHI),
  ANRIND = ifelse(is.na(ANRLO) & is.na(ANRHI) & LBNRIND=='ABNORMAL', 'H', ANRIND)
)


## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adlbh3 <- adlbh2 %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      str_detect(VISIT, "SCREENING") ~ 'Baseline',
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(case_when(
      VISIT == "SCREENING 1" ~ "0",
      str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", "")),
      TRUE ~ NA_character_
    ))
  ) # %>%
# derive_var_anrind()
## Derive a new record as a summary record (e.g. mean of the triplicates at each time point) ----
# advs <- advs %>%
#   derive_summary_records(
#     by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT, ADY),
#     filter = !is.na(AVAL),
#     analysis_var = AVAL,
#     summary_fun = mean,
#     set_values_to = vars(DTYPE = "AVERAGE")
#   )

# adlbc4 <- adlbc3 %>%
#   ## Calculate ONTRTFL ----
#   derive_var_ontrtfl(
#     start_date = ADT,
#     ref_start_date = TRTSDT,
#     ref_end_date = TRTEDT,
#     filter_pre_timepoint = AVISIT == "Baseline"
#   )

## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
# Also accommodates the ranges A1LO, A1HI
# adlbc4 <- adlbc3 %>%
#   derive_vars_merged(dataset_add = range_lookup, by_vars = vars(PARAMCD)) %>%
#   # Calculate ANRIND
#   derive_var_anrind()

## Derive baseline flags ----
adlbh4 <- adlbh3 %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, PARAMCD),
      order = vars(ADT, VISITNUM, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (VISITNUM==1)
  )

## Derive baseline information ----
adlbh5 <- adlbh4 %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = A1LO,
    new_var = BA1LO
  ) %>%
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = A1HI,
    new_var = BA1HI
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg() %>%
  mutate(BNRIND=format_bind(BASE, ANRLO, ANRHI),
         BNRIND = ifelse(BNRIND=='N' & LBNRIND=='ABNORMAL', 'H', BNRIND)
  )


adlbh6 <- adlbh5 %>%
  derive_var_analysis_ratio(numer_var = BASE,
                            denom_var = A1LO,
                            new_var = BR2A1LO
  ) %>%
  derive_var_analysis_ratio(
    numer_var = BASE,
    denom_var = A1HI,
    new_var = BR2A1HI
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = A1LO,
    new_var = R2A1LO
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = A1HI,
    new_var = R2A1HI
  ) # %>%
  # mutate(BR2A1LO=round(BR2A1LO,6),
  #        BR2A1HI=round(BR2A1HI,2),
  #        R2A1LO=round(R2A1LO,6),
  #        R2A1HI=round(R2A1HI,7)
  # )

# adlbcx <- adlbc7 %>% select(!!!negate_vars(adsl_vars))


## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
adlbh7 <- adlbh6 %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = AENTMTFL,
      by_vars = vars(STUDYID, USUBJID, PARAMCD),
      order = vars(ADT, VISITNUM),
      mode = "last"
    ),
    filter = (VISITNUM==12| (!is.na(AVISIT) & VISITNUM<12 & is.na(ABLFL)))
  ) %>% rowwise() %>%
  mutate(ALBTRVAL1= (1.5*A1HI)-LBSTRESN,
         ALBTRVAL2= LBSTRESN-(0.5*A1LO),
         # ALBTRVAL = ifelse(!is.na(ALBTRVAL1) & !is.na(ALBTRVAL2),max(across(c(ALBTRVAL1,ALBTRVAL2)), na.rm = T),NA_real_),
         CHG = ifelse(!is.na(ABLFL), NA, CHG)
  )


adlbh7 <- as.data.table(adlbh7)

adlbh7[, ALBTRVAL := pmax(ALBTRVAL1,ALBTRVAL2, na.rm = T)]

format_racen <- function(x) {
  case_when(
    x == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 6,
    x == 'ASIAN' ~ 3,
    x == 'BLACK OR AFRICAN AMERICAN' ~ 2,
    x == 'WHITE' ~ 1,
    TRUE ~ NA_real_
  )
}

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
adlbh8 <- adlbh7 %>%
  # Assign TRTA, TRTP
  # Create End of Treatment Record
restrict_derivation(
  derivation = derive_var_extreme_flag,
  args = params(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    order = vars(desc(ALBTRVAL), AVISITN),
    new_var = ANL01FL,
    mode = "first"
  ),
  filter = (is.na(ABLFL) & !is.na(AVISIT) & AVISITN<=24 & !is.na(ALBTRVAL))
) %>%
# filter(EOTFL == "Y") %>%
# mutate(
#   AVISIT = "End of Treatment",
#   AVISITN = 99
# ) %>%
# union_all(advs) %>%
# select(-EOTFL) %>%
mutate(
  TRTP = TRT01P,
  TRTPN = TRT01PN,
  TRTA = TRT01A,
  TRTAN = TRT01AN,
  RACEN=format_racen(RACE)
)

adlbh9 <- adlbh8 %>%
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    order = vars(ADT, AVISITN, AVAL),
    mode = "last",
    filter = (is.na(ABLFL) & AVISITN<=24 & !is.na(AVISIT)),
    set_values_to = vars(
      AVISIT = "End of Treatment",
      AVISITN = 99
    )
  ) %>%
  mutate(A1LO=ifelse(A1LO==0, NA_real_, A1LO),
         flg=ifelse(is.na(AVISITN) & AVISIT=='Baseline', 'Y', NA_character_)
         )

## Get ASEQ and AVALCATx and add PARAM/PARAMN ----
# advs <- advs %>%
#   # Calculate ASEQ
#   derive_var_obs_number(
#     new_var = ASEQ,
#     by_vars = vars(STUDYID, USUBJID),
#     order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN),
#     check_type = "error"
#   ) %>%
#   # Derive AVALCA1N and AVALCAT1
#   mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
#   derive_vars_merged(dataset_add = avalcat_lookup, by_vars = vars(PARAMCD, AVALCA1N)) %>%
#   # Derive PARAM and PARAMN
#   derive_vars_merged(dataset_add = select(param_lookup, -VSTESTCD), by_vars = vars(PARAMCD))


# Add all ADSL variables
# advs <- advs %>%
#   derive_vars_merged(
#     dataset_add = select(adsl, STUDYID, USUBJID, TRT01AN, TRT01PN),
#     by_vars = vars(STUDYID, USUBJID, TRT01AN, TRT01PN)
#   ) %>%
#   mutate(TRTPN=TRT01PN, TRTAN=TRT01AN, PCHG=janitor::round_half_up(PCHG,1),
#          CHG=janitor::round_half_up(CHG,1))

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

adlbh <- adlbh9 %>% filter(is.na(flg)) %>% drop_unspec_vars(adlbh_spec) %>%
  mutate(across(c(TRTSDT,TRTEDT,ADT),~ . - as.Date('1960-01-01')))

adlbh <- adlbh %>%
  check_variables(adlbh_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adlbh_spec) %>% # Orders the columns according to the spec
  sort_by_key(adlbh_spec) %>%
  xportr_label(adlbh_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(adlbh_spec) # %>%
  # xportr_format(adlbh_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADLBC")

# Save output ----
xportr_write(adlbh, "adam/adlbh.xpt")

# dir <- tempdir() # Change to whichever directory you want to save the dataset in
# saveRDS(advs, file = file.path(dir, "advs.rds"), compress = "bzip2")


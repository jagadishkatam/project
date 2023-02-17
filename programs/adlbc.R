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

sdtms <- c('dm','ds','ex','ae','lb','suppdm','suppds','sv','vs','qs','sc','supplb','suppae')

read_sdtm <- \(x){
  y <- haven::read_xpt(paste0("sdtm/",x,".xpt")) %>% convert_blanks_to_na()
  assign(x,y,envir = .GlobalEnv)
}

walk(sdtms, ~read_sdtm(.x))

# Lookup tables ----
metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adlbc_spec <- metacore %>% select_dataset('ADLBC')

supplb <- supplb %>% mutate(IDVARVAL=trimws(IDVARVAL))
lb1 <- metatools::combine_supp(lb, supplb) %>% as_tibble()

# lb1 <- create_var_from_codelist(lb, adlbc_spec, LBTESTCD, Term, decode_to_code = F)

paramcd <- adlbc_spec$codelist %>% as_tibble() %>% unnest(everything()) %>% filter(code_id=='PARAMCD_ADLBC') %>%
  transmute(LBTESTCD=code, PARAMCD=code, PARAM=decode)

paramn <- adlbc_spec$codelist %>% as_tibble() %>% unnest(everything()) %>% filter(code_id=='PARAMN_ADLBC') %>%
  transmute(PARAMN=as.numeric(code), PARAM=decode)

param_lookup <- paramcd %>%
  derive_vars_merged(
  dataset_add = paramn,
  by_vars = vars(PARAM)
)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(SUBJID, TRTSDT, TRTEDT, TRT01A, TRT01P,TRT01AN, TRT01PN, SITEID,AGE, AGEGR1N, RACE, RACEN, SEX, SAFFL, AGEGR1, COMP24FL, DSRAEFL)

adlbc1 <- lb1 %>% filter(LBCAT=='CHEMISTRY') %>%
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

adlbc2 <- adlbc1 %>%
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
    PARCAT1 = 'CHEM'
  )


## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adlbc3 <- adlbc2 %>%
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
adlbc4 <- adlbc3 %>%
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
adlbc5 <- adlbc4 %>%
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
  derive_var_pchg()


adlbc6 <- adlbc5 %>%
  derive_var_analysis_ratio(numer_var = AVAL,
                                  denom_var = BA1LO,
                                  new_var = BR2A1LO
                                  ) %>%
                     derive_var_analysis_ratio(
                                    numer_var = AVAL,
                                    denom_var = BA1HI,
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
                                 ) %>%
                    mutate(BR2A1LO=round(BR2A1LO,1),
                           BR2A1HI=round(BR2A1HI,1),
                           R2A1LO=round(R2A1LO,1),
                           R2A1HI=round(R2A1HI,1)
                           )

# adlbcx <- adlbc7 %>% select(!!!negate_vars(adsl_vars))


## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
adlbc7 <- adlbc6 %>%
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
  mutate(ALBTRVAL1= LBSTRESN-(1.5*A1HI),
         ALBTRVAL2= (0.5*A1LO) - LBSTRESN,
         ALBTRVAL = ifelse(!is.na(ALBTRVAL1) & !is.na(ALBTRVAL2),max(across(c(ALBTRVAL1,ALBTRVAL2)), na.rm = T),NA_real_),
         CHG = ifelse(!is.na(ABLFL), NA, CHG)
         )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
adlbc <- adlbc7 %>%
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
    TRTPN = TRT01PN,
    TRTA = TRT01A,
    TRTAN = TRT01AN,
    ANL01FL = NA_character_
  )

# adlbc <- adlbc %>%
#   derive_extreme_records(
#     by_vars = vars(STUDYID, USUBJID, PARAMCD),
#     order = vars(ADT, AVISITN, AVAL),
#     mode = "last",
#     filter = (ANL01FL == "Y" & AVISITN<=24),
#     set_values_to = vars(
#       AVISIT = "End of Treatment",
#       AVISITN = 99
#     )
#   )

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

adlbc <- adlbc %>% drop_unspec_vars(adlbc_spec)

adlbc <- adlbc %>%
  check_variables(adlbc_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adlbc_spec) %>% # Orders the columns according to the spec
  sort_by_key(adlbc_spec) %>%
  xportr_label(adlbc_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(adlbc_spec) %>%
  xportr_format(adlbc_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADLBC")

# Save output ----
xportr_write(adlbc, "adam/adlbc.xpt")

# dir <- tempdir() # Change to whichever directory you want to save the dataset in
# saveRDS(advs, file = file.path(dir, "advs.rds"), compress = "bzip2")


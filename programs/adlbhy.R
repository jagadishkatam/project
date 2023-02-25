# Name: ADLBHY
#
# Label: Analysis Dataset Lab Hy's Law
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
  # mutate(BR2A1LO=round(BR2A1LO,1),
  #        BR2A1HI=round(BR2A1HI,1),
  #        R2A1LO=round(R2A1LO,1),
  #        R2A1HI=round(R2A1HI,1)
  # )

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
         CHG = ifelse(!is.na(ABLFL), NA, CHG)
  )


# ALBTRVAL = ifelse(!is.na(ALBTRVAL1) & !is.na(ALBTRVAL2),max(across(c(ALBTRVAL1,ALBTRVAL2)), na.rm = T),NA_real_),

adlbc7 <- as.data.table(adlbc7)

adlbc7[, ALBTRVAL := pmin(ALBTRVAL1,ALBTRVAL2, na.rm = T)]


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

############################################################################################################################



adlbhy_spec <- metacore %>% select_dataset('ADLBHY')

# lb1 <- create_var_from_codelist(lb, adlbc_spec, LBTESTCD, Term, decode_to_code = F)

paramcd <- adlbhy_spec$codelist %>% as_tibble() %>% unnest(everything()) %>% filter(code_id=='PARAMCD_ADLBHY') %>%
  transmute(LBTESTCD=code, PARAMCD=code, PARAM=decode)

paramn <- adlbhy_spec$codelist %>% as_tibble() %>% unnest(everything()) %>% filter(code_id=='PARAMN_ADLBHY') %>%
  transmute(PARAMN=as.numeric(code), PARAM=decode)

param_lookup <- paramcd %>%
  derive_vars_merged(
    dataset_add = paramn,
    by_vars = vars(PARAM)
  )

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(SUBJID, TRTSDT, TRTEDT, TRT01A, TRT01P,TRT01AN, TRT01PN, SITEID,AGE, AGEGR1N, RACE, RACEN, SEX, SAFFL, AGEGR1, COMP24FL, DSRAEFL)

adlbhy1 <- adlbc %>% filter(PARAMCD %in% c('ALT','AST','BILI') & !is.na(AVISIT) & AVISITN<=24) # %>%
# Join ADSL with VS (need TRTSDT for ADY derivation)
# derive_vars_merged(
#   dataset_add = adsl,
#   new_vars = adsl_vars,
#   by_vars = vars(STUDYID, USUBJID)
# ) %>%
# ## Calculate ADT, ADY ----
# derive_vars_dt(
#   new_vars_prefix = "A",
#   dtc = LBDTC
# ) %>%
# derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlbhy2 <- adlbhy1 %>%
  convert_blanks_to_na() # %>%
# restrict_derivation(
#   derivation = derive_var_shift,
#   args = params(
#     new_var = SHIFT1,
#     from_var = BNRIND,
#     to_var = ANRIND
#   ),
#   filter = is.na(ABLFL) & !(PARAMCD %in% c('ALT','AST','BILI'))
# )


adlbhy3 <- adlbhy2 %>%
  derive_extreme_records(
    by_vars = vars(USUBJID, AVISITN),
    order = vars(AVAL),
    mode = "first",
    # filter = PARAMCD=='BILI' & AVAL>(1.5*A1HI) & !is.na(AVISIT),
    filter = PARAMCD=='BILI',
    set_values_to = vars(
      PARAMTYP = "DERIVED",
      PARAMCD = 'BILIHY',
      PARAM = 'Bilirubin 1.5 x ULN',
      PARAMN = 4,
      ADY = NA_real_,
      VISIT = NA_character_,
      VISITNUM = NA_real_,
      AVAL = 0,
      BASE = 0,
      A1LO = NA_real_,
      A1HI = NA_real_,
      PARCAT1 = 'HYLAW',
      # R2A1HI = NA_real_,
      R2A1LO = NA_real_,
      BR2A1LO = NA_real_,
      BR2A1HI = NA_real_
    )
  ) %>%
  derive_extreme_records(
    by_vars = vars(USUBJID,AVISITN),
    order = vars(AVAL),
    mode = "last",
    filter = PARAMCD %in% c('AST','ALT'),
    set_values_to = vars(
      PARAMTYP = "DERIVED",
      PARAMCD = 'TRANSHY',
      PARAM = 'Transaminase 1.5 x ULN',
      PARAMN = 5,
      ADY = NA_real_,
      VISIT = NA_character_,
      VISITNUM = NA_real_,
      AVAL = 0,
      BASE = 0,
      A1LO = NA_real_,
      A1HI = NA_real_,
      PARCAT1 = 'HYLAW',
      # R2A1HI = NA_real_,
      R2A1LO = NA_real_,
      BR2A1LO = NA_real_,
      BR2A1HI = NA_real_
    )
  )


format_shift1 <- function(x) {
  case_when(
    x == 'High to Normal' ~ 0,
    x == 'Normal to Normal' ~ 1,
    x == 'Normal to High' ~ 2,
    TRUE ~ NA_real_
  )
}

format_B <- function(x) {
  case_when(
    x == 0 ~ 'Normal',
    x == 1 ~ 'High',
    TRUE ~ NA_character_
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


BILIHY <- adlbhy3 %>% filter(PARAMCD %in% c('BILIHY'))  %>%
  select(-ADT, -BASE, -BNRIND, -ANRIND) %>%
  mutate(AVAL=ifelse(R2A1HI>1.5, 1, 0)) %>%
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>% mutate(BASE=ifelse(is.na(AVAL), NA_real_, BASE),
               BNRIND=format_B(BASE),
               ANRIND=format_B(AVAL)) %>%
  restrict_derivation(
    derivation = derive_var_shift,
    args = params(
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND
    ),
    filter = !(PARAMCD %in% c('ALT','AST','BILI'))
  )


TRANSHY <- adlbhy3 %>% filter(PARAMCD %in% c('TRANSHY')) %>%
  select(-ADT, -BASE, -BNRIND, -ANRIND) %>%
  mutate(AVAL=ifelse(R2A1HI>1.5, 1, 0)) %>%
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>% mutate(BASE=ifelse(is.na(AVAL), NA_real_, BASE),
               BNRIND=format_B(BASE),
               ANRIND=format_B(AVAL)) %>%
  restrict_derivation(
    derivation = derive_var_shift,
    args = params(
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND
    ),
    filter = !(PARAMCD %in% c('ALT','AST','BILI'))
  )


adlbhy4 <- adlbhy3 %>% filter(!(PARAMCD %in% c('BILIHY','TRANSHY')))

adlbhy5 <- bind_rows(adlbhy4,BILIHY,TRANSHY)

adlbhy <- adlbhy5 %>% drop_unspec_vars(adlbhy_spec) %>% group_by(USUBJID, PARAMCD, AVISITN) %>% slice_tail(n=1) %>%
  ungroup() %>%
  mutate(RACEN=format_racen(RACE)) %>%
  select(-PARAM,-PARAMN) %>%
  ## Add PARAMCD only - add PARAM etc later ----
derive_vars_merged_lookup(
  dataset_add = param_lookup,
  new_vars = vars(PARAMCD, PARAM, PARAMN),
  by_vars = vars(PARAMCD)
) %>%
  mutate(CRIT1=ifelse(!(PARAMCD %in% c('BILIHY','TRANSHY')), 'R2A1HI > 1.5', NA_character_),
         CRIT1FL=ifelse(!(PARAMCD %in% c('BILIHY','TRANSHY')) & R2A1HI>1.5, 'Y', NA_character_),
         CRIT1FL=ifelse(!is.na(CRIT1) & is.na(CRIT1FL), 'N', CRIT1FL),
         CRIT1FL=ifelse(is.na(AVAL), NA_character_, CRIT1FL),
         CRIT1FN=ifelse(!(PARAMCD %in% c('BILIHY','TRANSHY')) & R2A1HI>1.5, 1, NA_real_),
         CRIT1FN=ifelse(!is.na(CRIT1) & is.na(CRIT1FN), 0, CRIT1FN),
         CRIT1FN=ifelse(is.na(AVAL), NA_real_, CRIT1FN),
         SHIFT1N=ifelse(PARAMCD %in% c('BILIHY','TRANSHY'),format_shift1(SHIFT1),NA_real_),
         SHIFT1=ifelse(is.na(SHIFT1N),NA_character_,SHIFT1),
         R2A1HI=ifelse((PARAMCD %in% c('BILIHY','TRANSHY')), NA_real_, R2A1HI)
  )

hylaw1 <- adlbhy %>% filter(PARAMCD %in% c('BILIHY')) %>% arrange(USUBJID,AVISITN, AVAL) %>%
  group_by(USUBJID,AVISITN) %>% slice_head(n=1) %>% ungroup() %>% rename(AVAL_BILIHY=AVAL, BASE_BILIHY=BASE)

hylaw2 <- adlbhy %>% filter(PARAMCD %in% c('TRANSHY')) %>% arrange(USUBJID,AVISITN, AVAL) %>%
  group_by(USUBJID,AVISITN) %>% slice_head(n=1) %>% ungroup() %>% rename(AVAL_TRANSHY=AVAL, BASE_TRANSHY=BASE) %>%
  select(USUBJID, AVISITN, AVAL_TRANSHY, BASE_TRANSHY)

hylaw <- hylaw1 %>% inner_join(hylaw2, by=c('USUBJID','AVISITN')) %>%
  mutate(BASE_BILIHY=ifelse(is.na(BASE_BILIHY) & !is.na(BASE_TRANSHY), 0, BASE_BILIHY))

hylaw<- as.data.table(hylaw)

hylaw[, `:=`(AVAL = pmin(AVAL_BILIHY,AVAL_TRANSHY, na.rm = T), BASE = pmin(BASE_BILIHY,BASE_TRANSHY, na.rm = T))]

hylaw <- hylaw %>%
  mutate(PARAMTYP = "DERIVED",
         PARAM = 'Total Bili 1.5 x ULN and Transaminase 1.5 x ULN',
         PARAMCD = 'HYLAW',
         PARAMN = 6,
         BNRIND=format_B(BASE),
         ANRIND=format_B(AVAL)) %>%
  restrict_derivation(
    derivation = derive_var_shift,
    args = params(
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND
    ),
    filter = !is.na(PARAMCD)
  ) %>%
  mutate(SHIFT1N=format_shift1(SHIFT1),
         SHIFT1=ifelse(!is.na(SHIFT1) & is.na(SHIFT1N), NA_character_, SHIFT1))

adlbhy <- bind_rows(adlbhy, hylaw) %>% drop_unspec_vars(adlbhy_spec)

adlbhy <- adlbhy %>%
  check_variables(adlbhy_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adlbhy_spec) %>% # Orders the columns according to the spec
  sort_by_key(adlbhy_spec) %>%
  xportr_label(adlbhy_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(adlbhy_spec) %>%
  xportr_format(adlbhy_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADLBHY")


# Save output ----
xportr_write(adlbhy, "adam/adlbhy.xpt")


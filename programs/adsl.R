# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ex, ds
library(admiral) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)
library(purrr)
library(tidyverse)
library(xportr)
library(metacore)
library(metatools)
library(janitor)
# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

sdtms <- c('dm','ds','ex','ae','lb','suppdm','suppds','sv','vs','qs','sc','mh')

read_sdtm <- \(x){
  y <- haven::read_xpt(paste0("sdtm/",x,".xpt")) %>% convert_blanks_to_na()
  assign(x,y,envir = .GlobalEnv)
}

walk(sdtms, ~read_sdtm(.x))



# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

# User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.

# Grouping
format_racegr1 <- function(x) {
  case_when(
    x == "WHITE" ~ "White",
    x != "WHITE" ~ "Non-white",
    TRUE ~ "Missing"
  )
}

format_agegr1 <- function(x) {
  case_when(
    x < 65 ~ "<65",
    between(x, 65, 80) ~ "65-80",
    x > 80 ~ ">80",
    TRUE ~ "Missing"
  )
}

format_bmiblgr1 <- function(x) {
  case_when(
    x < 25 ~ "<25",
    25 <= x & x < 30 ~ "25-<30",
    x >= 30 ~ ">=30",
    TRUE ~ NA_character_
  )
}

format_agegr1n <- function(x) {
  case_when(
    x < 65 ~ 1,
    between(x, 65, 80) ~ 2,
    x > 80 ~ 3,
    TRUE ~ NA_real_
  )
}

format_durdsgr1 <- function(x) {
  case_when(
    x < 12 ~ "<12",
    x >= 12 ~ ">=12",
    TRUE ~ "Missing"
  )
}

format_region1 <- function(x) {
  case_when(
    x %in% c("CAN", "USA") ~ "NA",
    !is.na(x) ~ "RoW",
    TRUE ~ "Missing"
  )
}

format_lddthgr1 <- function(x) {
  case_when(
    x <= 30 ~ "<= 30",
    x > 30 ~ "> 30",
    TRUE ~ NA_character_
  )
}

format_racen <- function(x) {
  case_when(
    x == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 1,
    x == 'ASIAN' ~ 2,
    x == 'BLACK OR AFRICAN AMERICAN' ~ 3,
    x == 'WHITE' ~ 6,
    TRUE ~ NA_real_
  )
}

# EOSSTT mapping
format_eoxxstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    !(x %in% c("COMPLETED", "SCREEN FAILURE")) & !is.na(x) ~ "DISCONTINUED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "ONGOING"
  )
}

format_TRT01P <- function(x) {
  case_when(
    x == "Placebo" ~ 0,
    x == "Xanomeline High Dose" ~ 81,
    x == "Xanomeline Low Dose" ~ 54
  )
}

format_DCSREAS <- function(x) {
  case_when(
    x == "ADVERSE EVENT" ~ "Adverse Event",
    x == "STUDY TERMINATED BY SPONSOR" ~ "Sponsor Decision",
    x == "DEATH" ~ "Death",
    x == "WITHDRAWAL BY SUBJECT" ~ "Withdrew Consent",
    x == "PHYSICIAN DECISION" ~ "Physician Decision",
    x == "PROTOCOL VIOLATION" ~ "Protocol Violation",
    x == "LOST TO FOLLOW-UP" ~ "Lost to Follow-up",
    x == "LACK OF EFFICACY" ~ "Lack of Efficacy",
    T ~ NA_character_
  )
}

# Derivations ----
# impute start and end time of exposure to first and last respectively, do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )


ex_dt <- ex %>% group_by(USUBJID) %>% slice_tail(n=1) %>% select(USUBJID, EXENDTC) %>%
  mutate(TRTEDT_EXENDTC=as.Date(EXENDTC,'%Y-%m-%d')) %>% ungroup()


adsl <- dm %>% filter(ARMCD!='Scrnfail') %>%
  ## derive treatment variables (TRT01P, TRT01A) ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_adsl)
mutate(TRT01P = ARM, TRT01A = ARM) %>%
  ## derive treatment start date (TRTSDTM) ----
derive_vars_merged(
  dataset_add = ex_ext,
  filter_add = (EXDOSE > 0 |
                  (EXDOSE == 0 &
                     str_detect(EXTRT, "PLACEBO"))) &
    !is.na(EXSTDTM),
  new_vars = vars(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
  order = vars(EXSTDTM, EXSEQ),
  mode = "first",
  by_vars = vars(STUDYID, USUBJID)
) %>%
  ## derive treatment end date (TRTEDTM) ----
derive_vars_merged(
  dataset_add = ex_ext,
  filter_add = (EXDOSE > 0 |
                  (EXDOSE == 0 &
                     str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
  new_vars = vars(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
  order = vars(EXENDTM, EXSEQ),
  mode = "last",
  by_vars = vars(STUDYID, USUBJID)
) %>%
  ## Derive treatment end/start date TRTSDT/TRTEDT ----
derive_vars_dtm_to_dt(source_vars = vars(TRTSDTM, TRTEDTM))

## Disposition dates, status ----
# convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

# Screen fail date
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(SCRFDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE"
  ) %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(EOSDT = DSSTDT, DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  ) %>%
  # EOS status
  derive_var_disposition_status(
    dataset_ds = ds_ext,
    new_var = EOSSTT,
    status_var = DSDECOD,
    format_new_var = format_eoxxstt,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  ) %>%
  # Last retrieval date
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(FRVDT = DSSTDT),
    filter_add = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
  ) %>%
  # Derive Randomization Date
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(RANDDT = DSSTDT)
  ) %>%
  # Derive Randomization Date
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = VISITNUM>3 & DSCAT=='DISPOSITION EVENT',
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(DSDATE = DSSTDT)
  ) %>%
  # Death date - impute partial date to first day/month
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC,
    highest_imputation = "M",
    date_imputation = "first"
  ) %>%
  # Relative Day of Death
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  ) %>%
  # Elapsed Days from Last Dose to Death
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  ) %>%
  derive_vars_merged(
    dataset_add = ex_dt,
    by_vars = vars(USUBJID)
  ) %>%
  mutate(TRTEDT=case_when(
    is.na(TRTEDT) ~ DSDATE,
    is.na(TRTEDT_EXENDTC) ~ DSDATE,
    T ~ TRTEDT)) %>%
  ## derive treatment duration (TRTDURD) ----
derive_var_trtdurd()

## Last known alive date ----
ae_start_date <- date_source(
  dataset_name = "ae",
  date = AESTDT
)
ae_end_date <- date_source(
  dataset_name = "ae",
  date = AEENDT
)
lb_date <- date_source(
  dataset_name = "lb",
  date = LBDT,
  filter = !is.na(LBDT)
)
trt_end_date <- date_source(
  dataset_name = "adsl",
  date = TRTEDT
)

# impute AE start and end date to first
ae_ext <- ae %>%
  derive_vars_dt(
    dtc = AESTDTC,
    new_vars_prefix = "AEST",
    highest_imputation = "M"
  ) %>%
  derive_vars_dt(
    dtc = AEENDTC,
    new_vars_prefix = "AEEN",
    highest_imputation = "M"
  )

# impute LB date to first
lb_ext <- derive_vars_dt(
  lb,
  dtc = LBDTC,
  new_vars_prefix = "LB",
  highest_imputation = "M"
)

adsl <- adsl %>%
  derive_var_extreme_dt(
    new_var = LSTALVDT,
    ae_start_date, ae_end_date, lb_date, trt_end_date,
    source_datasets = list(ae = ae_ext, lb = lb_ext, adsl = adsl),
    mode = "last"
  ) %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = vars(STUDYID, USUBJID),
    new_var = SAFFL,
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  ) %>%
  ## Groupings and others variables ----
mutate(
  RACEGR1 = format_racegr1(RACE),
  AGEGR1 = format_agegr1(AGE),
  AGEGR1N = format_agegr1n(AGE),
  REGION1 = format_region1(COUNTRY),
  LDDTHGR1 = format_lddthgr1(LDDTHELD),
  TRT01PN=format_TRT01P(TRT01P),
  TRT01AN=format_TRT01P(TRT01A),
  DTH30FL = if_else(LDDTHGR1 == "<= 30", "Y", NA_character_),
  DTHA30FL = if_else(LDDTHGR1 == "> 30", "Y", NA_character_),
  DTHB30FL = if_else(DTHDT <= TRTSDT + 30, "Y", NA_character_),
  DOMAIN = NULL
)

visitnum <- derive_vars_dt(
  sv,
  new_vars_prefix = "VISIT1",
  dtc = SVSTDTC
) %>% filter(VISITNUM==1) %>% select(USUBJID, VISIT1DT)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(visitnum, USUBJID, VISIT1DT),
  by_vars = vars( USUBJID)
)

visnumen <- ds %>% filter(DSCAT=='DISPOSITION EVENT') %>%
  mutate(VISNUMEN=ifelse(DSTERM %in% c('PROTOCOL COMPLETED','ADVERSE EVENT') & VISITNUM==13, 12, VISITNUM)) %>%
  select(USUBJID, VISNUMEN)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(visnumen, USUBJID, VISNUMEN),
  by_vars = vars( USUBJID)
)

WEIGHT <- vs %>% filter(VSTESTCD=='WEIGHT' & VISITNUM==3) %>% group_by(USUBJID) %>% slice_tail(n=1) %>%
  mutate(WEIGHTBL=round(VSSTRESN,1)) %>% select(USUBJID, WEIGHTBL) %>% ungroup()

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(WEIGHT, USUBJID, WEIGHTBL),
  by_vars = vars( USUBJID)
)

HEIGHT <- vs %>% filter(VSTESTCD=='HEIGHT' & VISITNUM==1) %>%
  mutate(HEIGHTBL=round(VSSTRESN,1)) %>% select(USUBJID, HEIGHTBL)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(HEIGHT, USUBJID, HEIGHTBL),
  by_vars = vars( USUBJID)
)

trt <- sort(paste(adsl$SITEID,adsl$TRT01PN,sep='_'))

rle <- do.call(cbind,rle(trt)) %>% as_tibble() %>%
  tidyr::extract(values, into = c('SITEID','TRT01PN'), regex = '(.*)_(.*)', convert = T) %>%
  mutate(SITEGR1=ifelse(as.numeric(lengths)<3 | SITEID %in% c(715,717), '900', as.character(SITEID)),
         SITEID=as.character(SITEID),
         TRT01PN=as.numeric(TRT01PN)) %>%  select(-lengths)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(rle, SITEGR1, SITEID, TRT01PN),
  by_vars = vars( SITEID, TRT01PN)
)

MMSETOT <- qs %>% filter(QSCAT=="MINI-MENTAL STATE") %>% group_by(USUBJID) %>% summarise(MMSETOT=sum(QSSTRESN))

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(MMSETOT, USUBJID, MMSETOT),
  by_vars = vars( USUBJID)
)

EDUCLVL <- sc %>% filter(SCTESTCD=='EDLEVEL') %>% select(USUBJID, EDUCLVL=SCSTRESN)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(EDUCLVL, USUBJID, EDUCLVL),
  by_vars = vars( USUBJID)
)

dsdecod <- ds %>% filter(DSCAT=='DISPOSITION EVENT') %>% select(USUBJID, DCDECOD=DSDECOD)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(dsdecod, USUBJID, DCDECOD),
  by_vars = vars( USUBJID)
)

adsl <- adsl %>% mutate(DISCONFL=ifelse(DCDECOD!='COMPLETED', 'Y', NA_character_),
                        DSRAEFL=ifelse(DCDECOD=='ADVERSE EVENT','Y',NA_character_),
                        BMIBL=round(WEIGHTBL / ((HEIGHTBL/100)^2),1),
                        BMIBLGR1=format_bmiblgr1(BMIBL),
                        RACEN=format_racen(RACE))

DISONSDT <- mh %>% filter(MHTERM=="ALZHEIMER'S DISEASE") %>% mutate(DISONSDT=as.Date(MHSTDTC)) %>%
  select(USUBJID, DISONSDT)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(DISONSDT, USUBJID, DISONSDT),
  by_vars = vars( USUBJID)
)


adsl$DURDIS <- compute_duration(
  adsl$DISONSDT,
  adsl$VISIT1DT,
  in_unit = "days",
  out_unit = "months",
  floor_in = TRUE,
  add_one = TRUE,
  trunc_out = FALSE
)


adsl <- adsl %>%
  mutate(DURDIS=round(DURDIS,1),
         DURDSGR1=format_durdsgr1(DURDIS)
         )

COMPFL <- metatools::combine_supp(dm,suppdm) %>%
  transmute(COMP8FL=COMPLT8,COMP16FL=COMPLT16,COMP24FL=COMPLT24, USUBJID=USUBJID, EFFFL=EFFICACY, ITTFL=ITT)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(COMPFL, USUBJID, COMP8FL, COMP16FL, COMP24FL, ITTFL, EFFFL),
  by_vars = vars( USUBJID)) %>%
  mutate(COMP8FL=ifelse(is.na(COMP8FL),'N', COMP8FL),
             COMP16FL=ifelse(is.na(COMP16FL),'N', COMP16FL),
             COMP24FL=ifelse(is.na(COMP24FL),'N', COMP24FL),
             ITTFL=ifelse(is.na(ITTFL),'N', ITTFL),
             EFFFL=ifelse(is.na(EFFFL),'N', EFFFL)
)

dosconp <- ds %>% filter(DSCAT=='DISPOSITION EVENT') %>% derive_vars_dt(new_vars_prefix = "DISPX", dtc = DSSTDTC) %>%
  select(USUBJID, DISPXDT)

ex1 <- ex %>% derive_vars_dt(new_vars_prefix = "AST", dtc = EXSTDTC) %>%
  derive_vars_dt(new_vars_prefix = "AEN", dtc = EXENDTC) %>%
  select(USUBJID, VISIT, VISITNUM, EXDOSE, ASTDT, AENDT)

data1 <- ex1 %>% left_join(adsl %>% select(USUBJID, TRT01PN), by='USUBJID') %>%
  group_by(USUBJID) %>% mutate(row=row_number(), max=max(row)) %>% ungroup()

data2 <- data1 %>% left_join(dosconp, by='USUBJID') %>%
  mutate(AENDT=ifelse(row==max & is.na(AENDT) & !is.na(DISPXDT), DISPXDT, AENDT),
         AENDT=as.Date(AENDT, origin='1970-01-01')) %>%
  mutate(cum=(AENDT-ASTDT+1)*EXDOSE) %>% group_by(USUBJID) %>%
  summarise(CUMDOSE=as.numeric(sum(cum, na.rm = T))) %>%
  ungroup()

adsl <- adsl %>%
  derive_vars_merged(
  dataset_add = select(data2, USUBJID, CUMDOSE),
  by_vars = vars(USUBJID)
)

adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "RFEN",
    dtc = RFENDTC
  ) %>%
  mutate(DCSREAS=ifelse(DSTERM=='PROTOCOL ENTRY CRITERIA NOT MET', 'I/E Not Met', format_DCSREAS(DCDECOD)),
         AVGDD=round(CUMDOSE/TRTDURD,1)
         )

metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adsl_spec <- metacore %>% select_dataset('ADSL')

adsl <- adsl %>% drop_unspec_vars(adsl_spec)

# debugonce(xportr_type)
adsl <- adsl %>%
  check_variables(adsl_spec, dataset_name = "ADSL") %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adsl_spec, dataset_name = "ADSL") %>% # Orders the columns according to the spec
  sort_by_key(adsl_spec, dataset_name = "ADSL") %>%
  xportr_label(adsl_spec, domain = "ADSL") %>% # Assigns variable label from metacore specifications
  # xportr_type(adsl_spec, domain = "ADSL") %>%
  xportr_df_label(adsl_spec, domain = "ADSL") #%>%
  # xportr_format(adsl_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADSL")

# Save output ----
xportr_write(adsl, "adam/adsl.xpt")
# dir <- tempdir() # Change to whichever directory you want to save the dataset in
# saveRDS(adsl, file = file.path(dir, "adsl.rds"), compress = "bzip2")


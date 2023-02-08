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
    TRUE ~ "Missing"
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

adsl <- dm %>%
  ## derive treatment variables (TRT01P, TRT01A) ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_adsl)
mutate(TRT01P = ARM, TRT01A = ACTARM) %>%
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
derive_vars_dtm_to_dt(source_vars = vars(TRTSDTM, TRTEDTM)) %>%
  ## derive treatment duration (TRTDURD) ----
derive_var_trtdurd()

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
    new_vars = vars(EOSDT = DSSTDT),
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
  )

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
  TRT01PN=factor(TRT01P, labels = c(0,NA_real_,81, 54)),
  TRT01AN=factor(TRT01A, labels = c(0,NA_real_,81, 54)),
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

visnumen <- ds %>% filter(DSTERM=='PROTOCOL COMPLETED') %>%
  mutate(VISNUMEN=ifelse(VISITNUM==13, 12, VISITNUM)) %>% select(USUBJID, VISNUMEN)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(visnumen, USUBJID, VISNUMEN),
  by_vars = vars( USUBJID)
)

WEIGHT <- vs %>% filter(VSTESTCD=='WEIGHT' & VISITNUM==3) %>%  select(USUBJID, WEIGHTBL=VSSTRESN)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(WEIGHT, USUBJID, WEIGHTBL),
  by_vars = vars( USUBJID)
)


HEIGHT <- vs %>% filter(VSTESTCD=='HEIGHT' & VISITNUM==1) %>%  select(USUBJID, HEIGHTBL=VSSTRESN)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(HEIGHT, USUBJID, HEIGHTBL),
  by_vars = vars( USUBJID)
)

trt <- sort(paste(adsl$SITEID,adsl$TRT01PN,sep='_'))

rle <- do.call(cbind,rle(trt)) %>% as_tibble() %>%
  tidyr::extract(values, into = c('SITEID','TRT01PN'), regex = '(.*)_(.*)', convert = T) %>%
  mutate(SITEGR1=ifelse(as.numeric(lengths)<3, '999', as.character(SITEID)), SITEID=as.character(SITEID), TRT01PN=factor(TRT01PN)) %>%
  select(-lengths)

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
                        BMIBL=WEIGHTBL / ((HEIGHTBL/100)^2),
                        BMIBLGR1=format_bmiblgr1(BMIBL),
                        RACEN=format_racen(RACE))

DISONSDT <- mh %>% filter(MHTERM=="ALZHEIMER'S DISEASE") %>% mutate(DISONSDT=as.Date(MHSTDTC)) %>%
  select(USUBJID, DISONSDT)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(DISONSDT, USUBJID, DISONSDT),
  by_vars = vars( USUBJID)
)

adsl <- derive_vars_dy(adsl, reference_date = DISONSDT, source_vars = vars(DURDIS=VISIT1DT)) %>%
  mutate(DURDSGR1=format_durdsgr1(DURDIS))


COMPFL <- metatools::combine_supp(dm,suppdm) %>%
  transmute(COMP8FL=COMPLT8,COMP16FL=COMPLT16,COMP24FL=COMPLT24, USUBJID=USUBJID, EFFFL=EFFICACY, ITTFL=ITT)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(COMPFL, USUBJID, COMP8FL, COMP16FL, COMP24FL, ITTFL, EFFFL),
  by_vars = vars( USUBJID)
)


visit4 <- sv %>% filter(VISIT=='WEEK 4') %>% derive_vars_dt(new_vars_prefix = "VISIT4", dtc = SVENDTC) %>%
  select(USUBJID, VISIT4DT)

visit12 <- sv %>% filter(VISIT=='WEEK 12') %>% derive_vars_dt(new_vars_prefix = "VISIT12", dtc = SVENDTC) %>%
  select(USUBJID, VISIT12DT)

visit26 <- sv %>% filter(VISIT=='WEEK 26') %>% derive_vars_dt(new_vars_prefix = "VISIT26", dtc = SVENDTC) %>%
  select(USUBJID, VISIT26DT)

dosconp <- ds %>% filter(DSCAT=='DISPOSITION EVENT') %>% derive_vars_dt(new_vars_prefix = "DISPX", dtc = DSSTDTC) %>%
  select(USUBJID, DISPXDT)

cumdose <- visit4 %>% full_join(visit12, by='USUBJID') %>% full_join(visit26, by='USUBJID') %>%
  full_join(dosconp, by='USUBJID') %>% full_join(adsl %>% select(USUBJID, TRTSDT, TRTEDT, TRT01PN, TRTDURD), by='USUBJID') %>%
  mutate(first=ifelse(!(DISPXDT<=VISIT4DT),VISIT4DT-TRTSDT+1, TRTEDT-TRTSDT+1),
         second=ifelse(!(DISPXDT<=VISIT12DT),VISIT12DT-TRTSDT+1, TRTEDT-TRTSDT+1),
         third=ifelse((DISPXDT>VISIT12DT),TRTEDT-VISIT12DT+1, NA),
         CUMDOSE=as.numeric(TRT01PN)*rowSums(across(c(first,second,third)), na.rm=T),
         AVGDD=CUMDOSE/TRTDURD) %>%
  select(USUBJID, CUMDOSE, AVGDD)

adsl <- derive_vars_merged(
  adsl,
  dataset_add = select(cumdose, USUBJID, CUMDOSE, AVGDD),
  by_vars = vars( USUBJID)
)

adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "RFEN",
    dtc = RFENDTC
  )

metacore <- metacore::spec_to_metacore('/cloud/project/metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adsl_spec <- metacore %>% select_dataset('ADSL')

adsl <- adsl %>% drop_unspec_vars(adsl_spec) %>% mutate(DCSREAS=NA_character_)

adsl <- adsl %>%
  check_variables(adsl_spec) %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adsl_spec) %>% # Orders the columns according to the spec
  sort_by_key(adsl_spec) %>%
  xportr_label(adsl_spec) %>% # Assigns variable label from metacore specifications
  xportr_df_label(adsl_spec)

# Save output ----
xportr_write(adsl, "adam/adsl.xpt")
# dir <- tempdir() # Change to whichever directory you want to save the dataset in
# saveRDS(adsl, file = file.path(dir, "adsl.rds"), compress = "bzip2")


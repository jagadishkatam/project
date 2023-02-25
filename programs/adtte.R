# Name: ADTTE
#
# Label: AE Time To 1st Derm. Event Analysis
#
# Input: adsl, ae
library(admiral)# Contains example datasets from the CDISC pilot project
library(tidyverse)
library(lubridate)
library(stringr)
library(purrr)
library(xportr)
library(metacore)
library(metatools)

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
adsl <- haven::read_xpt("adam/adsl.xpt")
adae <- haven::read_xpt("adam/adae.xpt") %>% filter(CQ01NAM=='DERMATOLOGIC EVENTS' & AOCC01FL=='Y') %>%
  select(STUDYID, USUBJID, ASTDT, AENDT, TRTEMFL, AESEQ)


format_racen <- function(x) {
  case_when(
    x == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 6,
    x == 'ASIAN' ~ 3,
    x == 'BLACK OR AFRICAN AMERICAN' ~ 2,
    x == 'WHITE' ~ 1,
    TRUE ~ NA_real_
  )
}

# format_fmt <- function(x) {
#   case_when(
#     x == 'AGE' ~ '3.',
#     x == 'AGEGR1' ~ "$5.",
#     x == 'AGEGR1N' ~ '3.',
#     x == 'EVNTDESC' ~ '$25.',
#     x == 'PARAM' ~ '$32.',
#     x == 'PARAMCD' ~ '$4.',
#     x == 'RACE' ~ '$32.',
#     x == 'RACEN' ~ '3.',
#     x == 'SAFFL' ~ '$1.'
#   )
# }

metacore <- metacore::spec_to_metacore('metadata/specs.xlsx', where_sep_sheet = F, quiet = T)

adtte_spec <- metacore %>% select_dataset('ADTTE')

#adtte_spec$var_spec$format <- adtte_spec$var_spec %>% mutate(format=format_fmt(variable))

lstalv <- censor_source(
  dataset_name = "adsl",
  date = RFENDT,
  set_values_to = vars(
    EVNTDESC = "Study Completion Date",
    SRCDOM = "ADSL",
    SRCVAR = "RFENDT")
)

adtte1 <- derive_param_tte(
  dataset_adsl = adsl,
  start_date = TRTSDT,
  event_conditions = list(ae_event),
  censor_conditions = list(lstalv),
  source_datasets = list(adsl = adsl, adae = adae),
  set_values_to = vars(PARAMCD = "TTDE",
                       PARAM = "Time to First Dermatologic Event"
  )
)

adtte2 <- derive_vars_duration(
  adtte1,
  new_var = AVAL,
  start_date = STARTDT,
  end_date = ADT
)

adtte <- derive_vars_merged(
  adtte2,
  dataset_add = adsl,
  new_vars =vars(STUDYID, SITEID, USUBJID, AGE, AGEGR1, AGEGR1N, RACE, SEX, TRTSDT, TRTEDT, SAFFL, TRTDURD, TRT01P, TRT01A, TRT01AN),
  by_vars = vars(STUDYID, USUBJID)
) %>%
  mutate(TRTP=TRT01P, TRTA=TRT01A, TRTAN=TRT01AN, TRTDUR=TRTDURD,
         EVNTDESC = ifelse(CNSR==0,"Dematologic Event Occured",EVNTDESC),
         RACEN=format_racen(RACE)
  )

adtte <- adtte %>% drop_unspec_vars(adtte_spec)

adtte <- adtte %>%
  check_variables(adtte_spec, dataset_name = "ADTTE") %>% # Check all variables specified are present and no more
  # check_ct_data(adsl_spec, na_acceptable = T) %>% # Checks all variables with CT only contain values within the CT
  order_cols(adtte_spec, dataset_name = "ADTTE") %>% # Orders the columns according to the spec
  sort_by_key(adtte_spec, dataset_name = "ADTTE") %>%
  xportr_length(adtte_spec) %>%
  xportr_label(adtte_spec, domain = "ADTTE") %>% # Assigns variable label from metacore specifications
  xportr_df_label(adtte_spec, domain = "ADTTE")  %>%
   xportr_format(adtte_spec$var_spec %>% mutate_at(c("format"), ~replace_na(.,'')), domain = "ADTTE")

# Save output ----

xportr_write(adtte, "adam/adtte.xpt")


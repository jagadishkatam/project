library(haven)
library(admiral)
library(dplyr)
library(tidyr)
library(metacore)
library(metatools)
library(xportr)

dm <- read_xpt("sdtm/dm.xpt")
ex <- read_xpt("sdtm/ex.xpt")

xportr_write(adsl, "adam/adsl.xpt")

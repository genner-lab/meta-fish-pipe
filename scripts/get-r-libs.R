#!/usr/bin/env Rscript

libs <- sort(unique(unlist(lapply(list.files(path="scripts",pattern="\\.R"),function(x) grep("^suppressMessages\\(library\\(.+\\)\\)",readLines(here::here("scripts",x)),value=TRUE)))))
writeLines(libs,here::here("temp/libs.R"))
source(here::here("temp/libs.R"))
si <- sessionInfo()
date <- Sys.Date()
sink(here::here("temp/sessionInfo.txt"))
writeLines(paste("--------------------------------------------------\nR Package version list. Today's date is",date,"\n--------------------------------------------------\n\n"))
print(si)
sink()

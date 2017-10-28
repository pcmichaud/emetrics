clear all
capture log close
set more off

cd ~/emetrics/life

* using HMD data for Canada (both sexes) 2011
infile year age qx mx wx lx dx llx qqx ex using mortality.dat

gen logmx = log(mx)
reg logmx age

* example of a good series of do files (and output structure)

* start off by cleaning what is in memory and setting global options
clear all
capture log close
set more off

* go to the root of your project directory
cd ~/emetrics/stata-do

* setup a log file, in text form
mkdir logs
log using logs/amazingdo.txt, text replace

* put your raw data in a raw directory: never write there
mkdir raw
* run do-files which run particular sub-parts of your analysis
do do/prepare.do

use raw/data.dta, clear

* describe what your dataset contains
describe

* summarize variables
sum

* summarize with some details
sum, details

* tabulate data
tab age educ
* add some cell percentages for columns, take out freq
tab age educ, col nofreq
* do some means of earnings by age
tabstat earn, by(educ)
* look at medians and other quantiles
tabstat earn, by(educ) statistics(p10 p25 p50 p75 p90)






* clean up after your done (temp and other folders you do not need)
* here we are erasing everything to keep this archive clean, leave no file behind...
erase logs raw data tables figures temp

* close your log file
capture log close
exit

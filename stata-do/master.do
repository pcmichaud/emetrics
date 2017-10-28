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

do do/prepare.do


use raw/


* clean up after your done (temp and other folders you do not need)
* here we are erasing everything to keep this archive clean
erase logs raw data tables figures temp

* close your log file
capture log close
exit

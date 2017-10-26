* read simulated data from module

clear all
set more off

cd ~/emetrics/module

infile id age str1(married) consumption using output/simulatedpop.dat, clear

* describe variables
label var id "identifier"
label var age "age"
label var married "marital status (1=married)"
label var consumption "simulated consumption"
	gen marr = 1 if married=="T"
	replace marr = 0 if married=="F"
	label def marr 1 "married" 0 "not married"
	label values marr marr
	label var marr "marital status"
tsset id age

des

* summarize variables
sum 

* doing tabulations
tab married
tab age married, row nofreq

* doing stats by married status
tabstat consumption, by(married) statistics(p10 p25 p50 p75 p90 mean)

* doing a table of stats
table age married, content(mean consumption)

* doing nice graphs, first collapsing data
preserve
	collapse consumption, by(age marr)
	tsset marr age
	# d ;
	xtline consumption, overlay name(cons) xlabel(20(5)100) ytitle("consumption");
	# d cr
	graph export "figures/consumption.png", as(png) replace
restore



exit



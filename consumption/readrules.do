clear all
capture log close
set more off


cd ~/emetrics/consumption

infile cash1 cons1 cash5 cons5 cash10 cons10 using rules.dat, clear

foreach var of varlist * {
	replace `var' = `var'*1.0e-3
}
#d ;
twoway (line cons1 cash1) (line cons5 cash5) (line cons10 cash10), 
	 xtitle("cash-on-hand (in thousands)") ytitle("consumption (in thousands)")
	legend(label(1 lowest shock) label(2 mean shock) label(3 highest shock));
	#d cr


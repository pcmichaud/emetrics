clear all
capture log close
set more off

cd ~/emetrics/life

infile id age shock income wealth cash cons using simulated.dat, clear


label var id "id"
label var age "age"
label var shock "wage shock"
label var income "income"
label var wealth "wealth"
label var cash "cash-on-hand"
label var cons "consumption"
foreach var of varlist income wealth cash cons {
	replace `var' = `var'*1.0e-3
}

collapse cons income wealth cash, by(age)
drop if age>90
#d ;
twoway
	(line cons age)
	(line income age)
	(line wealth age), 
	legend(label(1 "consumption") label(2 "income") label(3 "wealth"))
	xtitle("age")
	ytitle("thousands");
#d cr
graph export "simulated.png", as(png) replace


clear all
capture log close
set more off

cd ~/emetrics/life

infile age shock wage cash cons value using rules.dat, clear


label var age "age"
label var shock "# wage shock"
label var wage "wage in dollars"
label var cash "cash-on-hand"
label var cons "consumption"
label var value "EPV of utility"

foreach var of varlist wage cash cons {
	replace `var' = `var'*1.0e-3
}
replace value = value * 1.0e6

local e = 5
#d ;
twoway (line cons cash if age==25&shock==`e') 
	   (line cons cash if age==45&shock==`e') 
	   (line cons cash if age==65&shock==`e')
	   (line cons cash if age==75&shock==`e'),
	   legend(label(1 "age 35") label(2 "age 45") 
	   label(3 "age 65") label(4 "age 75")) name(cons) nodraw;
#d cr	   

#d ;
twoway (line value cash if age==25&shock==`e') 
	   (line value cash if age==45&shock==`e') 
	   (line value cash if age==65&shock==`e')
	   (line value cash if age==75&shock==`e'),
	   legend(label(1 "age 35") label(2 "age 45") 
	   label(3 "age 65") label(4 "age 75")) name(value) nodraw;
#d cr

graph combine cons value, rows(2)
graph export rules.png, as(png) replace


	   

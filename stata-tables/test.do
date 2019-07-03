capture log close
clear all
set more off

cd ~/cedia/emetrics/stata-tables

set obs 10000

gen age = floor(runiform()*100)
gen male = runiform()<0.5
gen income = exp(log(34e3) + log(5)*rnormal() - 0.5*log(5)^2)

global vlist "age male income"
global labnames "Age Gender "Household Income""

* Package estout : net install st0085_2.pkg
eststo fem: estpost summarize age income if male==0
eststo men: estpost summarize age income if male==1
#d ;
esttab fem men using stats.tex, 
	replace main(mean %6.2f) aux(sd) mtitle("female" "male") nonumbers nonotes;
#d cr
* table full force...
local i = 1
file open table using "means.tex", write replace text
file write table "\begin{tabular}{lrrrr} " _n
file write table "\hline \hline "
file write table " & mean ($\mu$) & sd ($\sigma$) & min & max \\" _n
foreach var of varlist $vlist {
local lab : word `i' of $labnames
sum `var'
di "`lab'"
#d ;
file write table "`lab'"  " & " %7.3f (r(mean)) " & " %7.3f  (r(sd)) " & " 
%7.3f (r(min)) " & " %7.3f (r(max)) " \\"  _n ;
#d cr

local ++i
}
file write table "\hline \hline "
file write table "\end{tabular}" _n
file close table

* regression tables
eststo female: reg income age if male==0
eststo male: reg income age if male==1

esttab female male using reg.tex, se ar2 nonumbers nonotes mtitle("female" "male") replace





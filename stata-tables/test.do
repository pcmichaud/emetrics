capture log close
clear all
set more off

set obs 10000

gen age = floor(runiform()*100)
gen male = runiform()<0.5
gen income = exp(log(34e3) + log(5)*rnormal() - 0.5*log(5)^2)

global vlist "age male income"
global labnames "Age Gender "Household Income""

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

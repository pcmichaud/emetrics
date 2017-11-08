* sub do-file to prepare some data

* set number of observations
set obs 10000

* generate an id (could be used in some settings)
gen id = _n
* generate age
gen age = floor(runiform()*100)
* generate education
gen d = runiform()
gen educ = 1 if d < 0.33
replace educ = 2 if d >=0.33 & d<0.66
replace educ = 3 if d>=0.66
* generate a continuous earnings variable (when you use educ==2, you are creating a dummy equal to 1 if logical exp is true, zero if false)
gen earn = dexp(0.05*age - 0.001*age^2)*(20.0e3*(educ==1) + 40.0e3*(educ==2) + 60.0e3*(educ==3))*dexp(rnormal() - 0.5)
* if you want to to-code 







use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-severity.dta", clear
bysort patid: gen obs = _n
bysort patid: gen maxobs = _N
*keep if maxobs == 2
*keep if patid == 2885
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempSEV.dta", replace

use "getmatchedcohort-eczema-main-mhealth.dta", clear
*keep if patid == 2885
merge m:1 patid using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-age-sex-gp.dta", keep(match) 
gen date=indexdate
format date %td
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp.dta", replace

use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-BMI-all.dta", clear
*keep if patid == 2885
*rename dobmi date
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempBMI.dta", replace

use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-harmfulalcohol.dta", clear
*keep if patid == 2885
rename eventdate date
gen alc = 1
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempALC.dta", replace

use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-asthma.dta", clear
*keep if patid == 2885
rename eventdate date
gen asthma = 1
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempASTHMA.dta", replace

use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-sleep-all.dta", clear
*keep if patid == 2885
rename eventdate date
gen sleep = 1
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempSLEEP.dta", replace













use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp.dta", clear
* merge in setid
merge 1:1 patid indexdate using "getmatchedcohort-eczema-main-mhealth.dta" , keepusing(setid) nogen
drop _merge
/*
duplicates report patid setid
duplicates report patid indexdate
*/
* merge in BMI
merge m:1 patid using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempBMI.dta"

/*
	preserve
	use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempBMI.dta", clear
	l if patid==1167
	l in 1/5
	restore
 * add in BMI as the recorded value
bysort patid: egen bmi2 = median(bmi)
*/


* merge in Smok
drop _merge
merge m:1 patid using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-smoke-all.dta"
tab smokstatus, m

* merge in Ethnicitiy
drop _merge
merge m:1 patid using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/variables-ecz-ethnicity.dta"
tab eth5, m

* mrge in Carstairs & RUC
drop _merge
sort setid patid date

save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp-constant.dta" , replace

* Time updated variables 
use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp.dta", clear
keep setid patid indexdate date
append using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempALC.dta"
append using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempSEV.dta"
append using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempASTHMA.dta"
append using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/tempSLEEP.dta"
// steroids

drop obs maxobs
local varlist " "alc" "modsevere" "asthma" "sleep"" 
//local varlist "setid"
sort patid date
foreach xx of local varlist {
	cap drop test
	gen test = `xx'
	replace test = test[_n-1] if test == . & patid[_n-1] == patid[_n]
	replace test = 0 if test == .
	drop `xx'
	rename test `xx' 
}
// what about Time updated if event < indexdate? 

label define yesno 0 "no" 1 "yes" , replace
label define modsevere 0 "mild" 1 "moderate" 2 "severe" , replace
label val modsevere modsevere
label val alc asthma sleep yesno

gen dateTUC = date
format dateTUC %td
save "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp-timeupdated.dta", replace

* Merge the 2 together
use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp-constant.dta" , clear

preserve 
use "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp-timeupdated.dta" ,clear
l if patid==2292022
restore

merge m:m patid using "/Users/lsh1510922/Documents/Postdoc/2021_extract/out/temp-timeupdated.dta" 
sort setid patid date



// Persisting problems 
/* 
	- how to deal with date recorded < indexdate
	- merge m:m not sustainable. Lots of problems like TimeUpdated `date' value not merging across. Need start stop times
	- some patids in time updated file but not in master - huh?
	- only keeping one indexdate in  bmi code (but fine because merge on patid only) 
	
*/
	

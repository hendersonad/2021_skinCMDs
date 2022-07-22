//anxiety
use "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/out/variables-ecz-anxiety-all.dta", clear
append using "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/out/variables-pso-anxiety-all.dta"

collapse (count) patid, by(readcode)
merge 1:1 readcode using "/Users/lsh1510922/Documents/Postdoc/2021_SkinEpiExtract/codelists/medcodes-anxiety-nohistory.dta"
drop if _merge==2

gsort -patid 
local no20 =patid[20]
di `no20'
drop if patid < `no20'
keep readterm readcode patid  
order readterm 
l

export delim "/Users/lsh1510922/Documents/Postdoc/2021_skinCMDs/out/supplementary/top20_anxiety.csv", replace


//depression
use "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/out/variables-ecz-depression-all.dta", clear
append using "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/out/variables-pso-depression-all.dta"

collapse (count) patid, by(readcode)
merge 1:1 readcode using "/Users/lsh1510922/Documents/Postdoc/2021_SkinEpiExtract/codelists/medcodes-depression-nohistory.dta"
drop if _merge==2

gsort -patid 
local no20 =patid[20]
di `no20'
drop if patid < `no20'
keep readterm readcode patid  
order readterm 

export delim "/Users/lsh1510922/Documents/Postdoc/2021_skinCMDs/out/supplementary/top20_depression.csv", replace

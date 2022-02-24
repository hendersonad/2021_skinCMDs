version 15
clear all
capture log close

ssc install unique
* find path file location and run it

skinepipaths_v2

local files Patient 
local practot = 897
global filename "extract-${exposure}FlowDiagram"

* open log file
log using "${pathLogs}/${filename}", text replace

adopath + "${pathShareAdo}"
adopath + "${pathPrograms}"

dib "${exposure}", stars

// total with ${exposure} diagnosis code
qui {
use "$pathIn/Patient_extract_${ABBRVexp}_extract_1", clear
unique patid 
local ${ABBRVexp}Dx = r(unique)
}
di "total with ${exposure} diagnosis code = `${ABBRVexp}Dx'"

// total that do not meet exposed algorithm 
qui {
use ${pathOut}/${exposure}Exposed, clear
unique patid 
local ${ABBRVexp}Exposed = r(unique)
local diff = `${ABBRVexp}Dx'-`${ABBRVexp}Exposed'
}
di "total that meet ${exposure} algorithm = `${ABBRVexp}Exposed'"
di "total that DO NOT meet ${exposure} algorithm = `diff'"

// total that do not contribute follow up time
qui {
use  "${pathOut}/${exposure}Exposed-eligible-mhealth", clear
unique patid 
local ${ABBRVexp}Eligible = r(unique)
local diff = `${ABBRVexp}Exposed'-`${ABBRVexp}Eligible'
}

di "total that DO NOT contribute follow up time = `diff'"
di "total eligble for matching = `${ABBRVexp}Eligible'"

// total that do not have eligible match
qui{
use "${pathOut}/getmatchedcohort-${exposure}-main-mhealth", clear
keep if exposed == 1
unique patid 
local ${ABBRVexp}Matched = r(unique)
local diff = `${ABBRVexp}Eligible' - `${ABBRVexp}Matched'
}

di "total that DO NOT have a match = `diff'"
di "total matched for matching = `${ABBRVexp}Matched'"

log close

local user = c(username)
global user `user'
cd "C:/Users/$user/Dropbox/my_did_testing/"

capture mkdir juan
capture mkdir juan/tables
capture mkdir juan/graphs 
capture mkdir juan/dtafiles
capture mkdir juan/graphs/throw
tempfile temp 

/*
net install csdid, from ("https://raw.githubusercontent.com/friosavila/csdid_drdid/main/code/") replace
ssc install bacondecomp
ssc install did_multiplegt, replace
ssc install did_imputation, replace
ssc install event_plot, replace
net install stackedev, from ("https://raw.githubusercontent.com/joshbleiberg/stackedev/main/")
ssc install did2s
ssc install estout, replace
ssc install reghdfe
ssc install ftools
ssc install synth 
*/


use "Processed Data/actual state level turnout 2006-2018", clear
merge 1:1 state year using "Processed Data/voter id 2006 to 2018 fixed.dta"
keep if _m==3
drop _m 

drop if state=="Alabama" /* Alabama is a switcher */

qui tab state, gen(ss)
qui tab year, gen(yy)
forvalues state=1/49 {
	forvalues year=1/7 {
		qui gen sy_`state'_`year' = ss`state'*yy`year'
	}
}

qui encode state, gen(state_float)

gen event_stricty=.
forvalues x = 1/50{
	qui summ year if stricty==1 & state_float==`x'
	local y0 = r(min)
	local obs =r(N)
	if `obs'!=0{
		qui replace event_stricty = year - `y0' if state_float==`x'
	}
}

qui summ year
local y0 = r(max)
gen first_treat_strict =.
bysort state: egen tst = mean(stricty)

replace first_treat_strict =0 if tst==0 
forvalues x = 1/50{
	qui summ year if stricty==1 & state_float==`x'
	local y0 = r(min)
	qui replace first_treat_strict = `y0' if state_float==`x' & first_treat_strict==.
}


tab year first_treat_strict, missing
list state year stricty if first_treat_strict==., clean

tab year first_treat_strict

gen treat_states =""
levelsof state, local(levels) 
foreach x of local levels {
	qui summ year if state=="`x'" & tst>0
	qui replace treat_states="`x'" if state=="`x'" & tst>0 & r(N)!=0
}

gen never_treat = 1 if tst==0
replace never_treat = 0 if tst!=0

gen f_treat_nmiss = first_treat_strict /*  never-treated units as blank for treatment */
replace f_treat_nmiss =. if never_treat==1

save "dtafiles/grimmer_juan.dta", replace
********************************************************************************
*Callaway and Sant'Anna (2020)**************************************************
********************************************************************************
use "dtafiles/grimmer_juan.dta", clear

xtset state_float year 
bacondecomp VEPhighestoffice stricty, ddetail
gr export "juan/graphs/VEPhighestoffice_bacon.png", replace

csdid VEPhighestoffice stricty,  ivar(state_float)  time(year) gvar(first_treat_strict)
estat pretrend/* can reject null of pre-treatment = 0 */
estat simple
return list
display "No States Dropped, Estimate: " r(b)[1,1] " SE: " r(table)[2,1]

levelsof state, local(levels) 
foreach x of local levels {
	qui csdid VEPhighestoffice stricty if state!="`x'",  ivar(state_float)  time(year) gvar(first_treat_strict)
	estat pretrend
	qui estat simple
	display "State Dropped: `x' Estimate: " r(b)[1,1] " SE: " r(table)[2,1] " P_Val: " r(table)[4,1]
	display ""
}

csdid VEPhighestoffice stricty,  ivar(state_float)  time(year) gvar(first_treat_strict)
estat simple
local csdid_coe = round(r(b)[1,1],0.001)
local csdid_se = round(r(table)[2,1],0.001)
gen coe_1 = round(r(b)[1,1],0.001)
gen se_1 = round(r(table)[2,1],0.001)

*Chaisemartin & H..
did_multiplegt VEPhighestoffice state_float year stricty, placebo(1) breps(50) cluster(state_float)
local chaise_coe = round(e(effect_0),0.001)
local chaise_se = round(e(se_effect_0),0.001)
gen coe_2 = round(e(effect_0),0.001)
gen se_2 = round(e(se_effect_0),0.001)
gr export "juan/graphs/VEPhighestoffice_chaisemartin.png" , replace

*Borusyak 2021
* h did_imputation
did_imputation VEPhighestoffice state_float year first_treat_strict, fe(state_float) autosample
eststo reg2  
estout reg2 , cells(b(star fmt(5)) se(par fmt(5)))
local boy_coe = round(r(coefs)[1,1],0.001)
local boy_se = round(r(coefs)[1,2],0.001)
gen coe_3 = round(r(coefs)[1,1],0.001)
gen se_3 = round(r(coefs)[1,2],0.001)
/* https://github.com/borusyak/did_imputation/blob/main/five_estimators_example.do */

*Gardner 2021
did2s VEPhighestoffice, first_stage(state_float year) second_stage(first_treat_strict) treatment(first_treat_strict) cluster(state_float)
eststo reg4  
estout reg4 , cells(b(star fmt(5)) se(par fmt(5)))
local gar_coe = round(r(coefs)[1,1],0.001)
local gar_se = round(r(coefs)[1,2],0.001)
gen coe_4 = round(r(coefs)[1,1],0.001)
gen se_4 = round(r(coefs)[1,2],0.001)
*did2s VEPhighestoffice, first_stage(state_float year) second_stage(f_treat_nmiss) treatment(f_treat_nmiss) cluster(state_float)

********************************************************************************
*STACKED REGRESSIONS************************************************************
********************************************************************************
preserve
	replace event_stricty =0 if event_stricty==.
	qui tab event_stricty, gen(ee) 
	drop ee5 /* reference election prior */
	ren ee6 yr1
	
	stackedev VEPhighestoffice ee* yr1, cohort(f_treat_nmiss) time(year) never_treat(never_treat) clust_unit(state_float) unit_fe(state_float)
	eststo reg5  
	estout reg5 , cells(b(star fmt(5)) se(par fmt(5)))
	local stack_coe = round(r(coefs)[11,1],0.001)
	local stack_se = round(r(coefs)[11,2],0.001)
restore 

gen coe_5 = `stack_coe'
gen se_5 = `stack_se'


display ///
"Estimator: CSDID CHAISE BORUSYAK GARDNER STACKED" _n ///
"Beta: `csdid_coe' `chaise_coe' `boy_coe' `gar_coe' `stack_coe'" _n ///
"SE: `csdid_se' `chaise_se' `boy_se' `gar_se' `stack_se'" 

tempfile temp 
preserve
collapse (mean) coe* 
xpose , clear
ren v1 coefs
gen estimator = .
forvalues i=1/5{
	replace estimator = `i' if _n==`i'
}
save `temp', replace
restore

collapse (mean) se_*
xpose , clear 
ren v1 ses 
gen estimator = .
forvalues i=1/5{
	replace estimator = `i' if _n==`i'
}
merge 1:1 estimator using `temp'

gen ul = coefs+1.96*ses
gen ll = coefs-1.96*ses

twoway scatter coefs estimator, mc(blue) || ///
	rcap ul ll estimator , lc(blue) ||, ///
	xlabel(1 "CSDID" 2 "Chaisemartin" 3 "Borusyak" 4 "Gardner" 5 "Stacked") ///
	yline(0, lp(-) lc(black)) ti(CCSS State-Level) ///
	legend(off) xti("Various Estimators") yti("Estimate") ///
	graphregion(color(white)) bgcolor(white) 
	gr export "juan/graphs/ccse_estimations.png", replace

********************************************************************************
*DISAGGREGATED
********************************************************************************
use "dtafiles/grimmer_juan.dta", clear
ren state stname
merge m:1 stname using "state fips codes.dta"
keep if _m==3 
drop _m

ren st fips  
ren stname state 
sort fips year 

tempfile ops 
preserve
	use "C:\Users\\$user\Dropbox\DiD Testing Data\onepercentsample.dta" , clear
	encode broadethnic, gen(g_ethnicity)
	encode voters_gender, gen(g_gender)
	encode parties_des , gen(g_party)
	egen id = group(lal)
	drop broadethnic voters_gender parties_des lalvoterid
	gen fips = substr(geoid10, 1, 2)
	destring fips, replace
	drop geoid10
	sort fips year
	save `ops', replace
restore 

merge 1:m fips year using `ops'
keep if _m==3 
drop _m 

drop if state == "South Carolina"
xtset year id 


csdid general,  ivar(id)  time(year) gvar(first_treat_strict) cluster(state_float)
estat pretrend/* can reject null of pre-treatment = 0 */
estat simple 
local csdid_coe = round(r(b)[1,1],0.001)
local csdid_se = round(r(table)[2,1],0.001)

did_multiplegt general id year stricty, placebo(1) breps(50) cluster(id)
local chaise_coe = round(e(effect_0),0.001)
local chaise_se = round(e(se_effect_0),0.001)

did_imputation general id year first_treat_strict, fe(id) autosample
qui eststo reg2  
qui estout reg2 , cells(b(star fmt(5)) se(par fmt(5)))
local boy_coe = round(r(coefs)[1,1],0.001)
local boy_se = round(r(coefs)[1,2],0.001)

did2s general, first_stage(id year) second_stage(first_treat_strict) treatment(first_treat_strict) cluster(id)
eststo reg4  
estout reg4 , cells(b(star fmt(5)) se(par fmt(5)))
local gar_coe = round(r(coefs)[1,1],0.001)
local gar_se = round(r(coefs)[1,2],0.001)

"Estimator: CSDID CHAISE BORUSYAK GARDNER " _n ///
"Beta: `csdid_coe' `chaise_coe' `boy_coe' `gar_coe' " _n ///
"SE: `csdid_se' `chaise_se' `boy_se' `gar_se' " 

/*
merge 1:m year state using "cumulative_2006-2020.dta"
tab year if _m==2
keep if _m==3
drop _m 

*xtset state_float year 
*csdid VEPhighestoffice stricty,  ivar(state_float)  time(year) gvar(first_treat_strict)

bysort state year race gender educ: gen dup = cond(_N==1,0,_n)
summ dup 
drop dup 
*/


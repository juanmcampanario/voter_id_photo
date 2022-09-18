local user = c(username)
global user `user'
cd "C:/Users/$user/Dropbox/my_did_testing/"

capture mkdir juan
capture mkdir juan/tables
capture mkdir juan/graphs 
capture mkdir juan/dtafiles
capture mkdir juan/graphs/throw
capture mkdir juan/graphs/onepercent
capture mkdir juan/graphs/onepercent/throw

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

gen black = (g_eth==4) & g_eth!=. 
gen euro = (g_eth==2) & g_eth!=. 
gen asian = (g_eth==1) & g_eth!=. 
gen hisp = (g_eth==3) & g_eth!=. 

gen male =(g_gend==2) & g_gend!=.

gen dem = (inlist(g_party,9,18)) & g_party!=.
gen rep = (inlist(g_party,20,29)) & g_party!=.

collapse (mean) black euro asian hisp male dem rep general_ stricty event_stricty first_treat_strict state_float tst, by(state year)

gen treat_states =""
levelsof state, local(levels) 
foreach x of local levels {
	qui summ year if state=="`x'" & tst>0
	qui replace treat_states="`x'" if state=="`x'" & tst>0 & r(N)!=0
}

save "juan/dtafiles/onepercent_collapse.dta", replace

use "juan/dtafiles/onepercent_collapse.dta", clear
levelsof treat_states, local(levels)
foreach x of local levels {
	use "juan/dtafiles/onepercent_collapse.dta", clear
	drop if state=="South Carolina"

	xtset state_float year 
	display "`x'"
	keep if state=="`x'" | event_stricty==.
	qui summ state_float if  state=="`x'"
	local trunit = r(min)
	
	qui summ first_treat_strict if state=="`x'"
	local trperiod = r(mean)
	display `trperiod'
	
	local lno_treat = `trperiod'-2
	global controls ""
	forvalues y =2006(2)`lno_treat' {
		global controls $controls general_(`y') black(`y') euro(`y') dem(`y') rep(`y')
	} 
	display "$controls"
	display ""
	
	capture synth general_ $controls , trunit(`trunit') trperiod(`trperiod') nested allopt keep("juan/graphs/onepercent/throw/`x'") replace 
	if !_rc {
		preserve 
			use "juan/graphs/onepercent/throw/`x'.dta", clear
			twoway line _Y_treat _time, lp(solid) lc(black) || line _Y_synth _time, lp(_) lc(gray) ||, ///
			legend(order(1 "`x'" 2 "Synthetic `X'")) xlabel(2005(2)2019) yti("Turnout") xti("Year") ///
			xline(`trperiod', lp(_)) ti("`x'") ///
			graphregion(color(white)) bgcolor(white)
			gr export "juan/graphs/onepercent/`x'_synth.png", replace 
		restore
		
		global gr_gaps ""
		levelsof state_float, local(statelist)
		foreach i of local statelist {
			capture qui synth general_ $controls , trunit(`i') trperiod(`trperiod') keep("juan/graphs/onepercent/throw/state`i'") replace 
			if !_rc{
				preserve
					use "juan/graphs/onepercent/throw/state`i'.dta" ,clear
					keep _Y_treated _Y_synthetic _time 
					drop if _time==. 
					rename _time year 
					rename _Y_treated treat`i' 
					rename _Y_synthetic counterfact`i'
					gen gap`i'=treat`i'-counterfact`i' 
					sort year
					sleep 1000
					save "juan/graphs/onepercent/throw/state`i'.dta", replace 
				restore
				global gr_gaps $gr_gaps (line gap`i' year ,lp(solid)lw(vthin)lcolor(gray)) || 
		}
		}
		 
		use "juan/graphs/onepercent/throw/state`trunit'.dta", clear
		sort year
		save "juan/graphs/onepercent/throw/placebo_`trunit'.dta", replace

		foreach i of local statelist {
			merge year using "juan/graphs/onepercent/throw/state`i'.dta" 
			drop _merge 
			sort year 
			sleep 1000
			save "juan/graphs/onepercent/throw/placebo_`trunit'.dta", replace
		 }
		
		use "juan/graphs/onepercent/throw/placebo_`trunit'.dta", clear
		
		twoway $gr_gaps ///
		(line gap`trunit' year ,lp(solid)lw(thick)lcolor(black)) ||, ///
		yline(0, lpattern(shortdash) lcolor(black)) xline(`trperiod', lpattern(shortdash) ///
		lcolor(black)) xtitle("",si(small)) xlabel(#10) ///
		ytitle("Gap in Turnout prediction error", size(small)) subti("`x'") /// 
		legend(off)
		gr export "juan/graphs/onepercent/`x'_placebo.png", replace 
	}
}

**********************NON-PRES**************************************************

use "juan/dtafiles/onepercent_collapse.dta", clear
levelsof treat_states, local(levels)
foreach x of local levels {
	use "juan/dtafiles/onepercent_collapse.dta", clear
	drop if state=="South Carolina"
	
	xtset state_float year 
	display "`x'"
	keep if state=="`x'" | event_stricty==.
	qui summ state_float if  state=="`x'"
	local trunit = r(min)
	
	qui summ first_treat_strict if state=="`x'"
	local trperiod = r(mean)
	display `trperiod'
	
	local real_xline = `trperiod'
	
	if `trperiod'==2008 | `trperiod'==2012 | `trperiod'==2016 {
		local trperiod = `trperiod'+2 
	} 
	
	drop if inlist(year,2008,2012,2016)

	local lno_treat = `trperiod'-4
	global controls ""
	forvalues y =2006(4)`lno_treat' {
		global controls $controls general_(`y') black(`y') euro(`y') dem(`y') rep(`y')
	} 
	display "$controls"
	display ""
	
	capture synth general_ $controls , trunit(`trunit') trperiod(`trperiod') nested allopt keep("juan/graphs/onepercent/throw/`x'_nonpres") replace 
	if !_rc {
		preserve 
			use "juan/graphs/onepercent/throw/`x'_nonpres.dta", clear
			twoway line _Y_treat _time, lp(solid) lc(black) || line _Y_synth _time, lp(_) lc(gray) ||, ///
			legend(order(1 "`x'" 2 "Synthetic `X'")) xlabel(2005(2)2019) yti("Turnout") xti("Year") ///
			xline(`real_xline', lp(_)) ti("`x'") ///
			graphregion(color(white)) bgcolor(white)
			gr export "juan/graphs/onepercent/`x'_synth_nonpres.png", replace 
		restore
		
		global gr_gaps ""
		levelsof state_float, local(statelist)
		foreach i of local statelist {
			capture qui synth general_ $controls , trunit(`i') trperiod(`trperiod') keep("juan/graphs/onepercent/throw/state`i'_nonpres") replace 
			if !_rc{
				preserve
					use "juan/graphs/onepercent/throw/state`i'_nonpres.dta" ,clear
					keep _Y_treated _Y_synthetic _time 
					drop if _time==. 
					rename _time year 
					rename _Y_treated treat`i' 
					rename _Y_synthetic counterfact`i'
					gen gap`i'=treat`i'-counterfact`i' 
					sort year
					sleep 1000
					save "juan/graphs/onepercent/throw/state`i'_nonpres.dta", replace 
				restore
				global gr_gaps $gr_gaps (line gap`i' year ,lp(solid)lw(vthin)lcolor(gray)) || 
		}
		}
		 
		use "juan/graphs/onepercent/throw/state`trunit'_nonpres.dta", clear
		sort year
		save "juan/graphs/onepercent/throw/placebo_`trunit'_nonpres.dta", replace

		foreach i of local statelist {
			merge year using "juan/graphs/onepercent/throw/state`i'_nonpres.dta" 
			drop _merge 
			sort year 
			sleep 1000
			save "juan/graphs/onepercent/throw/placebo_`trunit'_nonpres.dta", replace
		 }
		
		use "juan/graphs/onepercent/throw/placebo_`trunit'_nonpres.dta", clear
		
		twoway $gr_gaps ///
		(line gap`trunit' year ,lp(solid)lw(thick)lcolor(black)) ||, ///
		yline(0, lpattern(shortdash) lcolor(black)) xline(`real_xline', lpattern(shortdash) ///
		lcolor(black)) xtitle("",si(small)) xlabel(#10) ///
		ytitle("Gap in Turnout prediction error", size(small)) subti("`x'") /// 
		legend(off)
		gr export "juan/graphs/onepercent/`x'_placebo_nonpres.png", replace 
	}
}

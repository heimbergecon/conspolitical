// FIGURE 3

// Response variables
local vars approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary 

* panel titles / y-axis labels
local varsnames " "Approval" "Government crises" "General strikes" "Demonstrations" "
local labels   " "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points""

* controls (one lag each): RGROWTH, RYIELD, REER, and the respective endogenous variable
local cntrls_approval L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval
local cntrls_GOVCRISIS_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  
local cntrls_STRIKE_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  
local cntrls_DEMONSTR_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  

		
// IVREG (binary state: plot according to 30th percentile) ***********************************************************************

cap drop recession
generate recession = .
sum OGAP, detail
	cap drop pOGAP
	_pctile OGAP, p(30)
    gen pOGAP = r(r1)
replace recession = 1 if OGAP < pOGAP & !missing(OGAP)
replace recession = 0 if OGAP > pOGAP & !missing(OGAP)
tab recession

* Endogenous regressors
capture drop diff_STRUCBAL_E diff_STRUCBAL_R
gen diff_STRUCBAL_R = recession * diff_STRUCBAL
gen diff_STRUCBAL_E = (1 - recession) * diff_STRUCBAL

* Instruments
capture drop TOTAL_R TOTAL_E
gen TOTAL_R = recession * TOTAL
gen TOTAL_E = (1 - recession) * TOTAL

xtset country_id year

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 0     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock_R diff_STRUCBAL_R
local shock_E diff_STRUCBAL_E

* options
global savefigs 1    
global verb qui     
cap gen t = _n 
cap gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
	
	* set controls
	local cntrls `cntrls_`var''
		
	* preallocate 
	qui gen biv`var'_R = .
	qui gen up90biv`var'_R = .
	qui gen lo90biv`var'_R = .
	qui gen up68biv`var'_R = .
	qui gen lo68biv`var'_R = .
	qui gen biv`var'_E = .
	qui gen up90biv`var'_E = .
	qui gen lo90biv`var'_E = .
	qui gen up68biv`var'_E = .
	qui gen lo68biv`var'_E = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		cap gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E)   `cntrls' i.year i.country_id, dkraay(2)
		}
		else if `estdiff' == 1 {
			`verb' ivreg2 F(`i').d`var' (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E)   `cntrls' i.year i.country_id,  dkraay(2)
		}
		else if `estdiff' == 2 {
			`verb' ivreg2 d`iname'`var' (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E)    `cntrls'  i.year i.country_id,  dkraay(1) partial(i.year i.country_id)
		}
		
		cap gen biv`var'h`iname'_R = _b[`shock_R']
		cap gen seiv`var'h`iname'_R = _se[`shock_R']
		qui replace biv`var'_R = biv`var'h`iname'_R if h==`i'
		qui replace up90biv`var'_R = biv`var'h`iname'_R + `z1'*seiv`var'h`iname'_R if h==`i'
		qui replace lo90biv`var'_R = biv`var'h`iname'_R - `z1'*seiv`var'h`iname'_R if h==`i'
		qui replace up68biv`var'_R = biv`var'h`iname'_R + `z2'*seiv`var'h`iname'_R if h==`i'
		qui replace lo68biv`var'_R = biv`var'h`iname'_R - `z2'*seiv`var'h`iname'_R if h==`i'
		
		cap gen biv`var'h`iname'_E = _b[`shock_E']
		cap gen seiv`var'h`iname'_E = _se[`shock_E']
		qui replace biv`var'_E = biv`var'h`iname'_E if h==`i'
		qui replace up90biv`var'_E = biv`var'h`iname'_E + `z1'*seiv`var'h`iname'_E if h==`i'
		qui replace lo90biv`var'_E = biv`var'h`iname'_E - `z1'*seiv`var'h`iname'_E if h==`i'
		qui replace up68biv`var'_E = biv`var'h`iname'_E + `z2'*seiv`var'h`iname'_E if h==`i'
		qui replace lo68biv`var'_E = biv`var'h`iname'_E - `z2'*seiv`var'h`iname'_E if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90biv`var'_E lo90biv`var'_E h, fcolor("$mblue%15")  lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68biv`var'_E lo68biv`var'_E h, fcolor("$mblue%40")  lcolor(mdred) lw(none) lpattern(solid)) ///
		(line biv`var'_E h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		 (rarea up90biv`var'_R lo90biv`var'_R h, fcolor("$mgreen%15") lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68biv`var'_R lo68biv`var'_R h, fcolor("$mgreen%40") lcolor(mdred) lw(none) lpattern(solid)) ///
		(line biv`var'_R h,  lcolor("$mgreen") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc("$mred") clw(vvthin)) if h<=`horizon', ///
		title("{bf:`varname'}", size(10) col(black) margin(b=2)) xtitle("Years", margin(t=3) size(7)) ///
		ytitle("`labname'", size(7))   ///
		xlabel(0(1)`horizon', labsize(7)) ///
				ylabel(,labsize(7)) /// 
		plotregion(color(white)) ///
		graphregion(color(white) )    name("irfgs_`var'_ts", replace) ///
						legend(order(3 6) label(3 "Upper regime") label(6 "Lower regime") size(7)region(lcolor(none))) ///
				saving("$FIGUREDIR/pirf_iv_binaryregimes_p30_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_binaryregimes_p30_`var'.pdf",  replace

			local ivar = `ivar'+1

		}		
		
	
	
foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}
	
// IVREG (binary state: splot according to 33th percentile) ***********************************************************************

cap drop recession
generate recession = .
sum OGAP, detail
	cap drop pOGAP
	_pctile OGAP, p(33)
    gen pOGAP = r(r1)
replace recession = 1 if OGAP < pOGAP & !missing(OGAP)
replace recession = 0 if OGAP > pOGAP & !missing(OGAP)
tab recession

* Endogenous regressors
capture drop diff_STRUCBAL_E diff_STRUCBAL_R
gen diff_STRUCBAL_R = recession * diff_STRUCBAL
gen diff_STRUCBAL_E = (1 - recession) * diff_STRUCBAL

* Instruments
capture drop TOTAL_R TOTAL_E
gen TOTAL_R = recession * TOTAL
gen TOTAL_E = (1 - recession) * TOTAL

xtset country_id year

* specs for local projection

local horizon = 5    // Impulse horizon
local estdiff = 0     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock_R diff_STRUCBAL_R
local shock_E diff_STRUCBAL_E


* options
global savefigs 1   
global verb qui     

cap gen t = _n 
cap gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
	
	* set controls
	local cntrls `cntrls_`var''
		
	* preallocate 
	qui gen biv`var'_R = .
	qui gen up90biv`var'_R = .
	qui gen lo90biv`var'_R = .
	qui gen up68biv`var'_R = .
	qui gen lo68biv`var'_R = .
	qui gen biv`var'_E = .
	qui gen up90biv`var'_E = .
	qui gen lo90biv`var'_E = .
	qui gen up68biv`var'_E = .
	qui gen lo68biv`var'_E = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		cap gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E)   `cntrls' i.year i.country_id, dkraay(2)
		}
		else if `estdiff' == 1 {
			`verb' ivreg2 F(`i').d`var' (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E)   `cntrls' i.year i.country_id,  dkraay(2)
		}
		else if `estdiff' == 2 {
			`verb' ivreg2 d`iname'`var' (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E)    `cntrls'  i.year i.country_id,  dkraay(1) partial(i.year i.country_id)
		}
		
		cap gen biv`var'h`iname'_R = _b[`shock_R']
		cap gen seiv`var'h`iname'_R = _se[`shock_R']
		qui replace biv`var'_R = biv`var'h`iname'_R if h==`i'
		qui replace up90biv`var'_R = biv`var'h`iname'_R + `z1'*seiv`var'h`iname'_R if h==`i'
		qui replace lo90biv`var'_R = biv`var'h`iname'_R - `z1'*seiv`var'h`iname'_R if h==`i'
		qui replace up68biv`var'_R = biv`var'h`iname'_R + `z2'*seiv`var'h`iname'_R if h==`i'
		qui replace lo68biv`var'_R = biv`var'h`iname'_R - `z2'*seiv`var'h`iname'_R if h==`i'
		
		cap gen biv`var'h`iname'_E = _b[`shock_E']
		cap gen seiv`var'h`iname'_E = _se[`shock_E']
		qui replace biv`var'_E = biv`var'h`iname'_E if h==`i'
		qui replace up90biv`var'_E = biv`var'h`iname'_E + `z1'*seiv`var'h`iname'_E if h==`i'
		qui replace lo90biv`var'_E = biv`var'h`iname'_E - `z1'*seiv`var'h`iname'_E if h==`i'
		qui replace up68biv`var'_E = biv`var'h`iname'_E + `z2'*seiv`var'h`iname'_E if h==`i'
		qui replace lo68biv`var'_E = biv`var'h`iname'_E - `z2'*seiv`var'h`iname'_E if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90biv`var'_E lo90biv`var'_E h, fcolor("$mblue%15")  lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68biv`var'_E lo68biv`var'_E h, fcolor("$mblue%40")  lcolor(mdred) lw(none) lpattern(solid)) ///
		(line biv`var'_E h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		 (rarea up90biv`var'_R lo90biv`var'_R h, fcolor("$mgreen%15") lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68biv`var'_R lo68biv`var'_R h, fcolor("$mgreen%40") lcolor(mdred) lw(none) lpattern(solid)) ///
		(line biv`var'_R h,  lcolor("$mgreen") lpattern(dash) lwidth(thick)) ///
		(line zero h, lcolor("$mred")   lpattern(solid) lwidth(medthick)) if h<=`horizon', ///
		title("{bf:`varname'}", size(10) col(black) margin(b=2)) xtitle("Years", margin(t=3) size(7)) ///
		ytitle("`labname'", size(7))   ///
		xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
		plotregion(color(white)) ///
		graphregion(color(white) )    name("irfgs_`var'_ts", replace) ///
						legend(order(3 6) label(3 "Upper regime") label(6 "Lower regime") size(7) region(lcolor(none))) ///
				saving("$FIGUREDIR/pirf_iv_binaryregimes_p33_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_binaryregimes_p33_`var'.pdf",  replace

			local ivar = `ivar'+1

		}	
		
		
foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}		
		


////////////////////////////////////////////////////////////////////////////////
// Figure combined
////////////////////////////////////////////////////////////////////////////////

*-----------------------------*
* p33
*-----------------------------*
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p33_approval.gph",            name(g1, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p33_GOVCRISIS_binary.gph",    name(g2, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p33_STRIKE_binary.gph",       name(g3, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p33_DEMONSTR_binary.gph",     name(g4, replace)

grc1leg2 g1 g2 g3 g4, cols(4) ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Recession.pdf", replace

*-----------------------------*
* p30
*-----------------------------*
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p30_approval.gph",            name(g1, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p30_GOVCRISIS_binary.gph",    name(g2, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p30_STRIKE_binary.gph",       name(g3, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p30_DEMONSTR_binary.gph",     name(g4, replace)

grc1leg2 g1 g2 g3 g4, cols(4) ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Recession_Robustness.pdf", replace

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


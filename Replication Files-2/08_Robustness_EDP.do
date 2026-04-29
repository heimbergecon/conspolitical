// Robustness

preserve

* Subset to EU sample
keep if ccode=="AUT" | ccode=="BEL" | ccode=="DEU" | ccode=="DNK" | ccode=="ESP" | ccode=="FIN" | ccode=="FRA" | ccode=="IRL" | ccode=="ITA" | ccode=="NLD" | ccode=="PRT" | ccode=="SWE"


// ---------------------------------------------------------------------------
// EDP DATA
merge m:1 ccode year using "$DATA\EDP_long.dta" 
drop _merge

gen EDP_TOTAL = 0
replace EDP_TOTAL = TOTAL if !missing(EDP_dummy) & EDP_dummy == 1

// LP SETUP
local vars approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary
local varsnames `" "Approval" "Government crises" "General strikes" "Demonstrations" "'
local labels    `" "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points" "'

local shock diff_STRUCBAL

sort country_id year
xtset country_id year, yearly

cap drop t h
gen t = _n
gen h = t-1
local horizon = 4

local CI1 = 0.10
local CI2 = 0.32
local z1  = abs(invnormal(`CI1'/2))
local z2  = abs(invnormal(`CI2'/2))


local cntrls_approval           L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval
local cntrls_GOVCRISIS_binary   L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_STRIKE_binary      L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_DEMONSTR_binary    L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER

local ivar = 1
foreach var in `vars' {

    local cntrls `cntrls_`var''

    cap drop biv`var' up90biv`var' lo90biv`var' up68biv`var' lo68biv`var'
    qui gen biv`var'     = .
    qui gen up90biv`var' = .
    qui gen lo90biv`var' = .
    qui gen up68biv`var' = .
    qui gen lo68biv`var' = .

    forvalues i=0/`horizon' {

        local iname `i'
        ivreg2 F(`i').`var' (diff_STRUCBAL = EDP_TOTAL) `cntrls' i.year i.country_id, dkraay(1)

        cap drop biv`var'h`iname' seiv`var'h`iname'
        gen biv`var'h`iname'  = _b[`shock']
        gen seiv`var'h`iname' = _se[`shock']

        qui replace biv`var'     = biv`var'h`iname' if h==`i'
        qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
        qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
        qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
        qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
    }

    cap drop zero
    gen zero=0

    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'

    tw ///
        (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none)) ///
        (rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none)) ///
        (line  biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
        (line  zero h,     lcolor("$mred")  lpattern(solid) lwidth(medthick)) ///
        if h<=`horizon', ///
        title("{bf:`varname'}", size(10) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`labname'", size(7)) ///
        xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white)) graphregion(color(white)) ///
        legend(off) ///
        name("ivU`ivar'", replace) ///
        saving("$FIGUREDIR/pirf_iv_`var'_EDP.gph", replace)

    graph export "$FIGUREDIR/pirf_iv_`var'_EDP.pdf", fontface($grfont) replace
    local ivar = `ivar'+1
}


graph use "$FIGUREDIR/pirf_iv_approval_EDP.gph",         name(g1, replace)
graph use "$FIGUREDIR/pirf_iv_GOVCRISIS_binary_EDP.gph", name(g2, replace)
graph use "$FIGUREDIR/pirf_iv_STRIKE_binary_EDP.gph",    name(g3, replace)
graph use "$FIGUREDIR/pirf_iv_DEMONSTR_binary_EDP.gph",  name(g4, replace)

grc1leg2  g1 g2 g3 g4, ///
	cols(4) ///
    loff ///
	ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Robustness_EDP.pdf", replace



// ---------------------------------------------------------------------------

// Fiscal balance
merge m:1 ccode year using "$DATA\FiscalBalance_IMF_long.dta" 
drop if _merge == 2
drop _merge

* Raw continuous distance
gen Distance = FiscalBalance + 3

* Lagged distance
bysort country (year): gen Distance_lag = Distance[_n-1]

* Binary breach indicator
gen EDP_breach = (FiscalBalance < -3)
gen TOTAL_EDP_breach = .
replace TOTAL_EDP_breach = TOTAL * EDP_breach if EDP_breach == 1 & !missing(EDP_breach)
replace TOTAL_EDP_breach = 0 if EDP_breach == 0 & !missing(EDP_breach)

* Main interaction instrument
gen Total_Distance = TOTAL * Distance_lag

* Asymmetric versions
gen Distance_below = Distance_lag if FiscalBalance < -3
replace Distance_below = 0 if FiscalBalance >= -3

gen Distance_above = Distance_lag if FiscalBalance >= -3
replace Distance_above = 0 if FiscalBalance < -3


// LP SETUP
local vars approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary
local varsnames `" "Approval" "Government crises" "General strikes" "Demonstrations" "'
local labels    `" "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points" "'

local shock diff_STRUCBAL

sort country_id year
xtset country_id year, yearly

cap drop t h
gen t = _n
gen h = t-1
local horizon = 4

local CI1 = 0.10
local CI2 = 0.32
local z1  = abs(invnormal(`CI1'/2))
local z2  = abs(invnormal(`CI2'/2))


local cntrls_approval           L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval
local cntrls_GOVCRISIS_binary   L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_STRIKE_binary      L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_DEMONSTR_binary    L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER

local ivar = 1
foreach var in `vars' {

    local cntrls `cntrls_`var''

    cap drop biv`var' up90biv`var' lo90biv`var' up68biv`var' lo68biv`var'
    qui gen biv`var'     = .
    qui gen up90biv`var' = .
    qui gen lo90biv`var' = .
    qui gen up68biv`var' = .
    qui gen lo68biv`var' = .

    forvalues i=0/`horizon' {

        local iname `i'
        ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL_EDP_breach) `cntrls' i.year i.country_id, dkraay(1)

        cap drop biv`var'h`iname' seiv`var'h`iname'
        gen biv`var'h`iname'  = _b[`shock']
        gen seiv`var'h`iname' = _se[`shock']

        qui replace biv`var'     = biv`var'h`iname' if h==`i'
        qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
        qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
        qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
        qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
    }

    cap drop zero
    gen zero=0

    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'

    tw ///
        (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none)) ///
        (rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none)) ///
        (line  biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
        (line  zero h,     lcolor("$mred")  lpattern(solid) lwidth(medthick)) ///
        if h<=`horizon', ///
        title("{bf:`varname'}", size(10) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`labname'", size(7)) ///
        xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white)) graphregion(color(white)) ///
        legend(off) ///
        name("ivU`ivar'", replace) ///
        saving("$FIGUREDIR/pirf_iv_`var'_EDP_breach.gph", replace)

    graph export "$FIGUREDIR/pirf_iv_`var'_EDP_breach.pdf", fontface($grfont) replace
    local ivar = `ivar'+1
}


graph use "$FIGUREDIR/pirf_iv_approval_EDP_breach.gph",         name(g1, replace)
graph use "$FIGUREDIR/pirf_iv_GOVCRISIS_binary_EDP_breach.gph", name(g2, replace)
graph use "$FIGUREDIR/pirf_iv_STRIKE_binary_EDP_breach.gph",    name(g3, replace)
graph use "$FIGUREDIR/pirf_iv_DEMONSTR_binary_EDP_breach.gph",  name(g4, replace)

grc1leg2  g1 g2 g3 g4, ///
	cols(4) ///
    loff ///
	ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Robustness_EDP_breach.pdf", replace

restore

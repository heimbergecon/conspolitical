// Macroeconomic effects

* Panel id/time
cap confirm variable country_id
if _rc {
    encode ccode, gen(country_id)
}
sort country_id year
xtset country_id year, yearly

* Transform real GDP as in your templates: 100*log_REALGDP
cap confirm variable log_REALGDP
if _rc {
    di as err "ERROR: log_REALGDP not found in dataset."
    exit 111
}
capture drop temp
gen double temp = 100*log_REALGDP
drop log_REALGDP
rename temp log_REALGDP


********************************************************************************
* RESPONSE VARIABLES + LABELS (ONLY TWO)
********************************************************************************
local vars log_REALGDP unemp
local varsnames `" "Real GDP" "Unemployment rate" "'
local labels    `" "in %" "in %-points" "'

********************************************************************************
* CONTROLS (one lag each): RGROWTH, RYIELD, REER, plus lag of y
* PLUS: lag of change in y (L1.D.y)
********************************************************************************
local cntrls_log_REALGDP  L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).log_REALGDP  L(1/1).D.log_REALGDP
local cntrls_unemp        L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).unemp        L(1/1).D.unemp

********************************************************************************
* IV LOCAL PROJECTIONS — LEVEL (2 outcomes)
********************************************************************************
local horizon  = 4
local estdiff  = 2          // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10
local CI2 = 0.32
local z1  = abs(invnormal(`CI1'/2))
local z2  = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL

global savefigs 1
global verb qui

* horizons index
cap drop t h
cap gen t = _n
cap gen h = t - 1

local ivar = 1
foreach var in `vars' {

    local cntrls `cntrls_`var''

    * preallocate storage for IRFs
    cap drop biv`var' up90biv`var' lo90biv`var' up68biv`var' lo68biv`var'
    qui gen double biv`var'      = .
    qui gen double up90biv`var'  = .
    qui gen double lo90biv`var'  = .
    qui gen double up68biv`var'  = .
    qui gen double lo68biv`var'  = .

    * create first-difference if needed (only used if estdiff==1)
    if `estdiff' > 0 {
        cap drop d`var'
        gen double d`var' = `var' - L.`var'
    }

    forvalues i = 0/`horizon' {

        local iname `i'

        * cumulative change helper (created regardless; used if estdiff==2)
        cap drop d`iname'`var'
        gen double d`iname'`var' = F`i'.`var' - L.`var'

        if `estdiff' == 0 {
            $verb ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)
        }
        else if `estdiff' == 1 {
            $verb ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)
        }
        else if `estdiff' == 2 {
            $verb ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)
        }

        * store coefficient + SE for shock
        cap drop biv`var'h`iname' seiv`var'h`iname'
        gen double biv`var'h`iname'  = _b[`shock']
        gen double seiv`var'h`iname' = _se[`shock']

        qui replace biv`var'      = biv`var'h`iname' if h==`i'
        qui replace up90biv`var'  = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
        qui replace lo90biv`var'  = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
        qui replace up68biv`var'  = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
        qui replace lo68biv`var'  = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
    }

    * plot helpers
    cap drop zero
    gen double zero = 0

    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'

    tw ///
        (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none) lpattern(solid)) ///
        (rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none) lpattern(solid)) ///
        (line  biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
        (line  zero h, lcolor("$mred")  lpattern(solid) lwidth(medthick)) ///
        if h<=`horizon', ///
        title("{bf:`varname'}", size(10) col(black) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`labname'", size(7)) ///
        xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white)) ///
        graphregion(color(white)) ///
        legend(off) ///
        name("iv_`var'", replace) ///
        saving("$FIGUREDIR/pirf_ivonly_`var'.gph", replace)

    graph export "$FIGUREDIR/pirf_ivonly_`var'.pdf", fontface($grfont) replace
    graph export "$FIGUREDIR/pirf_ivonly_`var'.jpg", width(4000) height(3000) replace

    local ivar = `ivar' + 1
}

********************************************************************************
* BASELINE
********************************************************************************
cap graph close _all

graph use "$FIGUREDIR/pirf_ivonly_log_REALGDP.gph", name(g1, replace)
graph use "$FIGUREDIR/pirf_ivonly_unemp.gph",       name(g2, replace)
grc1leg2 g1 g2, cols(2) ///
    loff ///
	ysize(3) xsize(6)
graph export "$FIGURECOM/PolEffects_MacroVars.pdf", replace


********************************************************************************
* State dependent effects
********************************************************************************

* create or reset recession safely
cap confirm variable recession
if _rc {
    gen byte recession = .
}
else {
    replace recession = .
}

* p33 cutoff stored as scalar
capture scalar drop OGAP_p33
_pctile OGAP, p(33)
scalar OGAP_p33 = r(r1)

replace recession = 1 if OGAP <  OGAP_p33 & !missing(OGAP)   // lower regime
replace recession = 0 if OGAP >= OGAP_p33 & !missing(OGAP)   // upper regime
tab recession

* Endogenous regressors (regime-specific shocks)
capture drop diff_STRUCBAL_E diff_STRUCBAL_R
gen double diff_STRUCBAL_R = recession * diff_STRUCBAL
gen double diff_STRUCBAL_E = (1 - recession) * diff_STRUCBAL

* Instruments (regime-specific)
capture drop TOTAL_R TOTAL_E
gen double TOTAL_R = recession * TOTAL
gen double TOTAL_E = (1 - recession) * TOTAL

********************************************************************************
* LP-IV SPECIFICATIONS
********************************************************************************
local horizon = 5      // Impulse horizon
local estdiff  = 2     // 2: cumulative

local CI1 = 0.10
local CI2 = 0.32
local z1  = abs(invnormal(`CI1'/2))
local z2  = abs(invnormal(`CI2'/2))

local shock_R diff_STRUCBAL_R
local shock_E diff_STRUCBAL_E

global savefigs 1
global verb qui        // leave empty to show regression output

* horizons index (safe reset)
cap drop t h
gen long t = _n
gen long h = t - 1

********************************************************************************
* LP-IV LOOP (TWO OUTCOMES)
********************************************************************************


local ivar = 1
foreach var in `vars' {

    * set controls
    local cntrls `cntrls_`var''

    * preallocate storage (drop first so reruns work)
    cap drop biv`var'_R up90biv`var'_R lo90biv`var'_R up68biv`var'_R lo68biv`var'_R
    cap drop biv`var'_E up90biv`var'_E lo90biv`var'_E up68biv`var'_E lo68biv`var'_E

    gen double biv`var'_R     = .
    gen double up90biv`var'_R = .
    gen double lo90biv`var'_R = .
    gen double up68biv`var'_R = .
    gen double lo68biv`var'_R = .

    gen double biv`var'_E     = .
    gen double up90biv`var'_E = .
    gen double lo90biv`var'_E = .
    gen double up68biv`var'_E = .
    gen double lo68biv`var'_E = .

    * create first-difference if needed (kept for completeness)
    if `estdiff' > 0 {
        cap drop d`var'
        gen double d`var' = `var' - L.`var'
    }

    forvalues i = 0/`horizon' {

        local iname `i'

        * cumulative change from t-1 to t+i
        cap drop d`iname'`var'
        gen double d`iname'`var' = F`i'.`var' - L.`var'

        if `estdiff' == 0 {
            $verb ivreg2 F(`i').`var' ///
                (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E) ///
                `cntrls' i.year i.country_id, dkraay(2)
        }
        else if `estdiff' == 1 {
            $verb ivreg2 F(`i').d`var' ///
                (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E) ///
                `cntrls' i.year i.country_id, dkraay(2)
        }
        else if `estdiff' == 2 {
            $verb ivreg2 d`iname'`var' ///
                (diff_STRUCBAL_R diff_STRUCBAL_E = TOTAL_R TOTAL_E) ///
                `cntrls' i.year i.country_id, ///
                dkraay(1) partial(i.year i.country_id)
        }

        * store coefficients + SEs (lower regime = R)
        cap drop biv`var'h`iname'_R seiv`var'h`iname'_R
        gen double biv`var'h`iname'_R  = _b[`shock_R']
        gen double seiv`var'h`iname'_R = _se[`shock_R']

        replace biv`var'_R     = biv`var'h`iname'_R if h==`i'
        replace up90biv`var'_R = biv`var'h`iname'_R + `z1'*seiv`var'h`iname'_R if h==`i'
        replace lo90biv`var'_R = biv`var'h`iname'_R - `z1'*seiv`var'h`iname'_R if h==`i'
        replace up68biv`var'_R = biv`var'h`iname'_R + `z2'*seiv`var'h`iname'_R if h==`i'
        replace lo68biv`var'_R = biv`var'h`iname'_R - `z2'*seiv`var'h`iname'_R if h==`i'

        * store coefficients + SEs (upper regime = E)
        cap drop biv`var'h`iname'_E seiv`var'h`iname'_E
        gen double biv`var'h`iname'_E  = _b[`shock_E']
        gen double seiv`var'h`iname'_E = _se[`shock_E']

        replace biv`var'_E     = biv`var'h`iname'_E if h==`i'
        replace up90biv`var'_E = biv`var'h`iname'_E + `z1'*seiv`var'h`iname'_E if h==`i'
        replace lo90biv`var'_E = biv`var'h`iname'_E - `z1'*seiv`var'h`iname'_E if h==`i'
        replace up68biv`var'_E = biv`var'h`iname'_E + `z2'*seiv`var'h`iname'_E if h==`i'
        replace lo68biv`var'_E = biv`var'h`iname'_E - `z2'*seiv`var'h`iname'_E if h==`i'
    }

    * plot helpers
    cap drop zero
    gen double zero = 0

    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'

    * Safe twoway
    #delimit ;
    twoway
        (rarea up90biv`var'_E lo90biv`var'_E h, fcolor("$mblue%15")  lwidth(none) lpattern(solid))
        (rarea up68biv`var'_E lo68biv`var'_E h, fcolor("$mblue%40")  lwidth(none) lpattern(solid))
        (line  biv`var'_E h, lcolor("$mblue") lpattern(dash) lwidth(thick))
        (rarea up90biv`var'_R lo90biv`var'_R h, fcolor("$mgreen%15") lwidth(none) lpattern(solid))
        (rarea up68biv`var'_R lo68biv`var'_R h, fcolor("$mgreen%40") lwidth(none) lpattern(solid))
        (line  biv`var'_R h,  lcolor("$mgreen") lpattern(dash) lwidth(thick))
        (line  zero h, lcolor("$mred")   lpattern(solid) lwidth(medthick))
        if h<=`horizon',
       title("{bf:`varname'}", size(10) col(black) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`labname'", size(7)) ///
        xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white))
        graphregion(color(white))
        name("irfgs_`var'_p33", replace)
        legend(order(3 6) label(3 "Upper regime") label(6 "Lower regime") size(7)  region(lcolor(none) fcolor(none)))
        saving("$FIGUREDIR/pirf_iv_binaryregimes_p33_`var'.gph", replace)
    ;
    #delimit cr

    graph export "$FIGUREDIR/pirf_iv_binaryregimes_p33_`var'.jpg", replace

    local ivar = `ivar' + 1
}

********************************************************************************
* State dependent effects
********************************************************************************
cap graph close _all

graph use "$FIGUREDIR/pirf_iv_binaryregimes_p33_log_REALGDP.gph", name(g1, replace)
graph use "$FIGUREDIR/pirf_iv_binaryregimes_p33_unemp.gph",       name(g2, replace)
grc1leg2 g1 g2, cols(2) ///
	ysize(3) xsize(6)
graph export "$FIGURECOM/PolEffects_MacroVars_Recession.pdf", replace


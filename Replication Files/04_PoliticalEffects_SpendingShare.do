// FIGURE 4

version 17.0
set more off

********************************************************************************
* RESPONSE VARIABLES
********************************************************************************

local vars approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary

********************************************************************************
* PANEL TITLES / Y-AXIS LABELS
********************************************************************************

local varsnames `" "Approval" "Government crises" "General strikes" "Demonstrations" "'

local labels `" "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points" "'

********************************************************************************
* CONTROLS
********************************************************************************

local cntrls_approval           L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval
local cntrls_GOVCRISIS_binary   L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_STRIKE_binary      L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_DEMONSTR_binary    L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER

********************************************************************************
* LP-IV
********************************************************************************

xtset country_id year, yearly

* specs for local projection
local horizon = 4

local CI1 = 0.10
local CI2 = 0.32
local z1 = abs(invnormal(`CI1'/2))
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local Z     TOTAL

* evaluation points
scalar s0   = 0
scalar s50  = 0.5
scalar s100 = 1

* options
global verb qui

di as txt "share_spend_cons values used for IRFs:"
di as txt "  s0   = " %9.4f s0   "  (fully tax-based)"
di as txt "  s50  = " %9.4f s50  "  (50/50 mix)"
di as txt "  s100 = " %9.4f s100 "  (fully spending-based)"

* -----------------------------
* Exogenous main effects
* -----------------------------
local het_ctrls share_spend_cons 

* -----------------------------
* Cleanup old graphs (optional)
* -----------------------------
cap erase "$FIGUREDIR/pirf_het3_approval.gph"
cap erase "$FIGUREDIR/pirf_het3_GOVCRISIS_binary.gph"
cap erase "$FIGUREDIR/pirf_het3_STRIKE_binary.gph"
cap erase "$FIGUREDIR/pirf_het3_DEMONSTR_binary.gph"
cap erase "$FIGUREDIR/iv_political_spendshare_0_50_100_4panel.jpg"

* -----------------------------
* Horizon storage dataset
* -----------------------------
tempfile irfstore
preserve
clear
set obs `= `horizon' + 1'
gen h = _n - 1
save `irfstore', replace
restore


////////////////////////////////////////////////////////////////////////////////
// Level change
////////////////////////////////////////////////////////////////////////////////

local ivar = 1
foreach var in `vars' {
    local cntrls `cntrls_`var''
    sort country_id year
    xtset country_id year, yearly
    preserve
    use `irfstore', clear
    gen double b_s0      = .
    gen double up90_s0   = .
    gen double lo90_s0   = .
    gen double up68_s0   = .
    gen double lo68_s0   = .
    gen double b_s50     = .
    gen double up90_s50  = .
    gen double lo90_s50  = .
    gen double up68_s50  = .
    gen double lo68_s50  = .
    gen double b_s100    = .
    gen double up90_s100 = .
    gen double lo90_s100 = .
    gen double up68_s100 = .
    gen double lo68_s100 = .
    gen double p_gamma   = .
    tempfile thisirf
    save `thisirf', replace
    restore

    local df = `horizon' + 1
    matrix b_gamma_`var'  = J(`df', 1, .)
    matrix se_gamma_`var' = J(`df', 1, .)

    forvalues i = 0/`horizon' {
        cap drop d`i'_`var'
        gen double d`i'_`var' = F`i'.`var' - L.`var'
        cap drop shock_x
        gen double shock_x = `shock' * share_spend_cons
        cap drop Z_x
        gen double Z_x = `Z' * share_spend_cons
       
	   qui ivreg2 F`i'.`var' ///
            (`shock' shock_x = `Z' Z_x) ///
            `cntrls' `het_ctrls' i.year i.country_id, dkraay(1)

        * Store γ_h for joint test
        local row = `i' + 1
        scalar tmp_b  = _b[shock_x]
        scalar tmp_se = _se[shock_x]
        matrix b_gamma_`var'[`row', 1]  = tmp_b
        matrix se_gamma_`var'[`row', 1] = tmp_se
        
        * Horizon-by-horizon p-value for γ_h = 0
        test shock_x = 0
        scalar p_g = r(p)

        * Implied IRFs at s0 / s50 / s100
        lincom `shock' + (s0)*shock_x
        scalar b0  = r(estimate)
        scalar se0 = r(se)
        lincom `shock' + (s50)*shock_x
        scalar b50  = r(estimate)
        scalar se50 = r(se)
        lincom `shock' + (s100)*shock_x
        scalar b100  = r(estimate)
        scalar se100 = r(se)

        preserve
        use `thisirf', clear
        replace b_s0      = b0               if h==`i'
        replace up90_s0   = b0 + `z1'*se0   if h==`i'
        replace lo90_s0   = b0 - `z1'*se0   if h==`i'
        replace up68_s0   = b0 + `z2'*se0   if h==`i'
        replace lo68_s0   = b0 - `z2'*se0   if h==`i'
        replace b_s50     = b50              if h==`i'
        replace up90_s50  = b50 + `z1'*se50 if h==`i'
        replace lo90_s50  = b50 - `z1'*se50 if h==`i'
        replace up68_s50  = b50 + `z2'*se50 if h==`i'
        replace lo68_s50  = b50 - `z2'*se50 if h==`i'
        replace b_s100    = b100               if h==`i'
        replace up90_s100 = b100 + `z1'*se100 if h==`i'
        replace lo90_s100 = b100 - `z1'*se100 if h==`i'
        replace up68_s100 = b100 + `z2'*se100 if h==`i'
        replace lo68_s100 = b100 - `z2'*se100 if h==`i'
        replace p_gamma   = scalar(p_g)               if h==`i'
        save `thisirf', replace
        restore
    }

    * Joint Wald Test
    matrix V_gamma = J(`df', `df', 0)
    forvalues i = 0/`horizon' {
        local row = `i' + 1
        scalar tmp_se = se_gamma_`var'[`row', 1]
        matrix V_gamma[`row', `row'] = tmp_se^2
    }
    matrix W_mat   = b_gamma_`var'' * inv(V_gamma) * b_gamma_`var'
    scalar W_stat  = W_mat[1,1]
    scalar p_joint = 1 - chi2(`df', W_stat)

    display " "
    display "============================================================"
    display " Variable: `var'"
    display " Joint Wald Test: H0: gamma_h = 0 for all h = 0,...,`horizon'"
    display " (i.e. spending share does not shift the IRF at any horizon)"
    display " Chi2(`df') = " %7.3f W_stat
    display " p-value    = " %7.4f p_joint
    display " NOTE: Diagonal approx. — cross-horizon covariance ignored"
    display "============================================================"



    preserve
    use `thisirf', clear
	
	 svmat b_gamma_`var',  names(gamma_coef)
    svmat se_gamma_`var', names(gamma_se)
    gen gamma_tstat = gamma_coef / gamma_se

    display " "
    display "============================================================"
    display " Variable: `var' — Horizon-by-horizon gamma_h results"
    display "============================================================"
    list h gamma_coef gamma_se gamma_tstat p_gamma if h <= `horizon', noobs clean
    display "============================================================"

    drop gamma_coef gamma_se gamma_tstat
	
	* Create significance markers sitting on each IRF line
	gen marker_s0   = b_s0   if p_gamma < 0.10 & h <= `horizon'
	gen marker_s50  = b_s50  if p_gamma < 0.10 & h <= `horizon'
	gen marker_s100 = b_s100 if p_gamma < 0.10 & h <= `horizon'
	
    gen zero = 0
    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'
    local pval_fmt : display %4.3f scalar(p_joint)

    twoway ///
        (rarea up90_s0 lo90_s0 h, fcolor("$mgreen%12") lwidth(none) lpattern(solid)) ///
        (rarea up68_s0 lo68_s0 h, fcolor("$mgreen%28") lwidth(none) lpattern(solid)) ///
        (line  b_s0      h, lcolor("$mgreen")  lpattern(solid) lwidth(thick)) ///
        (rarea up90_s50 lo90_s50 h, fcolor("$morange%12") lwidth(none) lpattern(solid)) ///
        (rarea up68_s50 lo68_s50 h, fcolor("$morange%28") lwidth(none) lpattern(solid)) ///
        (line  b_s50     h, lcolor("$morange") lpattern(dash)  lwidth(thick)) ///
        (rarea up90_s100 lo90_s100 h, fcolor("$mblue%12") lwidth(none) lpattern(solid)) ///
        (rarea up68_s100 lo68_s100 h, fcolor("$mblue%28") lwidth(none) lpattern(solid)) ///
        (line  b_s100    h, lcolor("$mblue")   lpattern(shortdash_dot) lwidth(thick)) ///
        (line  zero      h, lcolor("$mred")    lpattern(solid) lwidth(medthick)) ///
        if h<=`horizon', ///
        title("{bf:`varname'}", size(8) col(black) margin(b=2)) ///
        xtitle("Years", size(5)) ///
        ytitle("`labname'", size(5)) ///
        xlabel(0(1)`horizon', labsize(5)) ///
        ylabel(,labsize(5)) ///
        note("Joint Wald test: p = `pval_fmt'", size(5)) ///
        plotregion(color(white)) ///
        graphregion(color(white)) ///
        legend(order(3 "0% spending share" 6 "50% spending share" 9 "100% spending share") size(7) region(lcolor(none))) ///
        name("het3_`var'", replace) ///
        saving("$FIGUREDIR/pirf_het3_`var'.gph", replace)
    restore
    local ivar = `ivar' + 1
}


////////////////////////////////////////////////////////////////////////////////
// Figure combined
////////////////////////////////////////////////////////////////////////////////

graph use "$FIGUREDIR/pirf_het3_approval.gph",         name(g1, replace)
graph use "$FIGUREDIR/pirf_het3_GOVCRISIS_binary.gph", name(g2, replace)
graph use "$FIGUREDIR/pirf_het3_STRIKE_binary.gph",    name(g3, replace)
graph use "$FIGUREDIR/pirf_het3_DEMONSTR_binary.gph",  name(g4, replace)

grc1leg2 g1 g2 g3 g4, cols(4) ysize(3.5) xsize(12) 
graph export "$FIGURECOM/iv_political_spendshare_0_50_100_4panel.pdf", replace

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


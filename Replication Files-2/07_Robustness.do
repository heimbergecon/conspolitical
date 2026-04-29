// Robustness

version 16.0
set more off

preserve 

cap program drop run_iv_political_robust
program define run_iv_political_robust
    version 16.0

    syntax, SPECname(string) ///
        [ADDPDEBT(integer 0) ADDOGAP(integer 0) DROPREER(integer 0) KEEP0020(integer 0)]

    * LP + CI settings
    local horizon = 4
    local CI1 = 0.10
    local CI2 = 0.32
    local z90 = abs(invnormal(`CI1'/2))
    local z68 = abs(invnormal(`CI2'/2))

    * Robustness: restrict sample to 2000–2020 (inclusive)
    if `keep0020'==1 {
        keep if inrange(year, 2000, 2020)
        sort country_id year
        xtset country_id year, yearly
    }

    * Create binary outcomes (in %-points)
    cap drop GOVCRISIS_binary STRIKE_binary DEMONSTR_binary
    gen GOVCRISIS_binary = 0
    replace GOVCRISIS_binary = 100 if GOVCRISIS > 0

    gen STRIKE_binary = 0
    replace STRIKE_binary = 100 if STRIKE > 0

    gen DEMONSTR_binary = 0
    replace DEMONSTR_binary = 100 if DEMONSTR > 0

    * Horizon index variables (storage trick)
    cap drop t h
    gen t = _n
    gen h = t - 1

    * Outcomes + pretty labels
    local vars      approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary
    local varsnames `" "Approval" "Government crises" "General strikes" "Demonstrations" "'
    local labels    `" "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points" "'

    ********************************************************************************
    * Build baseline control block, then apply ONE modification per robustness spec
    ********************************************************************************
    local BASE_COMMON "L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER"
    if `dropreer'==1 {
        local BASE_COMMON : subinstr local BASE_COMMON "L(1/1).REER" "", all
    }

    * outcome-specific controls (baseline structure)
    local C_approval         "`BASE_COMMON' L(1/1).approval"
    local C_GOVCRISIS_binary "`BASE_COMMON'"
    local C_STRIKE_binary    "`BASE_COMMON'"
    local C_DEMONSTR_binary  "`BASE_COMMON'"

    * (1) add PDEBT
    if `addpdebt'==1 {
        cap confirm variable PDEBT
        if _rc {
            di as error "Spec `specname': PDEBT not found in dataset."
            exit 111
        }
        local C_approval         "`C_approval' L(1/1).PDEBT"
        local C_GOVCRISIS_binary "`C_GOVCRISIS_binary' L(1/1).PDEBT"
        local C_STRIKE_binary    "`C_STRIKE_binary' L(1/1).PDEBT"
        local C_DEMONSTR_binary  "`C_DEMONSTR_binary' L(1/1).PDEBT"
    }

    * (2) add output gap
    local ogapvar ""
    if `addogap'==1 {
        foreach cand in OGAP ogap OUTPUTGAP outputgap OUTPUT_GAP output_gap YGAP ygap GAP gap {
            cap confirm variable `cand'
            if !_rc {
                local ogapvar "`cand'"
                continue, break
            }
        }
        if ("`ogapvar'"=="") {
            di as error "Spec `specname': output gap variable not found (tried OGAP/OUTPUTGAP/OUTPUT_GAP/YGAP/GAP variants)."
            exit 111
        }
        local C_approval         "`C_approval' L(1/1).`ogapvar'"
        local C_GOVCRISIS_binary "`C_GOVCRISIS_binary' L(1/1).`ogapvar'"
        local C_STRIKE_binary    "`C_STRIKE_binary' L(1/1).`ogapvar'"
        local C_DEMONSTR_binary  "`C_DEMONSTR_binary' L(1/1).`ogapvar'"
    }

    ********************************************************************************
    * Run LP-IV for each outcome
    * IMPORTANT: variable names are SHORT and do NOT include `specname` (avoids r(198))
    ********************************************************************************
    local ivar = 1
    foreach y in `vars' {

        * short tag for names (<=32 chars)
        local ytag = cond("`y'"=="approval","appr", ///
                     cond("`y'"=="GOVCRISIS_binary","gcr", ///
                     cond("`y'"=="STRIKE_binary","str","dem")))

        local CNTRLS `C_`y''

        * preallocate IRF storage (short names)
        cap drop b_`ytag' up90_`ytag' lo90_`ytag' up68_`ytag' lo68_`ytag'
        qui gen b_`ytag'     = .
        qui gen up90_`ytag'  = .
        qui gen lo90_`ytag'  = .
        qui gen up68_`ytag'  = .
        qui gen lo68_`ytag'  = .

        forvalues i = 0/`horizon' {

            * LEVEL LP: dependent variable = F(i).y
            $verb ivreg2 F(`i').`y' ///
                (diff_STRUCBAL = TOTAL) ///
                `CNTRLS' i.year i.country_id, ///
                dkraay(1)

            * store coefficient + SE for shock
            cap drop b_`ytag'_h`i' se_`ytag'_h`i'
            gen b_`ytag'_h`i'  = _b[diff_STRUCBAL]
            gen se_`ytag'_h`i' = _se[diff_STRUCBAL]

            qui replace b_`ytag'     = b_`ytag'_h`i' if h==`i'
            qui replace up90_`ytag'  = b_`ytag'_h`i' + `z90'*se_`ytag'_h`i' if h==`i'
            qui replace lo90_`ytag'  = b_`ytag'_h`i' - `z90'*se_`ytag'_h`i' if h==`i'
            qui replace up68_`ytag'  = b_`ytag'_h`i' + `z68'*se_`ytag'_h`i' if h==`i'
            qui replace lo68_`ytag'  = b_`ytag'_h`i' - `z68'*se_`ytag'_h`i' if h==`i'
        }

        cap drop zero
        gen zero = 0

        local yttl : word `ivar' of `varsnames'
        local ylab : word `ivar' of `labels'

        * graph name: short + includes specname (fine; graph names can be short)
        local gname = "g_`specname'_`ytag'"

        tw ///
            (rarea up90_`ytag' lo90_`ytag' h, fcolor("$mblue%15") lwidth(none)) ///
            (rarea up68_`ytag' lo68_`ytag' h, fcolor("$mblue%40") lwidth(none)) ///
            (line  b_`ytag' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
            (line  zero h,     lcolor("$mred")  lpattern(solid) lwidth(medthick)) ///
            if h<=`horizon', ///
            title("{bf:`yttl'}" , size(10) col(black) margin(b=2)) ///
            xtitle("Years", size(7)) ///
            ytitle("`ylab'", size(7)) ///
            xlabel(0(1)`horizon', labsize(7)) ///
			ylabel(,labsize(7)) ///
            plotregion(color(white)) graphregion(color(white)) ///
            legend(off) ///
            name("`gname'", replace)

        graph export "$FIGUREDIR/irf_iv_pol_`specname'_`y'.pdf", fontface($grfont) replace
        graph save   "$FIGUREDIR/irf_iv_pol_`specname'_`y'.gph", replace

        local ivar = `ivar' + 1
    }

********************************************************************************
 * Combine plots
********************************************************************************
  
  grc1leg2  g_`specname'_appr g_`specname'_gcr ///
    g_`specname'_str  g_`specname'_dem, ///
	cols(4) ///
    loff ///
	ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Robustness_`specname'.pdf", replace
 
end

* (1) Add PDEBT as control variable
run_iv_political_robust, specname("addPDEBT") addpdebt(1)

* (2) Add output gap as control variable
run_iv_political_robust, specname("addOGAP") addogap(1)

* (3) Drop REER as control variable
run_iv_political_robust, specname("dropREER") dropreer(1)

* (4) Restrict time period to 2000–2020 (inclusive)
run_iv_political_robust, specname("keep2000_2020") keep0020(1)

restore

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

******************************************************************************
* Robustness: EU subsample
******************************************************************************

local vars      approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary
local varsnames `" "Approval" "Government crises" "General strikes" "Demonstrations" "'
local labels    `" "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points" "'

local shock diff_STRUCBAL

preserve
keep if ccode=="AUT" | ccode=="BEL" | ccode=="DEU" | ccode=="DNK" | ccode=="ESP" | ccode=="FIN" | ccode=="FRA" | ccode=="IRL" | ccode=="ITA" | ccode=="NLD" | ccode=="PRT" | ccode=="SWE"

sort country_id year
xtset country_id year, yearly

cap drop t h
gen t=_n
gen h=t-1
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
        ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)

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
        saving("$FIGUREDIR/pirf_iv_`var'_EU.gph", replace)

    graph export "$FIGUREDIR/pirf_iv_`var'_EU.pdf", fontface($grfont) replace
    local ivar = `ivar'+1
}
restore

graph use "$FIGUREDIR/pirf_iv_approval_EU.gph",         name(g1, replace)
graph use "$FIGUREDIR/pirf_iv_GOVCRISIS_binary_EU.gph", name(g2, replace)
graph use "$FIGUREDIR/pirf_iv_STRIKE_binary_EU.gph",    name(g3, replace)
graph use "$FIGUREDIR/pirf_iv_DEMONSTR_binary_EU.gph",  name(g4, replace)

grc1leg2  g1 g2 g3 g4, ///
	cols(4) ///
    loff ///
	ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Robustness_EU.pdf", replace

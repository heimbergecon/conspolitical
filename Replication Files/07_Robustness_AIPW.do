// Robustness (AIPW)

version 16.0
set more off

********************************************************************************
* REQUIRED VARIABLES CHECKS
********************************************************************************
foreach v in TOTAL_binary RGROWTH RYIELD REER debtgdp {
    cap confirm variable `v'
    if _rc {
        di as error "Required variable `v' not found in dataset."
        exit 111
    }
}

********************************************************************************
* HORIZON INDEX VARIABLES (storage trick)
********************************************************************************
cap drop t h
gen t = _n
gen h = t - 1

********************************************************************************
* CI SETTINGS (match baseline IV file)
********************************************************************************
local H = 4
local CI1 = 0.10
local CI2 = 0.32
local z90 = abs(invnormal(`CI1'/2))
local z68 = abs(invnormal(`CI2'/2))

********************************************************************************
* BASELINE CONTROLS (political baseline)
********************************************************************************
local cntrls_approval           L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval
local cntrls_GOVCRISIS_binary   L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_STRIKE_binary      L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER
local cntrls_DEMONSTR_binary    L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER

********************************************************************************
* PROGRAM: AIPW LP FOR ONE OUTCOME (LEVEL LHS = F(i).Y)
********************************************************************************
cap program drop aipw_lp_one_level
program define aipw_lp_one_level
    version 16.0
    syntax, Y(name) CNTRLS(string) ///
        H(integer) Z90(real) Z68(real) ///
        TITLE(string) YTITLE(string) ///
        TAG(string)

    * locals created by syntax are lowercase:
    * `y' `cntrls' `h' `z90' `z68' `title' `ytitle' `tag'

    * --- Propensity score ---
    capture drop pihat pihat0
    xi: probit TOTAL_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L.debtgdp L.TOTAL_binary i.country_id i.year

    * Short probit output filename
    cap which outreg2
    if _rc==0 {
        outreg2 using "$TABLEDIR/probit_AIPW_LVL_`tag'.xls", ///
            excel replace se dec(3) label
    }

    predict pihat0

    * truncate (diagnostics)
    gen pihat = pihat0
    replace pihat = .9 if pihat>.9 & pihat~=.
    replace pihat = .1 if pihat<.1 & pihat~=.

    * treatment + inverse weights
    capture drop a invwt
    gen a = TOTAL_binary
    gen invwt = a/pihat0 + (1-a)/(1-pihat0) if pihat0~=.

    * IRF storage (short names)
    cap drop b_`tag' se_`tag' up90_`tag' lo90_`tag' up68_`tag' lo68_`tag'
    qui gen b_`tag'     = .
    qui gen se_`tag'    = .
    qui gen up90_`tag'  = .
    qui gen lo90_`tag'  = .
    qui gen up68_`tag'  = .
    qui gen lo68_`tag'  = .

    forvalues i = 0/`h' {

        * LEVEL outcome at horizon i
        cap drop y`i'_`tag'
        gen y`i'_`tag' = F`i'.`y'

        reg y`i'_`tag' TOTAL_binary `cntrls' i.country_id i.year ///
            [pweight=invwt], cluster(country_id)

        gen samp = e(sample)

        predict mu0 if samp==1 & TOTAL_binary==0
        predict mu1 if samp==1 & TOTAL_binary==1
        replace mu0 = mu1 - _b[TOTAL_binary] if samp==1 & TOTAL_binary==1
        replace mu1 = mu0 + _b[TOTAL_binary] if samp==1 & TOTAL_binary==0

        * AIPW / DR1
        generate mdiff1 = (-(a-pihat0)*mu1/pihat0) - ((a-pihat0)*mu0/(1-pihat0))
        generate iptw   = (2*a-1)*y`i'_`tag'*invwt
        generate dr1    = iptw + mdiff1

        qui gen ATE = 1
        qui reg dr1 ATE, nocons cluster(country_id)

        cap drop bh`i'_`tag' seh`i'_`tag'
        gen bh`i'_`tag' = _b[ATE]

        sum dr1
        local dr1m = r(mean)
        gen Isq = (dr1-`dr1m')^2
        sum Isq
        gen seh`i'_`tag' = sqrt(r(mean)/r(N))

        qui replace b_`tag'     = bh`i'_`tag' if h==`i'
        qui replace se_`tag'    = seh`i'_`tag' if h==`i'
        qui replace up90_`tag'  = bh`i'_`tag' + `z90'*seh`i'_`tag' if h==`i'
        qui replace lo90_`tag'  = bh`i'_`tag' - `z90'*seh`i'_`tag' if h==`i'
        qui replace up68_`tag'  = bh`i'_`tag' + `z68'*seh`i'_`tag' if h==`i'
        qui replace lo68_`tag'  = bh`i'_`tag' - `z68'*seh`i'_`tag' if h==`i'

        * cleanup
        capture drop iptw Isq mdiff1 dr1 mu1 mu0 samp ATE
        capture scalar drop dr1m
    }

    cap drop zero
    gen zero = 0

    * graph name <=32 chars
    local gname = "gL_`tag'"

    tw ///
        (rarea up90_`tag' lo90_`tag' h, fcolor("$mblue%15") lwidth(none)) ///
        (rarea up68_`tag' lo68_`tag' h, fcolor("$mblue%40") lwidth(none)) ///
        (line  b_`tag' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
        (line  zero h,    lcolor("$mred")  lpattern(solid) lwidth(medthick)) ///
        if h<=`h', ///
        title("{bf:`title'}", size(10) col(black) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`ytitle'", size(7)) ///
        xlabel(0(1)`h', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white)) graphregion(color(white)) ///
        legend(off) ///
        name("`gname'", replace)

    graph export "$FIGUREDIR/irf_AIPW_LVL_`tag'.pdf", fontface($grfont) replace
    graph save   "$FIGUREDIR/irf_AIPW_LVL_`tag'.gph", replace
end

********************************************************************************
* RUN AIPW LEVEL IRFs (4 outcomes)
********************************************************************************
aipw_lp_one_level, ///
    y(approval) cntrls("`cntrls_approval'") ///
    h(`H') z90(`z90') z68(`z68') ///
    title("Approval") ytitle("in %-points") ///
    tag("ap")

aipw_lp_one_level, ///
    y(GOVCRISIS_binary) cntrls("`cntrls_GOVCRISIS_binary'") ///
    h(`H') z90(`z90') z68(`z68') ///
    title("Government crises") ytitle("change in probability in %-points") ///
    tag("gc")

aipw_lp_one_level, ///
    y(STRIKE_binary) cntrls("`cntrls_STRIKE_binary'") ///
    h(`H') z90(`z90') z68(`z68') ///
    title("General strikes") ytitle("change in probability in %-points") ///
    tag("st")

aipw_lp_one_level, ///
    y(DEMONSTR_binary) cntrls("`cntrls_DEMONSTR_binary'") ///
    h(`H') z90(`z90') z68(`z68') ///
    title("Demonstrations") ytitle("change in probability in %-points") ///
    tag("dm")

	
////////////////////////////////////////////////////////////////////////////////
// Figure combined
////////////////////////////////////////////////////////////////////////////////

grc1leg2 ///
    gL_ap gL_gc ///
    gL_st gL_dm, ///
	cols(4) ///
    loff ///
	ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Robustness_AIPW.pdf", replace

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
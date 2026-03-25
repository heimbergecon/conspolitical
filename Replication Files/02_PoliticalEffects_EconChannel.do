// FIGURE 2

// Response variables
local vars approval GOVCRISIS_binary STRIKE_binary DEMONSTR_binary 

* panel titles / y-axis labels
local varsnames " "Approval" "Government crises" "General strikes" "Demonstrations" "
local labels   " "in %-points" "change in probability in %-points" "change in probability in %-points" "change in probability in %-points""

// introduce leads
forvalues h = 0/4 {
    gen cum_rGDP_`h'   = F`h'.rGDP  - rGDP
    gen cum_unemp_`h'  = F`h'.unemp - unemp
}

// IVREG: controlling for leads ***********************************************

* controls (one lag each): RGROWTH, RYIELD, REER, and the respective endogenous variable
local cntrls_approval L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval 
local cntrls_GOVCRISIS_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  
local cntrls_STRIKE_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  
local cntrls_DEMONSTR_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  

* specs for local projection
local p = 0           // p = number of lags of endogenous variable
local ps = 0          // p = number of lags for shock
local horizon = 4     // Impulse horizon
local estdiff = 0      // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32      // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL

* options
global savefigs 1    /* 1: save to disk, 0: don't save */
global verb qui      /* leave empty if you want to display regression results */

* horizons index
cap drop t h
cap gen t = _n
cap gen h = t - 1


local ivar = 1
foreach var in `vars' {

    * set controls
    local cntrls `cntrls_`var''

    * preallocate storage for IRFs
    cap drop biv`var' up90biv`var' lo90biv`var' up68biv`var' lo68biv`var'
    qui gen biv`var'     = .
    qui gen up90biv`var' = .
    qui gen lo90biv`var' = .
    qui gen up68biv`var' = .
    qui gen lo68biv`var' = .

    * create first-difference if needed (not used when estdiff==2 except for estdiff==1 branch)
    if `estdiff' > 0 {
        cap drop d`var'
        gen d`var' = `var' - L.`var'
    }

    forvalues i = 0/`horizon' {

        local iname `i'

        * cumulative change from t-1 to t+i
        cap drop d`iname'`var'
        gen d`iname'`var' = F`i'.`var' - L.`var'

        if `estdiff' == 0 {
            ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL) `cntrls' F(0/`iname').RGROWTH F(0/`iname').D.unemp i.year i.country_id, dkraay(1)
        }
        else if `estdiff' == 1 {
            `verb' ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL) `cntrls'  F(0/`iname').RGROWTH F(0/`iname').D.unemp i.year i.country_id, dkraay(1)
        }
        else if `estdiff' == 2 {
            `verb' ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL) `cntrls'  F(0/`iname').RGROWTH F(0/`iname').D.unemp i.year i.country_id, dkraay(1)
        }

        * store coefficient + SE for shock
        cap drop biv`var'h`iname' seiv`var'h`iname'
        gen biv`var'h`iname'  = _b[`shock']
        gen seiv`var'h`iname' = _se[`shock']

        qui replace biv`var'     = biv`var'h`iname' if h==`i'
        qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
        qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
        qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
        qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
    }

    * plot helpers
    cap drop zero
    gen zero = 0

    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'

    tw ///
        (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none) lpattern(solid)) ///
        (rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none) lpattern(solid)) ///
        (line  biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
        (line  zero h, lcolor("$mred") lpattern(solid) lwidth(medthick)) ///  <-- red horizontal line at y=0
        if h<=`horizon', ///
        title("{bf:`varname'}", size(10) col(black) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`labname'",  size(7)) ///
        xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white)) ///
        graphregion(color(white)) ///
        legend(off) ///
        name("iv_`var'", replace) ///
        saving("$FIGUREDIR/pirf_ivonly_`var'_Channel.gph", replace)

    graph export "$FIGUREDIR/pirf_ivonly_`var'_Channel.pdf", fontface($grfont) replace

    local ivar = `ivar' + 1
}

foreach stat in se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}


// IVREG: Baseline ***********************************************************************

* controls (one lag each): RGROWTH, RYIELD, REER, and the respective endogenous variable
local cntrls_approval L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER L(1/1).approval 
local cntrls_GOVCRISIS_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER 
local cntrls_STRIKE_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER 
local cntrls_DEMONSTR_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  

* specs for local projection
local p = 0           // p = number of lags of endogenous variable
local ps = 0          // p = number of lags for shock
local horizon = 4     // Impulse horizon
local estdiff = 0      // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32      // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL

* options
global savefigs 1    /* 1: save to disk, 0: don't save */
global verb qui      /* leave empty if you want to display regression results */

* horizons index
cap drop t h
cap gen t = _n
cap gen h = t - 1


local ivar = 1
foreach var in `vars' {

    * set controls
    local cntrls `cntrls_`var''

    * preallocate storage for IRFs
    cap drop biv`var'_base up90biv`var' lo90biv`var' up68biv`var' lo68biv`var'
    qui gen biv`var'_base     = .
    qui gen up90biv`var' = .
    qui gen lo90biv`var' = .
    qui gen up68biv`var' = .
    qui gen lo68biv`var' = .

    * create first-difference if needed (not used when estdiff==2 except for estdiff==1 branch)
    if `estdiff' > 0 {
        cap drop d`var'
        gen d`var' = `var' - L.`var'
    }

    forvalues i = 0/`horizon' {

        local iname `i'

        * cumulative change from t-1 to t+i
        cap drop d`iname'`var'
        gen d`iname'`var' = F`i'.`var' - L.`var'

        if `estdiff' == 0 {
            `verb' ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)
        }
        else if `estdiff' == 1 {
            `verb' ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)
        }
        else if `estdiff' == 2 {
            `verb' ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL) `cntrls' i.year i.country_id, dkraay(1)
        }

        * store coefficient + SE for shock
        cap drop biv`var'h`iname' seiv`var'h`iname'
        gen biv`var'h`iname'_base  = _b[`shock']
        gen seiv`var'h`iname' = _se[`shock']

        qui replace biv`var'_base     = biv`var'h`iname'_base if h==`i'
        qui replace up90biv`var' = biv`var'h`iname'_base + `z1'*seiv`var'h`iname' if h==`i'
        qui replace lo90biv`var' = biv`var'h`iname'_base - `z1'*seiv`var'h`iname' if h==`i'
        qui replace up68biv`var' = biv`var'h`iname'_base + `z2'*seiv`var'h`iname' if h==`i'
        qui replace lo68biv`var' = biv`var'h`iname'_base - `z2'*seiv`var'h`iname' if h==`i'
    }

    * plot helpers
    cap drop zero
    gen zero = 0

    local varname : word `ivar' of `varsnames'
    local labname : word `ivar' of `labels'

    tw ///
        (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none) lpattern(solid)) ///
        (rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none) lpattern(solid)) ///
		(line  biv`var'_base h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
        (line  biv`var' h, lcolor("$mgreen") lpattern(longdash_dot) lwidth(thick)) ///
        (line  zero h, lcolor("$mred") lpattern(solid) lwidth(medthick)) ///  <-- red horizontal line at y=0
        if h<=`horizon', ///
        title("{bf:`varname'}", size(10) col(black) margin(b=2)) ///
        xtitle("Years", size(7)) ///
        ytitle("`labname'",  size(7)) ///
        xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
        plotregion(color(white)) ///
        graphregion(color(white)) ///
        legend(off) ///
        name("iv_`var'", replace) ///
        saving("$FIGUREDIR/pirf_ivonly_`var'_Channel_Base.gph", replace)

    graph export "$FIGUREDIR/pirf_ivonly_`var'_Channel_Base.pdf", fontface($grfont) replace

    local ivar = `ivar' + 1
}

////////////////////////////////////////////////////////////////////////////////
// Combined figure
////////////////////////////////////////////////////////////////////////////////

graph use "$FIGUREDIR/pirf_ivonly_approval_Channel_Base.gph", name(g1, replace)
graph use "$FIGUREDIR/pirf_ivonly_GOVCRISIS_binary_Channel_Base.gph", name(g2, replace)
graph use "$FIGUREDIR/pirf_ivonly_STRIKE_binary_Channel_Base.gph", name(g3, replace)
graph use "$FIGUREDIR/pirf_ivonly_DEMONSTR_binary_Channel_Base.gph", name(g4, replace)
set scheme s1color

grc1leg2 g1 g2 g3 g4, cols(4) ///
    loff ///
	ysize(3) xsize(12) 
graph export "$FIGURECOM/PolEffects_Channel_Base.pdf", replace

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
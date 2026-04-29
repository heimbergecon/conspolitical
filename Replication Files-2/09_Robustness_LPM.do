// Robustness LPM

* To rescale AME later
cap program drop rescale_ame
program rescale_ame, eclass
    matrix b = e(b) * 100
    matrix V = e(V) * 10000
    ereturn repost b=b V=V
end

// Introducy binary variables
drop GOVCRISIS_binary STRIKE_binary DEMONSTR_binary
gen GOVCRISIS_binary = 0
replace GOVCRISIS_binary = 1 if GOVCRISIS > 0
gen STRIKE_binary = 0
replace STRIKE_binary = 1 if STRIKE > 0
gen DEMONSTR_binary = 0
replace DEMONSTR_binary = 1 if DEMONSTR > 0

* Variables of interest
local vars GOVCRISIS_binary STRIKE_binary DEMONSTR_binary 

* Controls (one lag each)
local cntrls_GOVCRISIS_binary L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  
local cntrls_STRIKE_binary    L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  
local cntrls_DEMONSTR_binary  L(1/1).RGROWTH L(1/1).RYIELD L(1/1).REER  

foreach var in `vars' {
    local cntrls `cntrls_`var''
    
    * First stage
    cap drop v_hat
    reg diff_STRUCBAL TOTAL `cntrls' i.year i.country_id
    predict v_hat, residuals
    
    * Second stage: probit + AME
    forvalues h = 0/5 {
        probit F`h'.`var' diff_STRUCBAL `cntrls' i.year i.country_id v_hat, vce(cluster country_id)
        margins, dydx(diff_STRUCBAL) post
		
		rescale_ame  
		
        estimates store `var'_F`h'_probit
    }
    
    esttab `var'_F0_probit `var'_F1_probit `var'_F2_probit ///
       `var'_F3_probit `var'_F4_probit `var'_F5_probit  ///
       using "$TABLEDIR\probit_ame_`var'.csv", replace             ///
    b(3) se(3)                                           ///
    star(+ 0.32 * 0.10)                                  ///
    label                                                ///
    mtitles("h=0" "h=1" "h=2" "h=3" "h=4" "h=5")        ///
    stats(N, fmt(0) labels("Observations"))              ///
    nogaps compress
}

* Panel A: Government Crisis
esttab GOVCRISIS_binary_F0_probit GOVCRISIS_binary_F1_probit GOVCRISIS_binary_F2_probit ///
       GOVCRISIS_binary_F3_probit GOVCRISIS_binary_F4_probit GOVCRISIS_binary_F5_probit ///
       using "$TABLEDIR\probit_ame_combined.csv", replace                                          ///
    keep(diff_STRUCBAL)                                                                  ///
    b(3) se(3)                                                                           ///
    star(+ 0.32 * 0.10)                                                                  ///
    label                                                                                ///
    mtitles("h=0" "h=1" "h=2" "h=3" "h=4" "h=5")                                       ///
    prehead("Panel A: Government Crisis")                                                ///
    stats(N, fmt(0) labels("Observations"))                                              ///
    nogaps compress

* Panel B: Strikes
esttab STRIKE_binary_F0_probit STRIKE_binary_F1_probit STRIKE_binary_F2_probit          ///
       STRIKE_binary_F3_probit STRIKE_binary_F4_probit STRIKE_binary_F5_probit          ///
       using "$TABLEDIR\probit_ame_combined.csv", append                                           ///
    keep(diff_STRUCBAL)                                                                  ///
    b(3) se(3)                                                                           ///
    star(+ 0.32 * 0.10)                                                                  ///
    label                                                                                ///
    nomtitles nonumbers                                                                  ///
    prehead("Panel B: Strikes")                                                          ///
    stats(N, fmt(0) labels("Observations"))                                              ///
    nogaps compress

* Panel C: Demonstrations
esttab DEMONSTR_binary_F0_probit DEMONSTR_binary_F1_probit DEMONSTR_binary_F2_probit    ///
       DEMONSTR_binary_F3_probit DEMONSTR_binary_F4_probit DEMONSTR_binary_F5_probit    ///
       using "$TABLEDIR\probit_ame_combined.csv", append                                           ///
    keep(diff_STRUCBAL)                                                                  ///
    b(3) se(3)                                                                           ///
    star(+ 0.32 * 0.10)                                                                  ///
    label                                                                                ///
    nomtitles nonumbers                                                                  ///
    prehead("Panel C: Demonstrations")                                                   ///
    stats(N, fmt(0) labels("Observations"))                                              ///
    nogaps compress                                                                      ///
    note("Average marginal effects of diff_STRUCBAL reported."                          ///
         "Standard errors clustered by country in parentheses."                         ///
         "+ significant at 68% confidence level (a=0.32)."                              ///
         "* significant at 90% confidence level (a=0.10)."                              ///
         "Control function approach used to address endogeneity."                        ///
         "Comparable to LPM estimates in main results.")
		 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
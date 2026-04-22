/*==============================================================================
  Bitcoin as a Maturing Financial Asset:
  Evidence from Returns, Volatility, and Tail Risk (2016-2026)
  
  MASTER DO-FILE - Full Replication Code
  
  Author:  Mohamed Omar Ata Daoud Nazer
  Program: Bachelor in Economics, IE University
  Supervisor: Gael Sanchez Smith
  Date:    May 2026
  
  Description:
  This file replicates ALL quantitative analyses, tables, and figures reported
  in the thesis. It is organised in three parts:
    Part A - Core Bitcoin Analysis (Tables 1-8, Figures 1-6)
    Part B - Gold Benchmark Analysis
    Part C - Geopolitical Event Study (Tables 9-10, Figures 7-13)
  
  Working directory: C:\Users\mnazer.ieu2022\Documents\
  
  Data sources (all from Investing.com, placed in working directory):
    - Bitcoin_Historical_Data__4_.csv   (Bitcoin USD daily, via CoinMarketCap)
    - XAU_USD_Historical_Data.csv       (Gold XAU/USD daily)
    - S_P_500_Historical_Data__2_.csv   (S&P 500 daily)
  
  Software: Stata 17 or later
  Required user-written packages: estout (installed automatically below)
==============================================================================*/

clear all
set more off

cd "C:\Users\mnazer.ieu2022\Documents"

log close _all
log using "Thesis_RunLog_Final.txt", replace text


*###############################################################################
*  PART A: CORE BITCOIN ANALYSIS
*###############################################################################

*===============================================================================
* A1. IMPORT AND CLEAN BITCOIN DATA
*===============================================================================

import delimited "Bitcoin_Historical_Data__4_.csv", clear varnames(nonames) encoding("utf-8")

rename v1 date_raw
rename v2 price
rename v3 open
rename v4 high
rename v5 low
rename v6 volume_raw
rename v7 change_raw

drop in 1

gen date = date(date_raw, "MDY")
format date %td

destring price open high low, replace ignore(",")

gen change_pct = subinstr(change_raw, "%", "", .)
destring change_pct, replace

gen volume_num = volume_raw
replace volume_num = subinstr(volume_num, "K", "", .)
replace volume_num = subinstr(volume_num, "M", "", .)
replace volume_num = subinstr(volume_num, "B", "", .)
destring volume_num, replace
replace volume_num = volume_num * 1000       if strpos(volume_raw, "K") > 0
replace volume_num = volume_num * 1000000    if strpos(volume_raw, "M") > 0
replace volume_num = volume_num * 1000000000 if strpos(volume_raw, "B") > 0

sort date
tsset date, daily

*===============================================================================
* A2. GENERATE CORE VARIABLES
*===============================================================================

gen ln_price = ln(price)
gen ret      = ln_price - L.ln_price
gen ret_sq   = ret^2

gen regime2 = inrange(date, td(01jan2021), td(31dec2023))
gen regime3 = inrange(date, td(01jan2024), td(01jan2026))

gen etf_event = (date == td(10jan2024))
gen etf_win_1 = (date == td(10jan2024))
gen etf_win_3 = inrange(date, td(07jan2024), td(13jan2024))
gen etf_win_7 = inrange(date, td(03jan2024), td(17jan2024))

gen sample_main = !missing(ret, regime2, regime3, etf_win_3) & date <= td(01jan2026)

ssc install estout, replace

*===============================================================================
* A3. TABLE 1 - Full-Sample Descriptive Statistics
*===============================================================================

estpost summarize ret ret_sq price volume_num if sample_main==1
esttab using "Table1_Descriptives.rtf", replace ///
    cells("count mean sd min max") title("Descriptive Statistics")

*===============================================================================
* A4. TABLE 2 - Descriptive Statistics by Regime
*===============================================================================

estpost summarize ret ret_sq if sample_main==1 & regime2==0 & regime3==0
est store reg1
estpost summarize ret ret_sq if sample_main==1 & regime2==1
est store reg2
estpost summarize ret ret_sq if sample_main==1 & regime3==1
est store reg3
esttab reg1 reg2 reg3 using "Table2_ByRegime.rtf", replace ///
    cells("mean sd") mtitles("2016-2020" "2021-2023" "2024-2026")

*===============================================================================
* A5. TABLE 3 - Baseline Return Regressions
*===============================================================================

eststo clear
reg ret regime2 regime3 if sample_main==1, vce(robust)
eststo m1
reg ret regime2 regime3 etf_win_3 if sample_main==1, vce(robust)
eststo m2
esttab m1 m2 using "Table3_Returns.rtf", replace ///
    se star(* 0.10 ** 0.05 *** 0.01)

*===============================================================================
* A6. TABLE 4 - Volatility Regression (Squared Returns)
*===============================================================================

eststo clear
reg ret_sq regime2 regime3 if sample_main==1, vce(robust)
eststo v1
esttab v1 using "Table4_Volatility.rtf", replace ///
    se star(* 0.10 ** 0.05 *** 0.01)

*===============================================================================
* A7. TABLE 5 - ETF Event-Window Regression
*===============================================================================

eststo clear
reg ret etf_win_3 if sample_main==1, vce(robust)
eststo e1
reg ret regime2 regime3 etf_win_3 if sample_main==1, vce(robust)
eststo e2
esttab e1 e2 using "Table5_Event.rtf", replace ///
    se star(* 0.10 ** 0.05 *** 0.01)

*===============================================================================
* A8. TABLE 8 - ETF Event-Window Robustness (1-day, 3-day, 7-day)
*===============================================================================

eststo clear
reg ret etf_win_1 if sample_main==1, vce(robust)
eststo w1
reg ret etf_win_3 if sample_main==1, vce(robust)
eststo w3
reg ret etf_win_7 if sample_main==1, vce(robust)
eststo w7
esttab w1 w3 w7 using "Table8_ETF_Windows.rtf", replace ///
    se star(* 0.10 ** 0.05 *** 0.01)

*===============================================================================
* A9. ARCH-LM TEST AND GARCH(1,1) ESTIMATION
*===============================================================================

reg ret if sample_main==1
predict resid, resid
estat archlm

arch ret if sample_main==1, arch(1) garch(1) distribution(t)

predict h_t, variance
gen vol_garch = sqrt(h_t)

*===============================================================================
* A10. FIGURES 1-4 - Core Time-Series Plots (with corrected axis labels)
*
*  FIXES APPLIED:
*    - Added xtitle("Date") and descriptive ytitle() to all graphs
*    - Added explicit xlabel year labels to prevent x-axis date truncation
*    - Added graphregion(margin(r=4)) for right-margin breathing room
*    - Added width(1200) to all graph exports for resolution
*===============================================================================

* Figure 1: Daily Bitcoin Returns
tsline ret if sample_main==1, ///
    ytitle("Daily Log Return") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure1_Returns.png", replace width(1200)

* Figure 2: Squared Returns (Volatility Proxy)
tsline ret_sq if sample_main==1, ///
    ytitle("Squared Returns") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure2_Volatility.png", replace width(1200)

* Figure 2 (Zoomed): Squared Returns 2021-2026
tsline ret_sq if sample_main==1 & date >= td(01jan2021), ///
    ytitle("Squared Returns") xtitle("Date") ///
    xlabel(`=td(01jan2021)' "2021" `=td(01jan2022)' "2022" `=td(01jan2023)' "2023" `=td(01jan2024)' "2024" `=td(01jan2025)' "2025" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure2_Volatility_Zoom.png", replace width(1200)

* Figure 3: GARCH(1,1) Conditional Volatility
tsline vol_garch if sample_main==1, ///
    ytitle("GARCH Conditional Volatility") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure3_GARCH.png", replace width(1200)

* Figure 4: Bitcoin Price
tsline price if sample_main==1, ///
    ytitle("Bitcoin Price (USD)") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure4_Price.png", replace width(1200)

*===============================================================================
* A11. TABLE 6 - Tail-Risk Frequencies by Regime
*===============================================================================

gen abs_ret = abs(ret)
gen tail_3  = (abs_ret > 0.03) if sample_main==1 & !missing(ret)
gen tail_5  = (abs_ret > 0.05) if sample_main==1 & !missing(ret)
gen tail_10 = (abs_ret > 0.10) if sample_main==1 & !missing(ret)

preserve
    collapse (mean) tail_3 tail_5 tail_10 if sample_main==1 & regime2==0 & regime3==0
    gen regime = "2016-2020"
    tempfile t1
    save `t1'
restore
preserve
    collapse (mean) tail_3 tail_5 tail_10 if sample_main==1 & regime2==1
    gen regime = "2021-2023"
    tempfile t2
    save `t2'
restore
preserve
    collapse (mean) tail_3 tail_5 tail_10 if sample_main==1 & regime3==1
    gen regime = "2024-2026"
    tempfile t3
    save `t3'
restore

preserve
    use `t1', clear
    append using `t2'
    append using `t3'
    replace tail_3  = tail_3  * 100
    replace tail_5  = tail_5  * 100
    replace tail_10 = tail_10 * 100
    order regime
    list, noobs
    export delimited "Table6_TailRisk.csv", replace
restore

*===============================================================================
* A12. TABLES 7 & 7b - Drawdown Analysis by Regime
*===============================================================================

gen cummax = price if _n == 1
replace cummax = max(price, cummax[_n-1]) if _n > 1 & !missing(price)
gen drawdown = (price - cummax) / cummax

* Table 7: Max and mean drawdown by regime
preserve
    collapse (min) max_drawdown=drawdown (mean) mean_drawdown=drawdown ///
        if sample_main==1 & regime2==0 & regime3==0
    gen regime = "2016-2020"
    tempfile d1
    save `d1'
restore
preserve
    collapse (min) max_drawdown=drawdown (mean) mean_drawdown=drawdown ///
        if sample_main==1 & regime2==1
    gen regime = "2021-2023"
    tempfile d2
    save `d2'
restore
preserve
    collapse (min) max_drawdown=drawdown (mean) mean_drawdown=drawdown ///
        if sample_main==1 & regime3==1
    gen regime = "2024-2026"
    tempfile d3
    save `d3'
restore

preserve
    use `d1', clear
    append using `d2'
    append using `d3'
    order regime
    list, noobs
    export delimited "Table7_Drawdown.csv", replace
restore

* Table 7b: Percentage of days spent in drawdown
gen in_drawdown = (drawdown < 0) if sample_main==1 & !missing(drawdown)

preserve
    collapse (mean) pct_drawdown_days=in_drawdown ///
        if sample_main==1 & regime2==0 & regime3==0
    gen regime = "2016-2020"
    tempfile dd1
    save `dd1'
restore
preserve
    collapse (mean) pct_drawdown_days=in_drawdown ///
        if sample_main==1 & regime2==1
    gen regime = "2021-2023"
    tempfile dd2
    save `dd2'
restore
preserve
    collapse (mean) pct_drawdown_days=in_drawdown ///
        if sample_main==1 & regime3==1
    gen regime = "2024-2026"
    tempfile dd3
    save `dd3'
restore

preserve
    use `dd1', clear
    append using `dd2'
    append using `dd3'
    replace pct_drawdown_days = pct_drawdown_days * 100
    order regime
    list, noobs
    export delimited "Table7b_DrawdownDays.csv", replace
restore

*===============================================================================
* A13. FIGURES 5-6 - Rolling Volatility and Rolling Tail-Risk Frequency
*===============================================================================

* Figure 5: 30-Day Rolling Standard Deviation
gen roll_sd_30 = .
forvalues i = 30/`=_N' {
    quietly summ ret if _n >= `i' - 29 & _n <= `i' & sample_main==1
    if r(N) >= 20 {
        quietly replace roll_sd_30 = r(sd) in `i'
    }
}

twoway (line roll_sd_30 date if sample_main==1, lcolor(navy) lwidth(medium)), ///
    title("30-Day Rolling Volatility of Bitcoin Returns", size(medium)) ///
    ytitle("Rolling Standard Deviation") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure5_RollingVolatility.png", replace width(1200)

* Figure 6: 90-Day Rolling Frequency of |ret| > 5%
gen roll_tail_90 = .
forvalues i = 90/`=_N' {
    quietly summ tail_5 if _n >= `i' - 89 & _n <= `i' & sample_main==1
    if r(N) >= 60 {
        quietly replace roll_tail_90 = r(mean) in `i'
    }
}

twoway (line roll_tail_90 date if sample_main==1, lcolor(navy) lwidth(medium)), ///
    title("90-Day Rolling Frequency of Absolute Returns Above 5%", size(medium)) ///
    ytitle("Proportion of Days") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4))
graph export "Figure6_RollingTailRisk.png", replace width(1200)

*===============================================================================
* A14. FIGURE 7 - ETF Event-Window Price Chart
*===============================================================================

twoway (line price date if inrange(date, td(15nov2023), td(25feb2024)), ///
        lcolor(navy) lwidth(medium)) ///
    , xline(`=td(10jan2024)', lcolor(cranberry) lpattern(solid)) ///
    title("Bitcoin Price Around Spot ETF Approval", size(medium)) ///
    ytitle("Price (USD)") xtitle("Date") ///
    xlabel(`=td(26nov2023)' "26nov2023" `=td(17dec2023)' "17dec2023" `=td(08jan2024)' "08jan2024" `=td(29jan2024)' "29jan2024" `=td(19feb2024)' "19feb2024", labsize(small) angle(45)) ///
    graphregion(margin(r=4))
graph export "Figure7_ETF_Event.png", replace width(1200)

save "Thesis_Final_Clean.dta", replace


*###############################################################################
*  PART B: GOLD BENCHMARK ANALYSIS
*###############################################################################

*===============================================================================
* B1. IMPORT AND CLEAN GOLD DATA (XAU/USD)
*===============================================================================

import delimited "XAU_USD_Historical_Data.csv", clear varnames(nonames) encoding("utf-8")
drop in 1

gen date = date(v1, "MDY")
format date %td
destring v2, gen(gold_price) ignore(",")

sort date
tsset date, daily

gen gold_ln  = ln(gold_price)
gen gold_ret = gold_ln - L.gold_ln
gen gold_ret_sq = gold_ret^2

gen regime2 = inrange(date, td(01jan2021), td(31dec2023))
gen regime3 = inrange(date, td(01jan2024), td(01jan2026))

*===============================================================================
* B2. GOLD DESCRIPTIVE STATISTICS BY REGIME
*===============================================================================

preserve
    collapse (mean) mean_ret=gold_ret mean_sq=gold_ret_sq ///
             (sd) sd_ret=gold_ret sd_sq=gold_ret_sq ///
             (min) min_ret=gold_ret min_sq=gold_ret_sq ///
             (max) max_ret=gold_ret max_sq=gold_ret_sq ///
        if regime2==0 & regime3==0 & !missing(gold_ret)
    gen regime = "2016-2020"
    tempfile gd1
    save `gd1'
restore
preserve
    collapse (mean) mean_ret=gold_ret mean_sq=gold_ret_sq ///
             (sd) sd_ret=gold_ret sd_sq=gold_ret_sq ///
             (min) min_ret=gold_ret min_sq=gold_ret_sq ///
             (max) max_ret=gold_ret max_sq=gold_ret_sq ///
        if regime2==1 & !missing(gold_ret)
    gen regime = "2021-2023"
    tempfile gd2
    save `gd2'
restore
preserve
    collapse (mean) mean_ret=gold_ret mean_sq=gold_ret_sq ///
             (sd) sd_ret=gold_ret sd_sq=gold_ret_sq ///
             (min) min_ret=gold_ret min_sq=gold_ret_sq ///
             (max) max_ret=gold_ret max_sq=gold_ret_sq ///
        if regime3==1 & !missing(gold_ret)
    gen regime = "2024-2026"
    tempfile gd3
    save `gd3'
restore
preserve
    collapse (mean) mean_ret=gold_ret mean_sq=gold_ret_sq ///
             (sd) sd_ret=gold_ret sd_sq=gold_ret_sq ///
             (min) min_ret=gold_ret min_sq=gold_ret_sq ///
             (max) max_ret=gold_ret max_sq=gold_ret_sq ///
        if !missing(gold_ret)
    gen regime = "Full Sample"
    tempfile gd_all
    save `gd_all'
restore

preserve
    use `gd1', clear
    append using `gd2'
    append using `gd3'
    append using `gd_all'
    order regime
    list, noobs
    export delimited "Gold_Descriptives_ByRegime.csv", replace
restore

*===============================================================================
* B3. GOLD TAIL-RISK FREQUENCIES
*===============================================================================

gen gold_abs = abs(gold_ret)
gen gold_tail_3  = (gold_abs > 0.03) if !missing(gold_ret)
gen gold_tail_5  = (gold_abs > 0.05) if !missing(gold_ret)
gen gold_tail_10 = (gold_abs > 0.10) if !missing(gold_ret)

preserve
    collapse (mean) gold_tail_3 gold_tail_5 gold_tail_10 ///
        if regime2==0 & regime3==0 & !missing(gold_ret)
    gen regime = "2016-2020"
    tempfile gt1
    save `gt1'
restore
preserve
    collapse (mean) gold_tail_3 gold_tail_5 gold_tail_10 ///
        if regime2==1 & !missing(gold_ret)
    gen regime = "2021-2023"
    tempfile gt2
    save `gt2'
restore
preserve
    collapse (mean) gold_tail_3 gold_tail_5 gold_tail_10 ///
        if regime3==1 & !missing(gold_ret)
    gen regime = "2024-2026"
    tempfile gt3
    save `gt3'
restore
preserve
    collapse (mean) gold_tail_3 gold_tail_5 gold_tail_10 ///
        if !missing(gold_ret)
    gen regime = "Full Sample"
    tempfile gt_all
    save `gt_all'
restore

preserve
    use `gt1', clear
    append using `gt2'
    append using `gt3'
    append using `gt_all'
    replace gold_tail_3  = gold_tail_3  * 100
    replace gold_tail_5  = gold_tail_5  * 100
    replace gold_tail_10 = gold_tail_10 * 100
    order regime
    list, noobs
    export delimited "Gold_TailRisk.csv", replace
restore

*===============================================================================
* B4. GOLD DRAWDOWN ANALYSIS
*===============================================================================

gen gold_cummax = gold_price if _n == 1
replace gold_cummax = max(gold_price, gold_cummax[_n-1]) if _n > 1 & !missing(gold_price)
gen gold_drawdown = (gold_price - gold_cummax) / gold_cummax

preserve
    collapse (min) max_drawdown=gold_drawdown (mean) mean_drawdown=gold_drawdown ///
        if regime2==0 & regime3==0 & !missing(gold_drawdown)
    gen regime = "2016-2020"
    tempfile gdd1
    save `gdd1'
restore
preserve
    collapse (min) max_drawdown=gold_drawdown (mean) mean_drawdown=gold_drawdown ///
        if regime2==1 & !missing(gold_drawdown)
    gen regime = "2021-2023"
    tempfile gdd2
    save `gdd2'
restore
preserve
    collapse (min) max_drawdown=gold_drawdown (mean) mean_drawdown=gold_drawdown ///
        if regime3==1 & !missing(gold_drawdown)
    gen regime = "2024-2026"
    tempfile gdd3
    save `gdd3'
restore
preserve
    collapse (min) max_drawdown=gold_drawdown (mean) mean_drawdown=gold_drawdown ///
        if !missing(gold_drawdown)
    gen regime = "Full Sample"
    tempfile gdd_all
    save `gdd_all'
restore

preserve
    use `gdd1', clear
    append using `gdd2'
    append using `gdd3'
    append using `gdd_all'
    order regime
    list, noobs
    export delimited "Gold_Drawdown.csv", replace
restore

*===============================================================================
* B5. GOLD FIGURES - Price and Rolling Volatility
*===============================================================================

* Gold Price (USD)
twoway (line gold_price date, lcolor(gold) lwidth(medium)), ///
    title("Gold Price Over Time (USD per Ounce)", size(medium)) ///
    ytitle("Gold Price (USD/oz)") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4)) ///
    note("Source: Investing.com (XAU/USD)", size(vsmall))
graph export "Gold_Price_USD.png", replace width(1200)

* Gold 30-Day Rolling Volatility
gen gold_roll_sd = .
gen t_gold = _n

summ t_gold, meanonly
local tmax = r(max)

forvalues j = 30/`tmax' {
    local start = `j' - 29
    quietly summ gold_ret if t_gold >= `start' & t_gold <= `j' & !missing(gold_ret)
    if r(N) >= 20 {
        replace gold_roll_sd = r(sd) in `j'
    }
}

twoway (line gold_roll_sd date if !missing(gold_roll_sd), lcolor(gold) lwidth(medium)), ///
    title("30-Day Rolling Volatility of Gold Returns", size(medium)) ///
    ytitle("Rolling Standard Deviation") xtitle("Date") ///
    xlabel(`=td(01jan2016)' "2016" `=td(01jan2018)' "2018" `=td(01jan2020)' "2020" `=td(01jan2022)' "2022" `=td(01jan2024)' "2024" `=td(01jan2026)' "2026", labsize(small)) ///
    graphregion(margin(r=4)) ///
    note("Source: Author's calculations from Investing.com (XAU/USD)", size(vsmall))
graph export "Gold_RollingVolatility_USD.png", replace width(1200)


*###############################################################################
*  PART C: GEOPOLITICAL EVENT STUDY (THREE-ASSET COMPARISON)
*###############################################################################

*===============================================================================
* C1. IMPORT AND MERGE ALL THREE ASSETS
*===============================================================================

* --- Bitcoin ---
import delimited "Bitcoin_Historical_Data__4_.csv", clear varnames(nonames) encoding("utf-8")
drop in 1
gen date = date(v1, "MDY")
format date %td
destring v2, gen(btc_price) ignore(",")
sort date
tsset date, daily
gen ln_price = ln(btc_price)
gen btc_ret = ln_price - L.ln_price
keep date btc_price btc_ret
tempfile btc
save `btc'

* --- Gold ---
import delimited "XAU_USD_Historical_Data.csv", clear varnames(nonames) encoding("utf-8")
drop in 1
gen date = date(v1, "MDY")
format date %td
destring v2, gen(gold_price) ignore(",")
sort date
gen gold_ln = ln(gold_price)
gen gold_ret = gold_ln - gold_ln[_n-1] if _n > 1
keep date gold_price gold_ret
tempfile gold
save `gold'

* --- S&P 500 ---
import delimited "S_P_500_Historical_Data__2_.csv", clear varnames(nonames) encoding("utf-8")
drop in 1
gen date = date(v1, "MDY")
format date %td
destring v2, gen(spx_price) ignore(",")
sort date
gen spx_ln = ln(spx_price)
gen spx_ret = spx_ln - spx_ln[_n-1] if _n > 1
keep date spx_price spx_ret
tempfile spx
save `spx'

* --- Merge ---
use `btc', clear
merge 1:1 date using `gold', nogen
merge 1:1 date using `spx', nogen
sort date
gen t = _n
tsset t

save "Event_Merged.dta", replace

summ date, meanonly
local last_date_val = r(max)
local last_date_fmt : di %td r(max)

*===============================================================================
* C2. DEFINE GEOPOLITICAL EVENTS
*===============================================================================

local ev1_name "US-Iran Escalation"
local ev1_date = td(03jan2020)
local ev1_short "US-Iran"

local ev2_name "COVID Pandemic Declaration"
local ev2_date = td(11mar2020)
local ev2_short "COVID"

local ev3_name "Russia Invasion of Ukraine"
local ev3_date = td(24feb2022)
local ev3_short "Ukraine"

local ev4_name "US Regional Banking Crisis"
local ev4_date = td(09mar2023)
local ev4_short "SVB Crisis"

local ev5_name "Yen Carry Trade Unwinding"
local ev5_date = td(05aug2024)
local ev5_short "Yen Unwind"

local ev6_name "Trump Liberation Day"
local ev6_date = td(02apr2025)
local ev6_short "Liberation Day"

local ev7_name "Iran Conflict"
local ev7_date = td(28feb2026)
local ev7_short "Iran 2026"

local n_events = 7

*===============================================================================
* C3. COMPUTE ALL EVENT RETURNS
*===============================================================================

clear
set obs `n_events'

gen event_id = _n
gen str40 event_name = ""
gen str20 event_short = ""
gen str12 event_date_str = ""
gen event_date = .
gen btc_1d = .
gen gold_1d = .
gen spx_1d = .
gen btc_5d = .
gen gold_5d = .
gen spx_5d = .
gen btc_30d = .
gen gold_30d = .
gen spx_30d = .
gen btc_60d = .
gen gold_60d = .
gen spx_60d = .
gen str30 note_60d = ""

tempfile results_shell
save `results_shell'

use "Event_Merged.dta", clear

forvalues i = 1/`n_events' {
    
    local edate = `ev`i'_date'
    local edate_plus5 = `edate' + 5
    local edate_plus30 = `edate' + 30
    local edate_plus60 = `edate' + 60
    
    local is_incomplete = (`edate_plus60' > `last_date_val')
    
    * ===== 1-DAY RETURNS =====
    summ t if date >= `edate' & !missing(btc_ret), meanonly
    if r(N) > 0 {
        local btc_1d_`i' = btc_ret[r(min)] * 100
    }
    else {
        local btc_1d_`i' = .
    }
    
    summ t if date >= `edate' & !missing(gold_ret), meanonly
    if r(N) > 0 {
        local gold_1d_`i' = gold_ret[r(min)] * 100
    }
    else {
        local gold_1d_`i' = .
    }
    
    summ t if date >= `edate' & !missing(spx_ret), meanonly
    if r(N) > 0 {
        local spx_1d_`i' = spx_ret[r(min)] * 100
    }
    else {
        local spx_1d_`i' = .
    }
    
    * ===== 5-DAY RETURNS =====
    summ t if date >= `edate' & !missing(btc_price), meanonly
    local t0b = r(min)
    summ t if date >= `edate_plus5' & !missing(btc_price), meanonly
    local t1b = r(min)
    if !missing(`t0b') & !missing(`t1b') {
        local btc_5d_`i' = (btc_price[`t1b'] - btc_price[`t0b']) / btc_price[`t0b'] * 100
    }
    else {
        local btc_5d_`i' = .
    }
    
    summ t if date >= `edate' & !missing(gold_price), meanonly
    local t0g5 = r(min)
    summ t if date >= `edate_plus5' & !missing(gold_price), meanonly
    local t1g5 = r(min)
    if !missing(`t0g5') & !missing(`t1g5') {
        local gold_5d_`i' = (gold_price[`t1g5'] - gold_price[`t0g5']) / gold_price[`t0g5'] * 100
    }
    else {
        local gold_5d_`i' = .
    }
    
    summ t if date >= `edate' & !missing(spx_price), meanonly
    local t0s5 = r(min)
    summ t if date >= `edate_plus5' & !missing(spx_price), meanonly
    local t1s5 = r(min)
    if !missing(`t0s5') & !missing(`t1s5') {
        local spx_5d_`i' = (spx_price[`t1s5'] - spx_price[`t0s5']) / spx_price[`t0s5'] * 100
    }
    else {
        local spx_5d_`i' = .
    }
    
    * ===== 60-DAY RETURNS =====
    local end60 = min(`edate_plus60', `last_date_val')
    
    summ t if date >= `edate' & !missing(btc_price), meanonly
    local t0b60 = r(min)
    summ t if date >= `end60' & !missing(btc_price), meanonly
    local t1b60 = r(min)
    if !missing(`t0b60') & !missing(`t1b60') & `t1b60' > `t0b60' {
        local btc_60d_`i' = (btc_price[`t1b60'] - btc_price[`t0b60']) / btc_price[`t0b60'] * 100
    }
    else {
        local btc_60d_`i' = .
    }
    
    summ t if date >= `edate' & !missing(gold_price), meanonly
    local t0g60 = r(min)
    summ t if date >= `end60' & !missing(gold_price), meanonly
    local t1g60 = r(min)
    if !missing(`t0g60') & !missing(`t1g60') & `t1g60' > `t0g60' {
        local gold_60d_`i' = (gold_price[`t1g60'] - gold_price[`t0g60']) / gold_price[`t0g60'] * 100
    }
    else {
        local gold_60d_`i' = .
    }
    
    summ t if date >= `edate' & !missing(spx_price), meanonly
    local t0s60 = r(min)
    summ t if date >= `end60' & !missing(spx_price), meanonly
    local t1s60 = r(min)
    if !missing(`t0s60') & !missing(`t1s60') & `t1s60' > `t0s60' {
        local spx_60d_`i' = (spx_price[`t1s60'] - spx_price[`t0s60']) / spx_price[`t0s60'] * 100
    }
    else {
        local spx_60d_`i' = .
    }
    
    if `is_incomplete' {
        local note_`i' "As of `last_date_fmt'"
    }
    else {
        local note_`i' ""
    }
    
    local date_fmt : di %td `edate'
    di as text "`ev`i'_name' (`date_fmt'): BTC=" %6.1f `btc_60d_`i'' ///
       "  Gold=" %6.1f `gold_60d_`i'' "  SPX=" %6.1f `spx_60d_`i'' "  `note_`i''"
}

* Fill results dataset
use `results_shell', clear

forvalues i = 1/`n_events' {
    replace event_name = "`ev`i'_name'" in `i'
    replace event_short = "`ev`i'_short'" in `i'
    local edate = `ev`i'_date'
    local date_fmt : di %td `edate'
    replace event_date_str = "`date_fmt'" in `i'
    replace event_date = `edate' in `i'
    replace note_60d = "`note_`i''" in `i'
    
    replace btc_1d = `btc_1d_`i'' in `i'
    replace gold_1d = `gold_1d_`i'' in `i'
    replace spx_1d = `spx_1d_`i'' in `i'
    replace btc_5d = `btc_5d_`i'' in `i'
    replace gold_5d = `gold_5d_`i'' in `i'
    replace spx_5d = `spx_5d_`i'' in `i'
    replace btc_60d = `btc_60d_`i'' in `i'
    replace gold_60d = `gold_60d_`i'' in `i'
    replace spx_60d = `spx_60d_`i'' in `i'
}

format event_date %td
save "Event_Returns_Results.dta", replace

*===============================================================================
* C4. EXPORT TABLES 9 & 10
*===============================================================================

* --- RTF Table 9: 60-Day Returns ---
tempname fh
file open `fh' using "Table9_Event60d.rtf", write replace
file write `fh' "{\rtf1\ansi" _n
file write `fh' "{\b Table 9: S&P 500, Gold, and Bitcoin Through Geopolitical Events \emdash  60-Day Returns (%)}\line\line" _n
file write `fh' "\trowd\trgaph100" _n
file write `fh' "\cellx3200\cellx4400\cellx5500\cellx6600\cellx7700\cellx9000" _n
file write `fh' "{\b Event}\cell {\b Date}\cell {\b S&P 500}\cell {\b Gold}\cell {\b Bitcoin}\cell {\b Note}\cell\row" _n

forvalues i = 1/`n_events' {
    local ename = event_name[`i']
    local dstr = event_date_str[`i']
    local s60 : di %5.1f spx_60d[`i']
    local g60 : di %5.1f gold_60d[`i']
    local b60 : di %5.1f btc_60d[`i']
    local nt = note_60d[`i']
    file write `fh' "`ename'\cell `dstr'\cell `s60'\cell `g60'\cell `b60'\cell `nt'\cell\row" _n
}

file write `fh' "\line" _n
file write `fh' "{\i Source: Author's calculations from CoinMarketCap (Bitcoin) and Investing.com (S&P 500, XAU/USD).}\line" _n
file write `fh' "{\i Returns are simple percentage changes from closing price on event date to closing price 60 calendar days later.}\line" _n
file write `fh' "{\i Where the 60-day window extends beyond available data, returns are computed to the last available date.}" _n
file write `fh' "}" _n
file close `fh'

* --- RTF Table 10: Immediate Impact ---
tempname fh2
file open `fh2' using "Table10_EventImmediate.rtf", write replace
file write `fh2' "{\rtf1\ansi" _n
file write `fh2' "{\b Table 10: Immediate Impact of Geopolitical Events on Asset Returns (%)}\line\line" _n
file write `fh2' "\trowd\trgaph80" _n
file write `fh2' "\cellx2600\cellx3600\cellx4600\cellx5600\cellx6600\cellx7600\cellx8600" _n
file write `fh2' "{\b Event}\cell {\b SPX 1d}\cell {\b Gold 1d}\cell {\b BTC 1d}\cell {\b SPX 5d}\cell {\b Gold 5d}\cell {\b BTC 5d}\cell\row" _n

forvalues i = 1/`n_events' {
    local ename = event_short[`i']
    local s1 : di %5.2f spx_1d[`i']
    local g1 : di %5.2f gold_1d[`i']
    local b1 : di %5.2f btc_1d[`i']
    local s5 : di %5.2f spx_5d[`i']
    local g5 : di %5.2f gold_5d[`i']
    local b5 : di %5.2f btc_5d[`i']
    file write `fh2' "`ename'\cell `s1'\cell `g1'\cell `b1'\cell `s5'\cell `g5'\cell `b5'\cell\row" _n
}

file write `fh2' "\line" _n
file write `fh2' "{\i Source: Author's calculations. 1-day returns are log returns on event date (or nearest trading day).}\line" _n
file write `fh2' "{\i 5-day returns are simple returns from event date to event date + 5 calendar days.}" _n
file write `fh2' "}" _n
file close `fh2'

*===============================================================================
* C5. FIGURES 8-9 - Event Study Bar Charts
*===============================================================================

* Figure 8: 60-Day Cumulative Returns
graph bar spx_60d gold_60d btc_60d, ///
    over(event_short, sort(event_id) label(angle(45) labsize(small))) ///
    legend(order(1 "S&P 500" 2 "Gold" 3 "Bitcoin") position(6) rows(1)) ///
    title("60-Day Cumulative Returns After Geopolitical Events", size(medium)) ///
    subtitle("S&P 500, Gold, and Bitcoin (%)", size(small)) ///
    ytitle("Cumulative Return (%)") ///
    bar(1, color(navy)) ///
    bar(2, color(gold)) ///
    bar(3, color(orange)) ///
    yline(0, lcolor(black) lwidth(thin)) ///
    note("Source: Author's calculations from CoinMarketCap and Investing.com." ///
         "Iran 2026 return computed to last available date.", size(vsmall))
graph export "Figure8_Event60d.png", replace width(1200)

* Figure 9: 1-Day Returns
graph bar spx_1d gold_1d btc_1d, ///
    over(event_short, sort(event_id) label(angle(45) labsize(small))) ///
    legend(order(1 "S&P 500" 2 "Gold" 3 "Bitcoin") position(6) rows(1)) ///
    title("1-Day Returns on Geopolitical Event Dates", size(medium)) ///
    subtitle("S&P 500, Gold, and Bitcoin (%)", size(small)) ///
    ytitle("1-Day Return (%)") ///
    bar(1, color(navy)) ///
    bar(2, color(gold)) ///
    bar(3, color(orange)) ///
    yline(0, lcolor(black) lwidth(thin)) ///
    note("Source: Author's calculations from CoinMarketCap and Investing.com", size(vsmall))
graph export "Figure9_EventImmediate.png", replace width(1200)

*===============================================================================
* C6. FIGURES 10-13 - Indexed Price Path Graphs
*===============================================================================

use "Event_Merged.dta", clear

* --- Figure 10: COVID (March 11, 2020) ---
local ev_date = td(11mar2020)
local win_start = `ev_date' - 30
local win_end = `ev_date' + 60

summ t if date >= `ev_date' & !missing(btc_price), meanonly
local btc_base = btc_price[r(min)]
summ t if date >= `ev_date' & !missing(gold_price), meanonly
local gold_base = gold_price[r(min)]
summ t if date >= `ev_date' & !missing(spx_price), meanonly
local spx_base = spx_price[r(min)]

gen btc_idx = (btc_price / `btc_base') * 100 if date >= `win_start' & date <= `win_end'
gen gold_idx = (gold_price / `gold_base') * 100 if date >= `win_start' & date <= `win_end'
gen spx_idx = (spx_price / `spx_base') * 100 if date >= `win_start' & date <= `win_end'

twoway (line btc_idx date, lcolor(orange) lwidth(medthick)) ///
       (line gold_idx date, lcolor(gold) lwidth(medthick)) ///
       (line spx_idx date, lcolor(navy) lwidth(medthick)) ///
       if date >= `win_start' & date <= `win_end', ///
    xline(`ev_date', lcolor(cranberry) lpattern(dash)) ///
    title("Asset Prices Around COVID Pandemic Declaration", size(medium)) ///
    subtitle("Indexed to 100 on March 11, 2020", size(small)) ///
    ytitle("Indexed Price (Event Date = 100)") xtitle("Date") ///
    legend(order(1 "Bitcoin" 2 "Gold" 3 "S&P 500") position(6) rows(1)) ///
    note("Vertical line: WHO pandemic declaration (Mar 11, 2020)", size(vsmall))
graph export "Figure10_COVID_Window.png", replace width(1200)
drop btc_idx gold_idx spx_idx

* --- Figure 11: Russia-Ukraine (February 24, 2022) ---
local ev_date = td(24feb2022)
local win_start = `ev_date' - 30
local win_end = `ev_date' + 60

summ t if date >= `ev_date' & !missing(btc_price), meanonly
local btc_base = btc_price[r(min)]
summ t if date >= `ev_date' & !missing(gold_price), meanonly
local gold_base = gold_price[r(min)]
summ t if date >= `ev_date' & !missing(spx_price), meanonly
local spx_base = spx_price[r(min)]

gen btc_idx = (btc_price / `btc_base') * 100 if date >= `win_start' & date <= `win_end'
gen gold_idx = (gold_price / `gold_base') * 100 if date >= `win_start' & date <= `win_end'
gen spx_idx = (spx_price / `spx_base') * 100 if date >= `win_start' & date <= `win_end'

twoway (line btc_idx date, lcolor(orange) lwidth(medthick)) ///
       (line gold_idx date, lcolor(gold) lwidth(medthick)) ///
       (line spx_idx date, lcolor(navy) lwidth(medthick)) ///
       if date >= `win_start' & date <= `win_end', ///
    xline(`ev_date', lcolor(cranberry) lpattern(dash)) ///
    title("Asset Prices Around Russia's Invasion of Ukraine", size(medium)) ///
    subtitle("Indexed to 100 on February 24, 2022", size(small)) ///
    ytitle("Indexed Price (Event Date = 100)") xtitle("Date") ///
    legend(order(1 "Bitcoin" 2 "Gold" 3 "S&P 500") position(6) rows(1)) ///
    note("Vertical line: start of invasion (Feb 24, 2022)", size(vsmall))
graph export "Figure11_Ukraine_Window.png", replace width(1200)
drop btc_idx gold_idx spx_idx

* --- Figure 12: SVB / Banking Crisis (March 9, 2023) ---
local ev_date = td(09mar2023)
local win_start = `ev_date' - 30
local win_end = `ev_date' + 60

summ t if date >= `ev_date' & !missing(btc_price), meanonly
local btc_base = btc_price[r(min)]
summ t if date >= `ev_date' & !missing(gold_price), meanonly
local gold_base = gold_price[r(min)]
summ t if date >= `ev_date' & !missing(spx_price), meanonly
local spx_base = spx_price[r(min)]

gen btc_idx = (btc_price / `btc_base') * 100 if date >= `win_start' & date <= `win_end'
gen gold_idx = (gold_price / `gold_base') * 100 if date >= `win_start' & date <= `win_end'
gen spx_idx = (spx_price / `spx_base') * 100 if date >= `win_start' & date <= `win_end'

twoway (line btc_idx date, lcolor(orange) lwidth(medthick)) ///
       (line gold_idx date, lcolor(gold) lwidth(medthick)) ///
       (line spx_idx date, lcolor(navy) lwidth(medthick)) ///
       if date >= `win_start' & date <= `win_end', ///
    xline(`ev_date', lcolor(cranberry) lpattern(dash)) ///
    title("Asset Prices Around US Regional Banking Crisis", size(medium)) ///
    subtitle("Indexed to 100 on March 9, 2023", size(small)) ///
    ytitle("Indexed Price (Event Date = 100)") xtitle("Date") ///
    legend(order(1 "Bitcoin" 2 "Gold" 3 "S&P 500") position(6) rows(1)) ///
    note("Vertical line: SVB collapse (Mar 9, 2023)", size(vsmall))
graph export "Figure12_SVB_Window.png", replace width(1200)
drop btc_idx gold_idx spx_idx

* --- Figure 13: Iran Conflict (February 28, 2026) ---
local ev_date = td(28feb2026)
local win_start = `ev_date' - 30

summ date, meanonly
local win_end = r(max)
local end_fmt : di %td r(max)

summ t if date >= `ev_date' & !missing(btc_price), meanonly
if r(N) > 0 {
    local btc_base = btc_price[r(min)]
    summ t if date >= `ev_date' & !missing(gold_price), meanonly
    local gold_base = gold_price[r(min)]
    summ t if date >= `ev_date' & !missing(spx_price), meanonly
    local spx_base = spx_price[r(min)]

    gen btc_idx = (btc_price / `btc_base') * 100 if date >= `win_start' & date <= `win_end'
    gen gold_idx = (gold_price / `gold_base') * 100 if date >= `win_start' & date <= `win_end'
    gen spx_idx = (spx_price / `spx_base') * 100 if date >= `win_start' & date <= `win_end'

    twoway (line btc_idx date, lcolor(orange) lwidth(medthick)) ///
           (line gold_idx date, lcolor(gold) lwidth(medthick)) ///
           (line spx_idx date, lcolor(navy) lwidth(medthick)) ///
           if date >= `win_start' & date <= `win_end', ///
        xline(`ev_date', lcolor(cranberry) lpattern(dash)) ///
        title("Asset Prices Around Iran Conflict Escalation", size(medium)) ///
        subtitle("Indexed to 100 on February 28, 2026 (ongoing)", size(small)) ///
        ytitle("Indexed Price (Event Date = 100)") xtitle("Date") ///
        legend(order(1 "Bitcoin" 2 "Gold" 3 "S&P 500") position(6) rows(1)) ///
        note("Vertical line: Iran conflict escalation (Feb 28, 2026). Data through `end_fmt'.", size(vsmall))
    graph export "Figure13_Iran2026_Window.png", replace width(1200)
    drop btc_idx gold_idx spx_idx
}
else {
    di as error "Iran event date outside data range - Figure 13 skipped."
}


*###############################################################################
*  SUMMARY OF ALL OUTPUTS
*###############################################################################

di _newline
di as text "{hline 80}"
di as text "ALL OUTPUTS GENERATED"
di as text "{hline 80}"
di as text ""
di as text "  TABLES:"
di as text "    Table1_Descriptives.rtf"
di as text "    Table2_ByRegime.rtf"
di as text "    Table3_Returns.rtf"
di as text "    Table4_Volatility.rtf"
di as text "    Table5_Event.rtf"
di as text "    Table6_TailRisk.csv"
di as text "    Table7_Drawdown.csv"
di as text "    Table7b_DrawdownDays.csv"
di as text "    Table8_ETF_Windows.rtf"
di as text "    Table9_Event60d.rtf"
di as text "    Table10_EventImmediate.rtf"
di as text ""
di as text "  FIGURES:"
di as text "    Figure1_Returns.png"
di as text "    Figure2_Volatility.png"
di as text "    Figure2_Volatility_Zoom.png"
di as text "    Figure3_GARCH.png"
di as text "    Figure4_Price.png"
di as text "    Figure5_RollingVolatility.png"
di as text "    Figure6_RollingTailRisk.png"
di as text "    Figure7_ETF_Event.png"
di as text "    Figure8_Event60d.png"
di as text "    Figure9_EventImmediate.png"
di as text "    Figure10_COVID_Window.png"
di as text "    Figure11_Ukraine_Window.png"
di as text "    Figure12_SVB_Window.png"
di as text "    Figure13_Iran2026_Window.png"
di as text "    Gold_Price_USD.png"
di as text "    Gold_RollingVolatility_USD.png"
di as text ""
di as text "  GOLD BENCHMARKS:"
di as text "    Gold_Descriptives_ByRegime.csv"
di as text "    Gold_TailRisk.csv"
di as text "    Gold_Drawdown.csv"
di as text ""
di as text "  DATA FILES:"
di as text "    Thesis_Final_Clean.dta"
di as text "    Event_Merged.dta"
di as text "    Event_Returns_Results.dta"
di as text "{hline 80}"

log close

* END OF DO-FILE

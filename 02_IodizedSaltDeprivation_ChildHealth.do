****************************************************************************
* Title: Replication .do file for nutritional health analysis
* Manuscript: "Iodized salt has large impacts on child health and test scores in Ethiopia"
* Authors: Robel Alemu, Kibrom Tafere, Dawd Gashu, Edward J. M. Joy, 
*          Elizabeth H. Bailey, R. Murray Lark, Martin R. Broadley, William A. Masters
* Date: [Date]
*
* Description:
* This Stata .do file replicates the analysis relating to the impact of 
* war-induced iodized salt disruption on nutritional health outcomes. 
* It uses data from the 2000 and 2005 waves of the Ethiopia Demographic and Health Surveys (DHS) 
* matched to district-level measures of soil micronutrient content (iodine) by unique district ID.
* 
* Steps:
* 1. Load merged DHS data with district-level soil iodine data.
* 2. Generate child growth indicators (WAZ, HAZ, WHZ) using the zscore06 command 
*    (relative to the 2006 WHO reference population).
* 3. Convert birth dates from Ethiopian to Gregorian calendar.
* 4. Define a post-war indicator (children born after 1998) to identify cohorts 
*    exposed to iodized salt withdrawal in utero or infancy.
* 5. Apply DHS survey weights and clustering for robust standard errors.
* 6. Estimate linear regression models to examine how child growth outcomes vary 
*    by exposure to withdrawal of iodized salt (post-war) and local soil iodine availability, 
*    controlling for child, parental, and household characteristics, as well as 
*    district and survey-year fixed effects.
*
* Note:
* - The data file `DHS2000_2005ChildrenRecode_SoilCropIodineWoredaAggregate.dta` 
*   should contain merged DHS child-level records with district-level average 
*   soil iodine measures, along with household and parental characteristics.
* - Variables for child's age, sex, parental education, household wealth, and 
*   birth facility status should be pre-computed in the data.
* - The code below focuses on one-year birth windows as an example 
*   (e.g., children born in 1998 vs. 1997), as described in the manuscript.
****************************************************************************

clear all
set more off

*****************************************************************************
* 1. Load merged data
*****************************************************************************
use "DHS2000_2005ChildrenRecode_SoilCropIodineWoredaAggregate.dta", clear

*****************************************************************************
* 2. Prepare Anthropometric Z-Scores
*****************************************************************************
capture drop newm newh neww
gen newm = hw15 if hw15!=9 & hw15!=. & hw15!=0
gen newh = hw3/10 if hw3!=9999 & hw3!=.
gen neww = hw2/10 if hw2!=999 & hw2!=.

zscore06, a(hw1) s(b4) h(newh) w(neww) measure(newm) male(1) female(2)

* Remove biologically implausible z-score values
replace haz06 = . if haz06<-6 | haz06>6
replace waz06 = . if waz06<-6 | waz06>5
replace whz06 = . if whz06<-5 | whz06>5

*****************************************************************************
* 3. Convert Birth Dates from Ethiopian to Gregorian Calendar
*****************************************************************************
capture drop BirthYear_GC b3_gc
gen b3_gc = b3 + 92
gen BirthYear_GC = int((b3_gc-1)/12) + 1900

*****************************************************************************
* 4. Define Post-War Indicator
*****************************************************************************
gen post_war=0 if BirthYear_GC<1998 & BirthYear_GC!=.
replace post_war=1 if BirthYear_GC>=1998 & BirthYear_GC!=.
label define postwar_label 0 "Born before 1998" 1 "Born after 1998"
label values post_war postwar_label

*****************************************************************************
* 5. Adjusting for Complex Survey Design
*****************************************************************************
replace v005 = v005/1000000
svyset [pweight=v005], psu(v021) singleunit(centered)

*****************************************************************************
* 6. Regression Analysis (Example: 1-year birth window)
*****************************************************************************
* We examine outcomes among children born in 1997 and 1998 in rural areas:
* Dependent Variables: 
*   - waz06 (Weight-for-Age Z-score)
*   - haz06 (Height-for-Age Z-score)
*   - whz06 (Weight-for-Height Z-score)
*   - BelowAvgBirthWt (binary indicator for birthweight status)
*
* Key Regressors:
*   - post_war (1 if born after 1998)
*   - lnSimulated_OrgIodine (log of local soil iodine)
* Interaction:
*   - post_war#c.lnSimulated_OrgIodine tests how the effect differs for 
*     children born after iodized salt withdrawal.
*
* Controls include child's sex, age (and powers), birth order, parental education, 
* household wealth, and facility birth. District fixed effects are absorbed via 
* absorb(), and survey-year fixed effects are included with i.dhsyear.
*****************************************************************************

encode woreda_name, gen(WoredaCode_n)

preserve

svy, subpop(if BirthYear_GC>=1997 & BirthYear_GC<=1998 & rural==1): reg waz06 ///
    i.post_war##c.lnSimulated_OrgIodine female age_years age_years2 age_years3 ///
    HhldWealthScore MatEduc_Primary PatEduc_Primary HealthFac_Delivery i.birth_order ///
    i.dhsyear, absorb(WoredaCode_n)
est store reg_waz

svy, subpop(if BirthYear_GC>=1997 & BirthYear_GC<=1998 & rural==1): reg haz06 ///
    i.post_war##c.lnSimulated_OrgIodine female age_years age_years2 age_years3 ///
    HhldWealthScore MatEduc_Primary PatEduc_Primary HealthFac_Delivery i.birth_order ///
    i.dhsyear, absorb(WoredaCode_n)
est store reg_haz

svy, subpop(if BirthYear_GC>=1997 & BirthYear_GC<=1998 & rural==1): reg whz06 ///
    i.post_war##c.lnSimulated_OrgIodine female age_years age_years2 age_years3 ///
    HhldWealthScore MatEduc_Primary PatEduc_Primary HealthFac_Delivery i.birth_order ///
    i.dhsyear, absorb(WoredaCode_n)
est store reg_whz

svy, subpop(if BirthYear_GC>=1997 & BirthYear_GC<=1998 & rural==1): reg BelowAvgBirthWt ///
    i.post_war##c.lnSimulated_OrgIodine female HealthFac_Delivery ///
    MatEduc_Primary PatEduc_Primary HhldWealthScore i.birth_order i.dhsyear, absorb(WoredaCode_n)
est store reg_bw

outreg2 [reg_waz reg_haz reg_whz reg_bw] using "HealthOutcomes_lnSimulated_OrgIodine_1YrBand.xls", ///
replace stats(coef se) alpha(0.01, 0.05, 0.10) symbol(***,**,*) dec(2) sdec(3)

restore

*****************************************************************************
* Note on Results (Table S6)
*****************************************************************************
* Supplementary Table S6 presents similar placebo tests across different birth windows:
* - Panel A: 1-year window (e.g., 1998 vs. 1997)
* - Panel B: 2-year window (1998–99 vs. 1996–97)
* - Panel C: 3-year window (1998–2000 vs. 1995–97)
*
* These outcomes serve as placebo tests, where no significant effects are expected.
*****************************************************************************


*/////////////////////////////////////////////////////////////////////////////*
*/////////////////////////////////////////////////////////////////////////////*
*/////////////////////////////////////////////////////////////////////////////*


****************************************************************************
* Analysis for Figure 5: Kaplan-Meier Cumulative Incidence of Mortality
*
* Manuscript: "Iodized salt has large impacts on child health and test scores in Ethiopia"
*
* This section of the .do file demonstrates how we compute and plot the 
* cumulative mortality (failure) probabilities (Kaplan-Meier estimates) of under-five 
* children, distinguishing by their birth cohort (pre-war vs. post-war) and 
* district-level soil nutrient concentrations (iodine and selenium).
*
* Figure 5 shows two panels:
* - Panel A: Cumulative mortality by soil iodine levels (lowest vs. highest terciles)
* - Panel B: Cumulative mortality by soil selenium levels (lowest vs. highest quartiles)
*
* Children are followed from birth until death or censoring (alive at interview).
* The curves are adjusted for child and maternal characteristics, as well as district 
* fixed effects. The "sts graph" command with the "adjustfor()" option is used to produce 
* adjusted survival curves.
****************************************************************************

*---------------------------------------------------------------------------
* Step 1: Define a binary indicator of mortality (Deseased) from the "alive" variable.
* alive = 1 if child is alive at survey, 0 if child has died.
* We create Deseased = 1 for children who died, and 0 for those alive.
*---------------------------------------------------------------------------
capture drop Deseased
gen Deseased = 1 if alive==0 & alive!=. 
replace Deseased = 0 if alive==1 & alive!=.

codebook Deseased

*---------------------------------------------------------------------------
* Step 2: (Optional) Fit a logit model to see which children are included in estimation.
* This step is just to identify the estimation sample, not strictly needed for K-M curves.
*---------------------------------------------------------------------------
set more off
logit Deseased i.post_war##c.lnSimulated_OrgIodine i.zone_name_n

capture drop in_model_Deseased
gen in_model_Deseased = e(sample)
tab in_model_Deseased

*---------------------------------------------------------------------------
* Step 3: Convert the child's age at death (b7 in months) into years for analysis time.
* If the child is alive, we will use their current age at survey as analysis time.
*---------------------------------------------------------------------------
rename b7 AgeAtDeathMonths

capture drop AgeAtDeathYrs
gen AgeAtDeathYrs = AgeAtDeathMonths/12
hist AgeAtDeathYrs
codebook AgeAtDeathYrs

*---------------------------------------------------------------------------
* Step 4: Assign each child a unique ID and define the analysis time variable.
* ChildID_Unique is needed for survival analysis with multiple children.
*---------------------------------------------------------------------------
set scheme s1color
capture drop ChildID_Unique
sort caseid v007
gen ChildID_Unique = _n
browse caseid v007 ChildID_Unique

*---------------------------------------------------------------------------
* Step 5: Define AnalysisTime as the time at risk:
* - If child died, AnalysisTime = AgeAtDeathYrs
* - If child survived, AnalysisTime = child's age at survey (age_years)
*---------------------------------------------------------------------------
capture drop AnalysisTime
gen AnalysisTime = AgeAtDeathYrs if AgeAtDeathYrs != .
replace AnalysisTime = age_years if AgeAtDeathYrs == .
browse AnalysisTime AgeAtDeathYrs age_years alive

*---------------------------------------------------------------------------
* Step 6: Declare the data as survival-time data using stset.
* fail(alive==0) means failure = death. Children who are alive are censored.
* id(ChildID_Unique) identifies individuals.
*---------------------------------------------------------------------------
stset AnalysisTime, fail(alive==0) id(ChildID_Unique)

*---------------------------------------------------------------------------
* Step 7: Set graphical style preferences for clarity.
*---------------------------------------------------------------------------
grstyle init
grstyle set plain, horizontal grid dotted

*---------------------------------------------------------------------------
* Step 8: Restrict data to children born between 1995 and 2000 in rural areas 
* (the main sample of interest).
* We also encode the district (woreda_code) if needed.
*---------------------------------------------------------------------------
preserve 
capture drop woreda_code_n
encode woreda_code, gen(woreda_code_n)
keep if BirthYear_GC>=1995 & BirthYear_GC<=2000
keep if rural ==1

*---------------------------------------------------------------------------
* Step 9: The sts list and sts test commands can be used to display 
* basic Kaplan-Meier estimates and perform tests (log-rank, Wilcoxon, etc.) 
* for differences in survival curves across groups.
*
* Example: 
* sts list, at(0(.5)5) lists estimates at intervals of 0.5 years up to 5 years.
* sts test performs statistical tests comparing survival functions by groups.
*---------------------------------------------------------------------------
sts list, at(0(.5)5)

*---------------------------------------------------------------------------
* Step 10: Run sts test with log-rank test across strata defined by soil iodine 
* categories (e.g., SimuOrgIodine_3thVs1stQuart) and post-war indicator.
*---------------------------------------------------------------------------
sts test SimuOrgIodine_3thVs1stQuart, logrank strata(post_war_exact) detail

*---------------------------------------------------------------------------
* Step 11: Generate Kaplan-Meier survival curves adjusted for covariates.
* The "sts graph" command with "fail" plots cumulative incidence (failure) probability.
* by(...) specifies grouping variables for multiple curves.
* adjustfor(...) includes covariates to adjust the curves.
*---------------------------------------------------------------------------
sts graph, fail by(SimuOrgIodine_3thVs1stQuart post_war_exact) ///
    adjustfor(female MatEduc_Primary MatAge_1stBirth i.birth_order i.zone_name_n) ///
    xtitle("Time Since Birth (Years)") ytitle("Cumulative Mortality Incidence(%)") ///
    yla(0 "0" .05 "5" .10 "10" .15 "15" .20 "20" .25 "25" .30 "30" .35 "35", angle(0)) ///
    plot1opts(lpattern(dash) lcolor(red)) plot2opts(lpattern(solid) lcolor(blue)) ///
    play("SurvivalCurve_U5.grec") ///
    saving(3_KPSurvivalCurves_SimuIodine_Over3rdQuin, replace)

restore

*---------------------------------------------------------------------------
* Step 12: Repeat similar steps for soil selenium analysis:
* Adjusting strata for SimuSolSe_3thVs1stQuart and re-running sts test and sts graph.
*---------------------------------------------------------------------------
preserve 
capture drop woreda_code_n
encode woreda_code, gen(woreda_code_n)
keep if BirthYear_GC>=1995 & BirthYear_GC<=2000
keep if rural ==1

sts list, at(0(.5)5)
sts test SimuSolSe_3thVs1stQuart, logrank strata(post_war_exact) detail

sts graph, fail by(SimuSolSe_3thVs1stQuart post_war_exact) ///
    adjustfor(female MatEduc_Primary MatAge_1stBirth i.birth_order i.zone_name_n) ///
    xtitle("Time Since Birth (Years)") ytitle("Cumulative Mortality Incidence(%)") ///
    yla(0 "0" .05 "5" .10 "10" .15 "15" .20 "20" .25 "25" .30 "30" .35 "35", angle(0)) ///
    plot1opts(lpattern(dash) lcolor(red)) plot2opts(lpattern(solid) lcolor(blue)) ///
    play("SurvivalCurve_U5.grec") ///
    saving(3_KPSurvivalCurves_SimuSe_Over3rdQuin, replace)

restore

*---------------------------------------------------------------------------
* Explanation:
* The survival analysis code:
* - Creates a binary death indicator from the alive variable.
* - Converts birth dates and ages into a consistent "AnalysisTime".
* - Declares the data as survival-time data, with stset.
* - Uses sts list and sts test to inspect and test differences in survival.
* - Uses sts graph with adjustfor() to produce adjusted Kaplan-Meier curves 
*   illustrating cumulative mortality incidence over time.
*
* Figure 5 in the manuscript uses these results to show how cumulative mortality 
* differs by birth cohort (pre-war vs. post-war) and by terciles/percentiles 
* of soil iodine or selenium.
*---------------------------------------------------------------------------

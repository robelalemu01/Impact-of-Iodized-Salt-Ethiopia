****************************************************************************
* Title: "Iodine Transmission from Soil to Grains and Humans"
* Manuscript: "Iodized salt has large impacts on child health and test scores in Ethiopia"
* Authors: Robel Alemu, Kibrom Tafere, Dawd Gashu, Edward J. M. Joy, 
*          Elizabeth H. Bailey, R. Murray Lark, Martin R. Broadley, William A. Masters
* Date: This version was revised on 10 December 2024
*
* This .do file outlines steps for generating results for:
* - Table S3: Regression of ln(crop iodine) on soil iodine measures and controls.
* - Table S4: Regression of ln(urinary iodine) on district-level ln(grain iodine) and controls.
*
* The steps below load the datasets, prepare variables, run regressions,
* and export results. All variable definitions, transformations, and controls 
* are as described in the manuscript's supplementary notes.
****************************************************************************

************************************
* Steps for Analysis for Table S3
************************************

* Step 1: Set working directory and load dataset
cls
clear all
cd "/Users/XX/ProjectDocs_SharedWithTeam/DataFiles"
use GeoNutrition_EnvtalCovariatesData, clear

* Step 2: Label variables
label var Simulated_OrgIodine "Average Simulated Organic Iodine (µg/kg)"
label var Simulated_SolubleSe "Average Simulated Soluble Se (µg/kg)"
label var Adsorbed_Se "Adsorbed_Se (µg/kg)"
label var Soluble_Se "Soluble_Se (µg/kg)"
label var Organic_Se "Organic_Se (µg/kg)"
label var lnMeanEVI_2018 "Log Enhanced Vegetation Index, 2018"
label var lnAnnualPreciptation "Log Avg. Precipitation (1979-2015)"
label var lnAnnualTemp "Log Avg. Temperature (1979-2015)"
label define graintype 1 "Field Stock or Store" 2 "Standing Crop"

* Step 3: Run regressions for Table S3 using robust regression (rreg)
preserve
estimates clear

rreg lnCropIodine lnOrganicSoilIodine i.GrainSource_n i.Crop_n
est store reg1

rreg lnCropIodine lnOrganicSoilIodine lnAltitude lnAnnualTemp lnAnnualPreciptation ///
    slopedegrees lnDisCoastObs lnMeanEVI_2018 i.GrainSource_n i.Crop_n
est store reg2

rreg lnCropIodine lnSolubleSoilIodine i.GrainSource_n i.Crop_n
est store reg3

rreg lnCropIodine lnSolubleSoilIodine lnAltitude lnAnnualTemp lnAnnualPreciptation ///
    slopedegrees lnDisCoastObs lnMeanEVI_2018 i.GrainSource_n i.Crop_n
est store reg4

rreg lnCropIodine lnAdsobedSoilIodine i.GrainSource_n i.Crop_n
est store reg5

rreg lnCropIodine lnAdsobedSoilIodine lnAltitude lnAnnualTemp lnAnnualPreciptation ///
    slopedegrees lnDisCoastObs lnMeanEVI_2018 i.GrainSource_n i.Crop_n
est store reg6

rreg lnCropIodine SoilIodine_OrgSol i.GrainSource_n i.Crop_n
est store reg7

rreg lnCropIodine SoilIodine_OrgSol lnAltitude lnAnnualTemp lnAnnualPreciptation ///
    slopedegrees lnDisCoastObs lnMeanEVI_2018 i.GrainSource_n i.Crop_n
est store reg8

* Step 4: Export results for Table S3
outreg2 [reg1 reg2 reg3 reg4 reg5 reg6 reg7 reg8] using "SoilToGrainIodine_TableS3.xls", replace ///
stats(coef se) alpha(0.01,0.05,0.10) symbol(***,**,*) dec(2) sdec(3) ///
keep(lnOrganicSoilIodine lnSolubleSoilIodine lnAdsobedSoilIodine SoilIodine_OrgSol lnAltitude c.lnAnnualTemp lnAnnualPreciptation c.lnAnnualTemp#c.lnAnnualPreciptation c.lnAltitude#c.lnAnnualPreciptation slopedegrees lnDisCoastObs lnMeanEVI_2018)

restore

****************************************************************************
* Table S3. Determinants of iodine concentration in grain samples from Amhara, Oromia and Tigray regions of Ethiopia 
*
* Note: Data shown are coefficients from OLS regression of the natural logarithm 
* of raw crop iodine concentration on the corresponding variables shown. 
* The raw soil and grain concentrations were measured from similar locations 
* and matched with site-specific covariates as described in Supplementary Note A1. 
* All regressions include controls for ten types of grain and an indicator for 
* sample source (standing crop or storage). Robust standard errors are in parentheses, 
* significance levels: *** p<0.01, ** p<0.05, * p<0.1.
****************************************************************************


************************************
* Steps for Analysis for Table S4
************************************

* Step 1: Load data linking urinary iodine to district-level grain iodine
cd "/Users/XX/ProjectDocs_SharedWithTeam/DataFiles"
use EMNS_GeoNutrition_CorrelationEnvMN_HumanMN, clear

* Step 2: Create indicator for adequate iodized salt (>15 mg/kg)
capture drop AdequateIodine
gen AdequateIodine=0 if SaltIodineStatus<2 & SaltIodineStatus!=.
replace AdequateIodine=1 if SaltIodineStatus==2 & SaltIodineStatus!=.
label define AdequateIlabel 0 "Inadequate Salt Iodine, <15ppm" 1 "Adequate Salt Iodine, >15ppm"
label values AdequateIodine AdequateIlabel

* Step 3: Compute average simulated crop iodine (Teff, Wheat, Maize)
capture drop AvgSimulated_CropIodine
egen AvgSimulated_CropIodine = rmean(Simulated_TeffIodine Simulated_WheatIodine Simulated_MaizeIodine)
label var AvgSimulated_CropIodine "Avg. Simulated Crop Iodine (Teff, Wheat & Maize), µg/kg"
destring, replace

* Step 4: Create household wealth index via PCA
capture drop HhldWealthScore
pca electricit watch radio television mobile Fixed_tele Bicycle cart Car motorcycle SolarPanel if electricit!=., components(1)
predict HhldWealthScore score
xtile Wealth_Quint= HhldWealthScore, nq(5)
label var Wealth_Quint "Wealth Quintile: 5 categories"

* Step 5: Restrict sample to rural districts of Amhara, Oromia, Tigray and keep observations with grain iodine data
preserve
estimates clear
keep if REGIONNAME=="Amhara" | REGIONNAME=="Tigray" | REGIONNAME=="Oromia"
keep if lnAvgSimulated_CropIod!=.

* Step 6: Run linear regressions for Table S4
* Different models for reproductive-age women (WRA), school-aged children (SAC), and combined samples.
* Models vary in the inclusion of age (level & squared), infection status (CRP), wealth, sex, and iodized salt controls.

* WRA
reg lnUrineIodin lnAvgSimulated_CropIod if StudyCohort=="WRA", vce(cluster EA_ID)
est store reg1
reg lnUrineIodin lnAvgSimulated_CropIod Age_Year AgeYear_sq i.ElevAGPandCRP i.Wealth_Quint if StudyCohort=="WRA", vce(cluster EA_ID)
est store reg2
reg lnUrineIodin lnAvgSimulated_CropIod Age_Year AgeYear_sq i.ElevAGPandCRP i.Wealth_Quint AdequateIodine if StudyCohort=="WRA", vce(cluster EA_ID)
est store reg3

* SAC
reg lnUrineIodin lnAvgSimulated_CropIod if StudyCohort=="SAC", vce(cluster EA_ID)
est store reg4
reg lnUrineIodin lnAvgSimulated_CropIod Age_Year AgeYear_sq female i.ElevAGPandCRP i.Wealth_Quint if StudyCohort=="SAC", vce(cluster EA_ID)
est store reg5
reg lnUrineIodin lnAvgSimulated_CropIod Age_Year AgeYear_sq female i.ElevAGPandCRP i.Wealth_Quint AdequateIodine if StudyCohort=="SAC", vce(cluster EA_ID)
est store reg6

* Combined sample
reg lnUrineIodin lnAvgSimulated_CropIod, vce(cluster EA_ID)
est store reg7
reg lnUrineIodin lnAvgSimulated_CropIod Age_Year AgeYear_sq i.ElevAGPandCRP female i.Wealth_Quint, vce(cluster EA_ID)
est store reg8
reg lnUrineIodin lnAvgSimulated_CropIod Age_Year AgeYear_sq i.ElevAGPandCRP female i.Wealth_Quint AdequateIodine, vce(cluster EA_ID)
est store reg9

* Step 7: Export results for Table S4
outreg2 [reg1 reg2 reg3 reg4 reg5 reg6 reg7 reg8 reg9] using lnUrineIodin_lnOrganicIodine.xls, ///
replace stats(coef se) alpha(0.01,0.05,0.10) symbol(***,**,*) dec(2) sdec(3) ///
keep(lnAvgSimulated_CropIod i.AdequateIodine AdequateIodine Age_Year AgeYear_sq female i.ElevAGPandCRP i.Wealth_Quint)

restore

****************************************************************************
* Table S4. Association of urinary iodine level (log) in women and children with their district average grain iodine concentration (log)
*
* Note: Data shown are coefficients from a linear regression of an individual's 
* log-transformed urinary iodine level (µg/L) on the log of district-level average grain iodine. 
* District-level grain iodine is computed using a conditional simulation approach as described 
* in Supplementary Note A.3.B, aggregating iodine levels in Maize, Wheat, and Teff.
* The sample is restricted to rural districts of Amhara, Oromia, and Tigray. Models differ 
* in the inclusion of age, infection (CRP), wealth, sex, and iodized salt use. 
* Urinary iodine data are from the 2015 ENMS. Robust standard errors, clustered at the EA level, 
* are in parentheses. Significance: *** p<0.01, ** p<0.05, * p<0.1.
****************************************************************************

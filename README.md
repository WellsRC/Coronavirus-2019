# Impact of international travel and border control measures on the global spread of the novel 2019 Coronavirus outbreak 
 
## Data

DailyIncidenceOther-WHO-Dec62019=t0.mat - Date of symptom onset of international cases with travel history to China

DailyIncidenceWuhan-NEJM-Dec62019=t0.mat - Date of symptom onset of cases in Wuhan

ReportedIncidenceOther-WHO-Dec62019=t0.mat	- Date of report of international cases with travel history to China

ReportedIncidenceWuhan-Dec62019=t0.mat - Date of report of cases in Wuhan

ReportedIncidenceHubei-Dec62019=t0.mat	- Date of report of cases in Hubei (outside of Wuhan)

ReportedIncidenceChina-Dec62019=t0.mat	- Date of report of cases in China (outside of Wuhan and Hubei) 

ReportedTimeIncidence.mat - Date of report for all international cases with travel history to China

ArrivalToSymptomOnset.mat - Day index of the arrival of a country's first case by plane

TimetoHospitalJan1onward.mat -Digitized data for time from symptom onset to hospitalization (Jan. 1 2020, onward)

TimetoHospitaluptoDec31.mat	-Digitized data for time from symptom onset to hospitalization (before Jan. 1 2020)

TimetoMedVisitJan1onward.mat	-Digitized data for time from symptom onset to first medical visit (Jan. 1 2020, onward)

TimetoMedVisituptoDec31.mat	-Digitized data for time from symptom onset to first medical visit (before Jan. 1 2020)

Weight_Flights.mat - The flight weights based on the number of airports with for all flights in/out of China and all flights in/out of China except Wuhan

Weight_Flights_Routes.mat - The flight weights based on the number of routes in/out of China and all flights in/out of China except Wuhan

IncidenceData.m - Returns incidence data for Wuhan, Hubei, China, and international cases with travel history to China. The time of symptom onset is approximated by sampling the time to first medical visit (Jan 1. onward) for cases with only date of report.

Exported_Case_Arrival_Final.xlsx - This excel file contains the arrival dates of cases into the various contries and was used to generate Table S5

## Distributions

IncubationDist.m - For a specified average duration of the incubation period, the function returns the probability distribution and S samples from the distribution. 

TimeHospitalDec31.m - Returns S samples from the distribution for the time from symptom onset to hospitalization (before Jan. 1 2020)

TimeHospitalJan1.m	- Returns S samples from the distribution for the time from symptom onset to hospitalization (Jan. 1 2020 onward)

TimeMedDec31.m	- Returns S samples from the distribution for the time from symptom onset to first medical visit (before Jan. 1 2020)

TimeMedJan1.m	- Returns S samples from the distribution for the time from symptom onset to first medical visit (Jan. 1 2020 onward)

## Fitting

CPIPT.m - Estimates the probability of travel used in the analysis using the time from symptom onset to first medical visit

CPHPT.m - Estimates the probability of travel used in the analysis using the time from symptom onset to hospitalization

CPIPTM.m - Estimates the probability of travel used in the analysis using the time from symptom onset to first medical visit for various amounts of data

## Analysis

DPCAT.m - Runs the analysis needed for Figure1 using the time from symptom onset to first medical visit

DPHAT.m - Runs the analysis needed for Figure1 using the time from symptom onset to hospitalization

WT.m - Computes the cumulative risk of the importation under a variety of travel weights (Used for Figure 2)

WTC.m - Computes the cumulative risk of the importation under the country specific travel weights

TimeAfterArival.m - Computes the time from arrival by plane to symptom onset, as well as the time from symptom onset to first transmission event

NotScreened.m - Computes the probability of a case in their incubation period travelling under various levels of contact tracing

HealthSurveyScreen.m - Computes the probability of identifying a case travelling in the incubation period based on asking the time from their last exposure.

WTR.m - Computes the cumulative risk of the importation under a variety of travel weights based on the number of routes (Used for Figure S4)

WTCR.m - Computes the cumulative risk of the importation under the country specific travel weights based on the number of routes

## Output

OutputProb.m - Outputs the results from all the analysis in text form in the command line

Figure1.m - Generates Figure 1 in the manuscript

Figure2.m - Generates Figure 2 in the manuscript

Figure3.m - Generates Figure 3B-D in the manuscript

TableTimeEstimateCountry.m - Generates Table 1 results

FigureS1.m - Generates Figure S1

FigureS2.m - Generates Figure S2

FigureS3.m - Generates Figure S3

FigureS4.m - Generates Figure S4

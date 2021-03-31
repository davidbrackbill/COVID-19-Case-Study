########################
# SIRDC model for U.S. counties
#
# Author: Hanmo Li & Mengyang Gu 
#
# Reference: Robust estimation of SARS-CoV-2 epidemic at US counties (https://arxiv.org/abs/2010.11514)
# COVID-19 dashboard: https://covid19-study.pstat.ucsb.edu/
########################

library(deSolve)
library(stringr)
library(zoo)
library(ggplot2)
library(mFilter)
library(RobustGaSP)
library(matrixStats)
library(dplyr)

source(file = "C:/Users/David/Desktop/Pstat 120c/120c Report 2/covid-19-20201128/data_and_functions/functions_covid19.R")

set.seed(1)

# select the state name, its abbreviation and the county name
state_name = "California"
state_name_short = "CA"
county_name_selected = "Imperial"

# define the start and end date
start_date_ini = as.Date("2020-03-21")
end_date = as.Date("2020-10-1")

# define the length of the prediction period
prediction_length = 21

# define the length of the training period
fitted_days_beta=180 ###we will recent 180 days to learn parameter 

# predefined the epidemiological parameters. DON'T change them unless you are confident to do so.
gamma = 0.2
theta = 0.1
delta = 0.0066

# download death toll and test data from JHU dataset and COVID-19 tracking project, respectively
base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

# download death and confirmed data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))
us_confirm = read.csv(paste0(base,confirm,"US.csv"))

covid19_project_url = "https://api.covidtracking.com/v1/states/daily.csv"

# get state level test data from COVID-19 tracking project
covid_19_project = read.csv(covid19_project_url)
covid_19_project$date = as.Date(as.character(covid_19_project$date), "%Y %m %d")

file_pow_param_path = "data_and_functions/"
file_sudden_death_path = "data_and_functions/"

# counties that have sudden death increase
county_sudden_death = read.csv(paste0(file_sudden_death_path,"counties with sudden death increase.csv" ) )
county_sudden_death$state = as.character(county_sudden_death$state)

# read the file for optim power parameter in US states
us_pow_param = read.csv(paste0(file_pow_param_path, "power parameter in US states.csv"))
us_pow_param$state_name = as.character(us_pow_param$state_name)

# clean the death data
us_death = clean_us_death(us_death)

# extract the dates in the death data
col_names = colnames(us_death)
all_dates = as.Date(paste0(str_sub(col_names[13:dim(us_death)[2]], 2, -1), "20"), tryFormats = c("%m.%d.%Y"))

# calculate the death pre million
death_rate = function(x){x/us_death$Population*10^6}
us_death_rate = us_death %>%
  dplyr::mutate_at(vars(starts_with("X")), death_rate)

n_ini = as.numeric(end_date - start_date_ini) + 1

# get the optim power parameter in US states
power_parameter = as.numeric(us_pow_param$pow_param[us_pow_param$state_name == state_name])

# get state level test data from COVID-19 tracking project
state_test = covid_19_project %>%
  dplyr::filter(state == state_name_short) %>%
  dplyr::select(date, totalTestResultsIncrease, positiveIncrease)

# clean the test data
state_test$totalTestResultsIncrease[state_test$totalTestResultsIncrease<0] = abs(state_test$totalTestResultsIncrease[state_test$totalTestResultsIncrease<0])
state_test$positiveIncrease[state_test$positiveIncrease<0] = abs(state_test$positiveIncrease[state_test$positiveIncrease<0])

state_test$daily_test_avg = rep(0, dim(state_test)[1])
state_test$daily_positive_avg = rep(0, dim(state_test)[1])

state_test$daily_test_avg = data_seven_day_smoothing(state_test$totalTestResultsIncrease)
state_test$daily_positive_avg = data_seven_day_smoothing(state_test$positiveIncrease)

state_test$PositiveRate = state_test$daily_positive_avg / state_test$daily_test_avg

# change NAs to the positive rate in the following days 
for (i in 2:(dim(state_test)[1])) {
  if(is.na(state_test$PositiveRate[i]) | is.infinite(state_test$PositiveRate[i])){
    state_test$PositiveRate[i] = state_test$PositiveRate[i-1]
  }
}

state_test$PositiveRate[is.na(state_test$PositiveRate)] = state_test$PositiveRate[which(!is.na(state_test$PositiveRate))[1]]

# change negative rate to positive
state_test$PositiveRate = abs(state_test$PositiveRate)
# change positive rate larger tha 1 to 1
state_test$PositiveRate[state_test$PositiveRate>1] = 1 / state_test$PositiveRate[state_test$PositiveRate>1]

state_death = us_death %>%
  dplyr::filter(Province_State == state_name)

state_confirmed = us_confirm %>%
  dplyr::filter(Province_State == state_name)

#########################
# Part 1: Data Cleaning
#########################

# select the county-level death and confirmed cases
county_death_raw = state_death %>%
  dplyr::filter(Admin2 == county_name_selected)

county_confirm_raw = state_confirmed %>%
  dplyr::filter(Admin2 == county_name_selected)

# the real start date is defined as the date that the state confirmed cases is no less than 5 for the first time
start_date=all_dates[which(as.numeric(county_confirm_raw[12:length(county_confirm_raw)])>=5)[1]]

if (start_date < start_date_ini){
  start_date = start_date_ini
}

# the length of real training period
n = as.numeric(end_date - start_date) + 1

# get county data from real start_date to end_date

death_with_county_names = get_output_same_time_zone(
  data_type = "death",
  state_name = state_name,
  start_date = start_date,
  state_name_short = state_name_short,
  duration = n ,
  training_length = n,
  criterion_death = 2,
  smoothness = F
)

death_with_county_names_all = get_output_same_time_zone(
  data_type = "death",
  state_name = state_name,
  start_date = start_date,
  state_name_short = state_name_short,
  duration = n+prediction_length ,
  training_length = n,
  criterion_death = 2,
  smoothness = F
)






county_names = death_with_county_names[[2]]

if(!(county_name_selected %in% county_names)){
  stop("Sorry, we cannot analysis the county you select since its death toll is less than 2.")
}

each_index = which(county_names == county_name_selected)

state_confirmed_selected = state_confirmed %>%
  dplyr::filter(Admin2 %in% county_names)

state_confirmed_selected = state_confirmed_selected[,12:dim(state_confirmed_selected)[2]]

state_confirmed_selected = state_confirmed_selected[, which(all_dates %in% seq.Date(start_date, end_date, by=1))]

SIRDC_county_population = us_death %>%
  filter(Admin2 == county_name_selected, Province_State == state_name) %>%
  dplyr::select(Population)

# county population
N = SIRDC_county_population$Population

# county death and confirmed cases from real start date to end date
death_selected = death_with_county_names[[1]][each_index, ]
confirm_selected = as.numeric(state_confirmed_selected[each_index,])

# eliminate decreasing trend in the county death data
for( i in 1:(n-1)){
  if (death_selected[i+1]<death_selected[i]){
    death_selected[i+1] = death_selected[i]
  }
}

# eliminate decreasing trend in the county confirm data
for( i in 1:(n-1)){
  if (confirm_selected[i+1]<confirm_selected[i]){
    confirm_selected[i+1] = confirm_selected[i]
  }
}

# eliminate the sudden increase of death in several prespecified counties
if(paste0(state_name_short, county_name_selected) %in%   paste0(county_sudden_death$state, county_sudden_death$county)){
  daily_death = diff(death_selected)
  sudden_increase_index = which(daily_death==max(daily_death))[1]
  daily_death[(sudden_increase_index-13):(sudden_increase_index-1)] = daily_death[(sudden_increase_index-13):(sudden_increase_index-1)] + daily_death[sudden_increase_index]/14
  daily_death[sudden_increase_index] = daily_death[sudden_increase_index]/14
  death_selected = c(death_selected[1], cumsum(daily_death) + death_selected[1])
}

# smooth the confirm cases
daily_confirm_selected = diff(confirm_selected)
daily_confirm_selected_avg = data_seven_day_smoothing(daily_confirm_selected)
unadjusted_confirm_selected_smoothed = c(confirm_selected[1], cumsum(daily_confirm_selected_avg) + confirm_selected[1])

daily_positive_rate_selected = rep(0,n)

# extract the positive rate on selected date
daily_positive_rate_selected = rev(state_test$PositiveRate[state_test$date>=(start_date) & state_test$date<=end_date])

# calculate the weights
c_t_seq =  daily_positive_rate_selected^(power_parameter)

# adjust the confirmed cases using positive_rate^{power_parameter}
daily_confirm_selected = diff(confirm_selected)
daily_adjusted_confirm_smoothed = data_seven_day_smoothing(daily_confirm_selected * c_t_seq[2:n] )
confirm_selected_smoothed = c(c_t_seq[1] * confirm_selected[1], cumsum(daily_adjusted_confirm_smoothed) + c_t_seq[1] * confirm_selected[1])

# eliminate decreasing trend in the smoothed county confirmed cases
for(i in 2:n){
  if (confirm_selected_smoothed[i]<confirm_selected_smoothed[i-1]){
    confirm_selected_smoothed[i] = confirm_selected_smoothed[i-1]
  }
}

#########################
# Part 2: Optimization
#########################

# try constrained optimization
ui = matrix(c(1,0,-1,0,1,-1), 3,2)
ci = matrix(c(0,0, -((confirm_selected_smoothed[1]*N)/confirm_selected_smoothed[n] -  death_selected[1]-0.1)  ))

I_0_ini = 1000
R_0_ini = 1000

if((I_0_ini + R_0_ini)> abs(ci[3])){
  I_0_ini = abs(ci[3])/4
  R_0_ini = abs(ci[3])/4
}

# do optimization
m_approx = constrOptim(
  c(I_0_ini, R_0_ini),
  loss_approx_beta,
  NULL,
  ui = ui,
  ci = ci,
  death_cases = death_selected,
  confirmed_cases = confirm_selected_smoothed,
  unadjusted_confirm_selected_smoothed = unadjusted_confirm_selected_smoothed,
  N_population = N,
  trianing_length = n,
  fixed_global_params = c(gamma, theta, delta),
  penalty = F,
  fitted_days= fitted_days_beta,
  weight_loss=T
)

I_0 = m_approx$par[1]
R_0 = m_approx$par[2]

#########################
# Part 3: Parameter Estimation: S_t, I_t, R_t, D_t, C_t and beta_t
#########################

ratio = confirm_selected_smoothed[1]/(I_0+ R_0+ death_selected[1])

ratio_real = confirm_selected[2]/(I_0+ R_0+ death_selected[1])

# use ratio to adjust smoothed confirmed cases and get susceptible cases S_t
estimated_confirm = confirm_selected_smoothed/ratio

# make sure the adjusted confirmed cases is at least larger than the observed confirm cases
for(i in 1:length(estimated_confirm)){
  if (estimated_confirm[i] < unadjusted_confirm_selected_smoothed[i]){
    estimated_confirm[i] = unadjusted_confirm_selected_smoothed[i]
  }
}

S_t_seq = N - estimated_confirm

# initial values of each compartment
init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)

param_record_approx_for_beta = matrix(0, 5, n) # 5 rows: S_t, I_t, R_t, D_t, C_t
param_record_approx_for_beta[,1] = init_for_beta
param_record_approx_for_beta[1,] = S_t_seq

# record the value of transmission rate
approx_beta_seq = rep(0, n-1)

# iterative approach for calculating the seq of compartments in SIRDC
for (i in 1:(n-1)){
  S_t_1 = param_record_approx_for_beta[1,i]
  S_t_2 = param_record_approx_for_beta[1,i+1]
  I_t_1 = param_record_approx_for_beta[2,i]
  R_t_1 = param_record_approx_for_beta[3,i]
  D_t_1 = param_record_approx_for_beta[4,i]
  C_t_1 = param_record_approx_for_beta[5,i]
  
  if(I_t_1<1){
    I_t_1 = 1
  }
  
  beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N, gamma=gamma)
  
  I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N) - gamma)
  R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
  D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
  C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
  
  param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
  approx_beta_seq[i] = beta_t_1_2$root
}

# calculate the smoothed transmission rate using 7 day average
approx_beta_seq_smoothed = rollapply(approx_beta_seq, width = 7, by = 1, FUN = mean, align = "left")


#########################
# Part 4: Prediction: use extended beta to calculate death prediction
#########################

param_record_approx_all = matrix(0, 5, n+prediction_length) # 5 rows: S_t, I_t, R_t, D_t, C_t
param_record_approx_all[,1] = init_for_beta

approx_beta_seq_all = c(approx_beta_seq, rep(approx_beta_seq[n-1], prediction_length))

for (i in 1:(n+prediction_length-1)){
  S_t_1 = param_record_approx_all[1,i]
  I_t_1 = param_record_approx_all[2,i]
  R_t_1 = param_record_approx_all[3,i]
  D_t_1 = param_record_approx_all[4,i]
  C_t_1 = param_record_approx_all[5,i]
  
  beta_t_1_2 = approx_beta_seq_all[i]
  
  if(I_t_1<1){
    I_t_1 = 1
  }
  
  S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma)
  I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma)
  R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
  D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
  C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
  
  param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
  
}

#########################
# Part 5: Figures
#########################

date_seq_beta = seq.Date(start_date, end_date-1, by=1)
date_seq_beta_smoothed = seq.Date(start_date+3, end_date-1-3, by=1)

y_limit_R_t = c(min(approx_beta_seq/0.2), max(approx_beta_seq/0.2))

##plot the effective reproduction number 7-day average
plot(approx_beta_seq/0.2~date_seq_beta, type="p",ylim = y_limit_R_t, ylab = "R_t", xlab = "Date", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(approx_beta_seq_smoothed/0.2~date_seq_beta_smoothed, col="blue")
abline(h=1, lty=2,col="red")
abline(h=0,lty=1,col="black")
legend("topright", legend = c("Daily Approximated R_t", "7-day Averaged Approximated R_t"), pch = c(1,NA), lty=c(NA,1), col=c("black", "blue"))

date_seq = seq.Date(start_date, end_date, by=1)

plot(N-param_record_approx_for_beta[1,]~date_seq, type="l", col="blue", xlab = "Date", ylab = "Infective Cases", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(confirm_selected~date_seq, type="l", col="red")
legend("topleft", legend = c("Estimated Confirmed Cases", "Observed Confirmed Cases"), lty = c(1,1), col = c("red", "blue"))

##this is the estimated percentage between  people infected and confirmed
confirm_selected[n]/(N-param_record_approx_for_beta[1,n])


###fitted death 
ylimit_death = c(min(param_record_approx_for_beta[4,], death_selected), max(param_record_approx_for_beta[4,], death_selected))

plot(param_record_approx_for_beta[4,]~date_seq,ylim = ylimit_death, type="l", col="blue", xlab = "Date", ylab = "Death Cases", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(death_selected~date_seq, col = "red")
legend("topleft", legend = c("Observed death toll", "Estimated death toll"), lty = c(1,1), col = c("red", "blue"))

###fitted death and forecast death 
date_seq_all = seq.Date(start_date, end_date+prediction_length, by=1)

ylimit_death_all = c(min(param_record_approx_all[4,], death_selected), max(param_record_approx_all[4,], death_selected))

plot(param_record_approx_all[4,]~date_seq_all,ylim = ylimit_death_all, type="l", col="blue", xlab = "Date", ylab = "Death Cases", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(death_selected~date_seq, col = "red")
lines(date_seq_all[(n+1):(n+prediction_length-1) ],death_with_county_names_all[[1]][each_index,(n+1):(n+prediction_length-1) ], 
      col = "red",lty=2)
abline(v = end_date+1, lty=2)
legend("topleft", legend = c("Observed death toll", "Estimated death toll","Held-out death toll"), lty = c(1,1,2), col = c("red", "blue",'red'))



###work on some simulation if the infectious period decreases
gamma_new = 1/4.75 ###suppose it changes from 5 day to 4.75 days



param_record_approx_for_beta_new = matrix(0, 5, n) # 5 rows: S_t, I_t, R_t, D_t, C_t
param_record_approx_for_beta_new[,1] = init_for_beta
param_record_approx_for_beta_new[1,] = S_t_seq

# record the value of transmission rate
# approx_beta_seq_new = rep(0, n-1)
# we should fix the beta when we change the gamma parameter

# iterative approach for calculating the seq of compartments in SIRDC
for (i in 1:(n-1)){
  S_t_1 = param_record_approx_for_beta_new[1,i]
  I_t_1 = param_record_approx_for_beta_new[2,i]
  R_t_1 = param_record_approx_for_beta_new[3,i]
  D_t_1 = param_record_approx_for_beta_new[4,i]
  C_t_1 = param_record_approx_for_beta_new[5,i]
  
  beta_t_1_2 = approx_beta_seq[i]
  
  if(I_t_1<1){
    I_t_1 = 1
  }
  
  S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_new)
  I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_new)
  R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_new/(2+theta)*(I_t_1+I_t_2)
  D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
  C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
  
  param_record_approx_for_beta_new[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
}

# calculate the smoothed transmission rate using 7 day average
# approx_beta_seq_smoothed_new = rollapply(approx_beta_seq_new, width = 7, by = 1, FUN = mean, align = "left")


###plot the simulated death 
plot(param_record_approx_for_beta_new[4,]~date_seq,ylim = ylimit_death, type="l", col="blue", xlab = "Date", ylab = "Death Cases", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(death_selected~date_seq, col = "red")
legend("topleft", legend = c("Observed death toll", "Simulated death toll"), lty = c(1,1), col = c("red", "blue"))

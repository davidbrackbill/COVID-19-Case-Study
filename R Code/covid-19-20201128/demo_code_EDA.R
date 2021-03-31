library(dplyr) ##please update it to the latest version 
library(stringr)
library(zoo)
library(ggplot2)
library(urbnmapr)

##some functions we wrote for the analysis
source(file = "data_and_functions/functions_covid19.R")



##get the data from JHU CSSE, which contain the death and confirmed cases at county-level
base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

###up-to-date US death and confirmed cases
us_death = read.csv(paste0(base,death,"US.csv"))
us_confirm = read.csv(paste0(base,confirm,"US.csv"))

##dimension
dim(us_death)
dim(us_confirm)
##get the column names of the data 
col_names = colnames(us_death)
all_dates = as.Date(paste0(str_sub(col_names[13:dim(us_death)[2]], 2, -1), "20"), tryFormats = c("%m.%d.%Y"))
##these are the dates from the data set 
all_dates

##this site provide the confirmed state-level total test from which you can get the positive rate  
covid19_project_url = "https://api.covidtracking.com/v1/states/daily.csv"

covid_19_project = read.csv(covid19_project_url)

covid_19_project$date = as.Date(as.character(covid_19_project$date), "%Y %m %d")
##note that this date starts from  the current day
covid_19_project$date


##nation level analysis
##nation level daily confirmed cases, plot and 
nation_death = us_death %>%
  dplyr::filter(Country_Region == "US")

###observed confirmed cases
nation_confirmed = us_confirm %>%
  dplyr::filter(Country_Region == "US")


###the first row of the state death data set and you can see the format
nation_death[1,]
###get the cumulative death toll and confirmed cases
nation_death_sum = apply(nation_death[,12:dim(nation_death)[2]], 2, sum)

nation_confirmed_sum = apply(nation_confirmed[,12:dim(nation_confirmed)[2]], 2, sum)

##get  dates you want to analyze  
start_date = as.Date("2020-3-21")

end_date = as.Date("2020-10-21")




##the deaths and confirmed cases for the state on the selected dates
nation_death_selected = nation_death_sum[1 + which(all_dates %in% seq.Date(start_date, end_date, by=1))]
nation_confirmed_selected = nation_confirmed_sum[which(all_dates %in% seq.Date(start_date, end_date, by=1))]

nation_death_selected=as.numeric(nation_death_selected)
nation_confirmed_selected=as.numeric(nation_confirmed_selected)

##plot cumulative confirmed cases and death
date_selected=seq.Date(start_date, end_date, by=1)
par(mfrow=c(1,2))
plot(date_selected,nation_confirmed_selected,xlab='date',ylab='cumulative observed confirmed cases',type='l')
plot(date_selected,nation_death_selected,xlab='date',ylab='cumulative death toll',type='l')
dev.off() ##close it

##daily increase between each date
daily_date_selected=date_selected[2:length(date_selected)]

##let's get the daily confirmed cases 
nation_confirmed_selected_daily=nation_confirmed_selected[2:length(nation_confirmed_selected)]-nation_confirmed_selected[1:(length(nation_confirmed_selected)-1)]
##create a data frame
daily_confirmed_nation_df = data.frame(date = daily_date_selected, value = nation_confirmed_selected_daily)

##let's get the daily death cases 
nation_death_selected_daily=nation_death_selected[2:length(nation_death_selected)]-nation_death_selected[1:(length(nation_death_selected)-1)]
##create a data frame
daily_death_nation_df = data.frame(date = daily_date_selected, value = nation_death_selected_daily)

###daily confirmed cases in US
daily_confirmed_nation_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#ff8540", width = 1) +
  ylab("Daily confirmed cases in  US")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off()
###daily death in US
daily_death_nation_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#0000FF", width = 1) +
  ylab("Daily Death cases in US")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off()

##let's obtain a seven-day average of the smoothed version of the confirmed cases and deaths

nation_confirmed_selected_daily_avg = data_seven_day_smoothing(nation_confirmed_selected_daily)

daily_confirmed_nation_smoothed_df = data.frame(date = daily_date_selected, value = nation_confirmed_selected_daily_avg)

nation_death_selected_daily_avg = data_seven_day_smoothing(nation_death_selected_daily)

daily_death_nation_smoothed_df = data.frame(date = daily_date_selected, value = nation_death_selected_daily_avg)

##plot the smoothed version 

###daily confirmed cases in US
daily_confirmed_nation_smoothed_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#ff8540", width = 1) +
  ylab("7-day averaged daily confirmed cases in  US")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off()
###daily death in US
daily_death_nation_smoothed_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#0000FF", width = 1) +
  ylab("7-day averaged daily death cases in US")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off()

##get the test positive rates 

############this is state-level positive rate
###deal with positive rate 


nation_test = covid_19_project %>%  
  dplyr::select(date, state,totalTestResultsIncrease, positiveIncrease)



nation_test_aggregated = nation_test %>%
  group_by(date) %>%
  summarise_each(funs(sum), positiveIncrease, totalTestResultsIncrease)


nation_test_aggregated$positiveIncrease_7_day_avg = data_seven_day_smoothing(nation_test_aggregated$positiveIncrease)

nation_test_aggregated$totalTestResultsIncrease_7_day_avg = data_seven_day_smoothing(nation_test_aggregated$totalTestResultsIncrease)

nation_test_aggregated$positive_rate = nation_test_aggregated$positiveIncrease_7_day_avg / nation_test_aggregated$totalTestResultsIncrease_7_day_avg


# ##reverse the sequence because it start from the current date
# nation_daily_test_selected=rev(us_test_PositiveRateus_test_PositiveRate[nation_test_aggregated$date>=(start_date) & nation_test_aggregated$date<=end_date])

###let's smooth it and get the seven day average 
us_test_daily_test_smoothed = data_seven_day_smoothing(nation_test_aggregated$totalTestResultsIncrease)
us_test_daily_positive_smoothed  = data_seven_day_smoothing(nation_test_aggregated$positiveIncrease)


##note that the following sequences start from the latest day
us_test_PositiveRate_smoothed = us_test_daily_positive_smoothed / us_test_daily_test_smoothed

us_test_PositiveRate_smoothed_selected=us_test_PositiveRate_smoothed[nation_test_aggregated$date>=(start_date) & nation_test_aggregated$date<=end_date]

###plot  the smoothed positive rates 
plot(date_selected, us_test_PositiveRate_smoothed_selected,type='l',xlab='date',ylab='-day averaged daily positive rate')


##you can also apply the above analysis  to the state-level observations


#####county level analysis 
##let's look at Santa Barbara 
state_name = "California"
state_name_short = "CA"
county_name = "Santa Barbara"
##get the death and confirmed cases
county_death = us_death%>%
  filter(Admin2 == county_name, Province_State == state_name) %>%
  select(starts_with("x"))
county_confirmed = us_confirm %>%
  filter(Admin2 == county_name, Province_State == state_name) %>%
  select(starts_with("x"))

county_death_sum = apply(county_death, 2, sum)

county_confirmed_sum = apply(county_confirmed, 2, sum)


county_death_selected = county_death_sum[which(all_dates %in% seq.Date(start_date, end_date, by=1))]
county_confirmed_selected = county_confirmed_sum[which(all_dates %in% seq.Date(start_date, end_date, by=1))]

county_death_selected=as.numeric(county_death_selected)
county_confirmed_selected=as.numeric(county_confirmed_selected)


##plot the data
##There is a jump in death on July 31, see the news: 
##https://www.ksby.com/news/coronavirus/santa-barbara-co-announces-28-previously-unreported-covid-19-related-deaths-discovered-in-data-review
par(mfrow=c(1,2))
plot(date_selected,county_death_selected,xlab='date',ylab='cumulative death toll',main=county_name,type='l')
plot(date_selected,county_confirmed_selected,xlab='date',ylab='cumulative observed confirmed cases',main=county_name,type='l')
dev.off() ##close it 

##You may plot the daily confirmed by histogram

county_confirmed_selected_daily=county_confirmed_selected[2:length(county_confirmed_selected)]-county_confirmed_selected[1:(length(county_confirmed_selected)-1)]
daily_date_selected=date_selected[2:length(date_selected)]

daily_confirmed_county_df = data.frame(date = daily_date_selected, value = county_confirmed_selected_daily)

###plot the daily confirmed cases
##you can save it as a png
#png(filename = paste0(file_path, "US_daily_confirmed_cases.png"), width = 900, height = 600)
daily_confirmed_county_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#ff8540", width = 1) +
  ylab("Daily confirmed cases in  Santa Barbara")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off() ##close it

##let's get a 7-day average of daily cases
county_confirmed_selected_daily_avg = data_seven_day_smoothing(county_confirmed_selected_daily)
daily_confirmed_county_avg_df = data.frame(date = daily_date_selected, value = county_confirmed_selected_daily_avg)

##let's plot the 7-day average of daily confirmed cases
daily_confirmed_county_avg_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#ff8540", width = 1) +
  ylab("Daily confirmed cases in Santa Barbara")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off() ##close it

##you can also plot the daily death 
county_death_selected_daily=county_death_selected[2:length(county_death_selected)]-county_death_selected[1:(length(county_death_selected)-1)]
#daily_date_selected=date_selected[2:length(date_selected)]

county_death_selected_daily_avg = data_seven_day_smoothing(county_death_selected_daily)
daily_death_county_avg_df = data.frame(date = daily_date_selected, value = county_death_selected_daily_avg)

##daily smoothed death
daily_death_county_avg_df %>%
  ggplot(aes(x=date, y=value)) +
  geom_bar(stat = 'identity', color="white",  fill="#0000FF", width = 1) +
  ylab("Daily deaths in  Santa Barbara")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))
dev.off() ##close it


###Let's see whether we can make a map for CA about the confirmed cases over county popultation
### at a particular date

# set the date for map the end date you selected before
date_for_map = end_date

##cleaning
us_confirm_death_clean = clean_JHU_data_for_map(us_confirm, us_death)

us_confirm_clean = us_confirm_death_clean[[1]]
us_death_clean = us_confirm_death_clean[[2]]

# calculate the daily confirmed cases for all US counties
us_daily_confirm_clean = us_confirm_clean
us_daily_confirm_clean[,12] = NA
us_daily_confirm_clean[,13:dim(us_confirm_clean)[2]] = t(apply(us_confirm_clean[12:dim(us_confirm_clean)[2]], 1, diff))

##you can make the confirmed cases to be zero if it is smaller than zero
#us_daily_confirm_clean[(us_daily_confirm_clean)<0] = 0

# calculate the daily death toll for all US counties
us_daily_death_clean = us_death_clean
us_daily_death_clean[,13] = NA
us_daily_death_clean[,14:dim(us_death_clean)[2]] = t(apply(us_death_clean[13:dim(us_death_clean)[2]], 1, diff))

##you can make the daily death cases to be zero if it is smaller than zero
#us_daily_death_clean[us_daily_death_clean<0] = 0

# select the daily confirm cases for CA counties
state_daily_confirm = us_daily_confirm_clean%>%
  filter(Province_State == state_name)

# This shape file contains the coordinates for county boundaries
state_shape_data = counties %>%
  filter(state_name == state_name)

# get the population for CA counties
state_population = us_daily_death_clean %>%
  filter(Province_State == state_name)%>%
  dplyr::select(Population)
state_population = as.numeric(state_population$Population)

# extract the daily confirmed cases on the selected date date_for_map
state_daily_confirm_selected = state_daily_confirm[, c(1:11, 11+ which(all_dates == date_for_map))]

# the Ratio = daily confirmed cases / population
state_daily_confirm_rate_selected = state_daily_confirm_selected
state_daily_confirm_rate_selected[,12] = state_daily_confirm_selected[,12]/state_population
state_daily_confirm_rate_selected[,12][state_daily_confirm_rate_selected[,12]==0]=NA
colnames(state_daily_confirm_rate_selected)[12] = "Ratio"

##This is  the county and rate
state_daily_confirm_rate_selected[,c(11,12)]

# joint the daily confirmed cases with the shape file
state_daily_confirm_rate_selected_joint <- left_join(state_daily_confirm_rate_selected, state_shape_data, by = "county_fips")

# find the lower and upper limits of the ratio
range(state_daily_confirm_rate_selected$Ratio*100, na.rm = T)

state_daily_confirm_rate_selected_joint %>%
  ggplot(aes(long, lat, group = group, fill = Ratio*100)) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "grey90")+
  geom_polygon(col = "black") +
  coord_map(projection = "albers", lat0 = 10, lat1 = 45) +
  labs(fill = expression("Ratio (%)")) +
  ggtitle(paste0("Daily confirmed cases / population in ", state_name, ", ", date_for_map))+
  xlab("lon") +ylab("lat" )
dev.off() ##close it
#limits

###Let's look at 7-day daily average versus county population 

# calculate the 7-day averaged daily confirm cases
state_daily_confirm_avg = state_daily_confirm
state_daily_confirm_avg[, 13:dim(state_daily_confirm)[2]] = t(apply(state_daily_confirm[,13:dim(state_daily_confirm)[2]], 1, data_seven_day_smoothing))

# extract the averaged daily confirmed cases on the selected date
state_daily_confirm_avg_selected = state_daily_confirm_avg[, c(1:11, 11+ which(all_dates == date_for_map))]

# the Ratio = averaged daily confirmed cases / population
state_daily_confirm_rate_avg_selected = state_daily_confirm_avg_selected
state_daily_confirm_rate_avg_selected[, 12] = state_daily_confirm_avg_selected[, 12]/state_population
state_daily_confirm_rate_avg_selected[,12][state_daily_confirm_rate_avg_selected[,12]==0]=NA
colnames(state_daily_confirm_rate_avg_selected)[12] = "Ratio"

## county and rate
state_daily_confirm_rate_avg_selected[,c(11,12)]

index_largest=which(state_daily_confirm_rate_avg_selected[,c(12)]==max(state_daily_confirm_rate_avg_selected[,c(12)]))
state_daily_confirm_rate_avg_selected[index_largest,c(11,12)]

# joint the averaged daily confirmed cases with the shape file
state_daily_confirm_rate_avg_selected_joint <- left_join(state_daily_confirm_rate_avg_selected, state_shape_data, by = "county_fips")


# find the lower and upper limits of the ratio
range(state_daily_confirm_rate_avg_selected$Ratio*100, na.rm = T)

state_daily_confirm_rate_avg_selected_joint %>%
  ggplot(aes(long, lat, group = group, fill = Ratio*100)) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "grey90")+
  geom_polygon(col = "black") +
  coord_map(projection = "albers", lat0 = 10, lat1 = 45) +
  labs(fill = expression("Ratio (%)")) +
  ggtitle(paste0("7-day averaged daily confirmed cases / population in ", state_name, ", ", date_for_map))+
  xlab("lon") +ylab("lat" )
dev.off() ##close it

##
###Let's see whether we can plot the whole nation (state-level)

# get the population for all US states
us_state_population = us_daily_death_clean %>%
  group_by(Province_State) %>%
  summarise(state_population = sum(Population), .groups = 'drop')

# get the daily confirm cases for all US states
us_state_daily_confirm = us_daily_confirm_clean %>%
  group_by(Province_State) %>%
  summarise(across(.fns=sum,.cols= starts_with("x")), .groups = 'drop')

# extract the daily confirmed cases for all US states on the selected date
us_state_daily_confirm_selected = us_state_daily_confirm[, c(1, 1+which(all_dates==date_for_map))]

# Ratio = 7-day averaged daily confirmed cases / population
us_state_daily_confirm_rate_selected = us_state_daily_confirm_selected
us_state_daily_confirm_rate_selected[,2] = us_state_daily_confirm_selected[,2] / us_state_population$state_population
us_state_daily_confirm_rate_selected[,2][us_state_daily_confirm_rate_selected[,2] == 0] = NA
us_state_daily_confirm_rate_selected$state_name = us_state_daily_confirm_rate_selected$Province_State
colnames(us_state_daily_confirm_rate_selected)[2] = "Ratio"

# joint the daily confirmed cases with the shape file for all US states
us_state_daily_confirmed_rate_selected_joint <- left_join(us_state_daily_confirm_rate_selected, states, by = "state_name")

range(us_state_daily_confirm_rate_selected$Ratio*100, na.rm=T)

us_state_daily_confirmed_rate_selected_joint %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (Ratio*100)), show.legend = T, col="black") +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "grey90", trans = "log10"
                      ###, limits = c(0.001, 0.125) ##you may set the limit you want to plot
                      )+
  coord_map() +
  labs(fill = expression("Ratio (%)")) +
  ggtitle(paste0("daily confirmed cases / population in the U.S.", ", ", date_for_map))
dev.off() ##close it

###Let's look at 7-day daily average versus state population (state-level)

# calculate the 7-day averaged daily confirm cases
us_state_daily_confirm_avg = us_state_daily_confirm
us_state_daily_confirm_avg[, 3:dim(us_state_daily_confirm)[2]] = t(apply(us_state_daily_confirm[, 3:dim(us_state_daily_confirm)[2]], 1, data_seven_day_smoothing))

# extract the averaged daily confirmed cases for all US states on the selected date
us_state_daily_confirm_avg_selected = us_state_daily_confirm_avg[, c(1, 1+which(all_dates==date_for_map))]

# Ratio = 7-day averaged daily confirmed cases / population
us_state_daily_confirm_rate_avg_selected = us_state_daily_confirm_avg_selected
us_state_daily_confirm_rate_avg_selected[,2] = us_state_daily_confirm_avg_selected[,2] / us_state_population$state_population
us_state_daily_confirm_rate_avg_selected[,2][us_state_daily_confirm_rate_avg_selected[,2] == 0] = NA
us_state_daily_confirm_rate_avg_selected$state_name = us_state_daily_confirm_rate_avg_selected$Province_State
colnames(us_state_daily_confirm_rate_avg_selected)[2] = "Ratio"

# joint the averaged daily confirmed cases with the shape file for all US states
us_state_daily_confirmed_rate_avg_selected_joint <- left_join(us_state_daily_confirm_rate_avg_selected, states, by = "state_name")

range(us_state_daily_confirm_rate_avg_selected$Ratio*100, na.rm=T)

us_state_daily_confirmed_rate_avg_selected_joint %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (Ratio*100)), show.legend = T, col="black") +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "grey90", trans = "log10"
                      ###, limits = c(0.001, 0.125)
                      )+
  coord_map() +
  labs(fill = expression("Ratio (%)")) +
  ggtitle(paste0("7-day averaged daily confirmed cases / population in the U.S.", ", ", date_for_map))
dev.off() ##close it

###Let's see whether we can plot the whole nation (county-level)

# extract the population for all US counties
nation_population = us_daily_death_clean$Population

# extract the daily confirmed cases for all US counties on the selected date
nation_daily_confirmed_selected = us_daily_confirm_clean[, c(1:11, 11+which(all_dates==date_for_map))]

# Ratio = daily confirmed cases / population
nation_daily_confirmed_rate_selected = nation_daily_confirmed_selected
nation_daily_confirmed_rate_selected[,12] = nation_daily_confirmed_selected[,12] / nation_population
nation_daily_confirmed_rate_selected[,12][nation_daily_confirmed_rate_selected[,12] == 0] = NA
colnames(nation_daily_confirmed_rate_selected)[12] = "Ratio"

# joint the daily confirmed cases with the shape file for all US counties
nation_daily_confirmed_rate_selected_joint <- left_join(nation_daily_confirmed_rate_selected, counties, by = "county_fips")

range(nation_daily_confirmed_rate_selected$Ratio*100, na.rm=T)

nation_daily_confirmed_rate_selected_joint %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (Ratio*100)), show.legend = T) +
  geom_polygon(
    data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),
    fill = NA, color = 'black', size = 1
  ) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "grey90", trans = "log10"
                      ###, limits = c(0.00038, 1.158)
                      )+
  coord_map() +
  labs(fill = expression("Ratio (%)")) +
  ggtitle(paste0("daily confirmed cases / population in the U.S.", ", ", date_for_map))
dev.off() ##close it

###Let's look at 7-day daily average versus county population (county-level)

# calculate the 7-day averaged daily confirm cases
us_daily_confirm_clean_avg = us_daily_confirm_clean
us_daily_confirm_clean_avg[, 13:dim(us_daily_confirm_clean)[2]] = t(apply(us_daily_confirm_clean[, 13:dim(us_daily_confirm_clean)[2]], 1, data_seven_day_smoothing))

# extract the averaged daily confirmed cases for all US counties on the selected date
nation_daily_confirmed_avg_selected = us_daily_confirm_clean_avg[, c(1:11, 11+which(all_dates==date_for_map))]

# Ratio = 7-day averaged daily confirmed cases / population
nation_daily_confirmed_rate_avg_selected = nation_daily_confirmed_avg_selected
nation_daily_confirmed_rate_avg_selected[,12] = nation_daily_confirmed_rate_avg_selected[,12] / nation_population
nation_daily_confirmed_rate_avg_selected[,12][nation_daily_confirmed_rate_avg_selected[,12] == 0] = NA
colnames(nation_daily_confirmed_rate_avg_selected)[12] = "Ratio"

# joint the averaged daily confirmed cases with the shape file for all US counties
nation_daily_confirmed_rate_avg_selected_joint <- left_join(nation_daily_confirmed_rate_avg_selected, counties, by = "county_fips")

range(nation_daily_confirmed_rate_avg_selected$Ratio*100, na.rm=T)

nation_daily_confirmed_rate_avg_selected_joint %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (Ratio*100)), show.legend = T) +
  geom_polygon(
    data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),
    fill = NA, color = 'black', size = 1
  ) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "grey90", trans = "log10"
                      ###, limits = c(0.00038, 1.158)
                      )+
  coord_map() +
  labs(fill = expression("Ratio (%)")) +
  ggtitle(paste0("7-day averaged daily confirmed cases / population in the U.S.", ", ", date_for_map))
dev.off() ##close it





library(here)
library(data.table)
library(lubridate)
library(stringr)
library(dplyr)
library(sp)
library(raster)
library(ggplot2)
library(suncalc)
library(tidyr)
library(purrr)
library(sf)
library(matlib)
library(tmap)


#Conglomeration of code for chosen figures

#read in data
Comp_temp <- read.csv(file = here::here("data", "cimis", "eco_and_cimis.csv"))

#csvs don't keep things as date objects so I need to fix that
fixdate <- function(csv){
  csv$date <- ymd(csv$date)
  csv$dt <- ymd_hms(csv$dt, tz = "America/Los_Angeles")
  csv
}

Comp_temp <- fixdate(Comp_temp)
biting_df <- fixdate(read.csv(here::here("data", "June-Sept_2018-2020", "tars_biting_df.csv")))
transmission_df <- fixdate(read.csv(here::here("data", "June-Sept_2018-2020", "tars_transmission_df.csv")))
temperature_df <- fixdate(read.csv(here::here("data", "June-Sept_2018-2020", "temperature_df.csv")))
lst_df <- fixdate(read.csv(here::here("data", "June-Sept_2018-2020", "lst_df.csv")))

sunify <- function(df){
  sunrise_sunset <- function(row){
    
    getSunlightTimes(date = ymd(row[[1]]),
                     lat = 36.7378, 
                     lon = -119.7871, #fresno coordinates
                     keep = c("sunrise", "sunset"),
                     tz = "America/Los_Angeles")
  }
  
  sunrise_set <- c()
  for(i in 1:nrow(df)){
    sunrise_set <- rbind(sunrise_set, sunrise_sunset(df[i,]))
  }
  
  df <- merge(df, sunrise_set, by = "date") #merge sunrise/sunset data 
  df
}

add_tod <- function(df, timenames = c("Dawn", "Dusk", "Day"), elsecat = "Night", ints = list(c(-1, 2), c(-1, 3), c(2, -1)), suncats = list(c("sunrise", "sunrise"), c("sunset", "sunset"), c("sunrise", "sunset"))){
  df$tod <- elsecat
  
  for(i in 1:length(timenames)){
    int <- ints[[i]]
    suncat <- suncats[[i]]
    time <- hour(df$dt) + minute(df$dt)/60
    
    df$tod <-  ifelse(time > hour(df[,suncat[1]]) + minute(df[,suncat[1]])/60 + int[1] & time < hour(df[,suncat[2]]) + minute(df[,suncat[2]])/60 + int[2], 
                      timenames[i], 
                      df$tod)
  }
  
  df
}

biting_df <- add_tod(sunify(biting_df))
transmission_df <- add_tod(sunify(transmission_df))
temperature_df <- add_tod(sunify(temperature_df))


#plotted tarsalis biting and trasnmission rates: 
T = c(seq(0,100,by=0.1))

tx_tars = - (2.94*10^-3) * T * (T - 11.3) * (T - 41.9)
bx_tars = (1.67*10^-4) * T * (T- 2.3) * (32.0 - T)^(1/2)

plot(0,col="white",xlim=c(0,55),ylim=c(-10,50),
     xlab="Temperature (C)", ylab = "")
points(T,tx_tars,pch=16,col="blue")
points(T,bx_tars*100,pch=16,col="red")
lines(T,bx_tars*100,pch=16,col="red")
legend(25, 45, legend=c("WNV transmission probability", "Biting rate x 100"),
       col=c("blue", "red"), lty = 1)
title("Culex tarsalis biting and West Nile Virus (WNV) transmission probability 
      as a function of temperature")

#visusal scatterplot of chosen model
ggplot(Comp_temp, aes(x =  ECOSTRESS, y =Air_Temp)) + 
  geom_point(aes(color = dayness, alpha = .2)) +
  # stat_function(fun = function(x){6.906e-01 + 1.432e+00*x -2.213e-02*x^2 + 1.473e-04*x^3}, show.legend = TRUE) + 
  stat_function(fun = function(x){x}, show.legend = TRUE) + 
  labs(title = "Land surface temperature vs. air temperature in the San Joaquin Valley", 
       subtitle = "June - September, 2018 - 2020", 
       caption = "Air temperature values from CIMIS weather stations.
       Land surface temperature from ECOSTRESS") + 
  ylab("Air Temperature (C)") + 
  xlab("Land Surface Temperature (C)") + 
  ylim(min = 0, max = 45) + 
  xlim(min = 0, max = 60) + 
  scale_colour_gradient2(low = "blue", high = "red") + 
  labs(color = "Hours into the day") + 
  guides(alpha = FALSE) +
  theme_dark() + 
  facet_grid(.~ifelse(vegetated, "Vegetation", "No vegetation"))

#comparative time series temperature vs lst

#day
ggplot(rbind(dplyr::select(temperature_df, names(lst_df)), lst_df), 
       aes(x = hour(dt) + minute(dt)/60 , y = mean, color = ifelse(correction == "corrected", "Modeled Air Temperature", "ECOSTRESS LST"))) + 
  geom_point(aes(alpha = .2)) +
  geom_smooth(method = "lm", formula = y ~ I(sin((2*pi*(x))/24)) + I(cos((2*pi*(x))/24)), se = FALSE) + 
  xlab("Hour of day") +
  ylab("Temperature (C)") + 
  labs(title = "ECOSTRESS LST and Modeled Air Temperature Daily Pattern", 
       subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
  guides(alpha = FALSE) + 
  ylim(min = 0, max = 50) + 
  theme(legend.title = element_blank())

#month
ggplot(rbind(dplyr::select(temperature_df, names(lst_df)), lst_df), 
       aes(x = as.factor(month(dt)), y = mean, color = ifelse(correction == "corrected", "Modeled Air Temperature", "ECOSTRESS LST"))) + 
  geom_boxplot(aes(alpha = .2)) +
  xlab("Month") +
  ylab("Temperature (C)") + 
  labs(title = "ECOSTRESS LST and Modeled Air Temperature by Month", 
       subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
  ylim(min = 0, max = 50) + 
  guides(alpha = FALSE) + 
  theme(legend.title = element_blank())

#year
ggplot(rbind(dplyr::select(temperature_df, names(lst_df)), lst_df), 
       aes(x = as.factor(year(dt)), y = mean, color = ifelse(correction == "corrected", "Modeled Air Temperature", "ECOSTRESS LST"))) + 
  geom_boxplot(aes(alpha = .2)) +
  xlab("Year") +
  ylab("Temperature (C)") + 
  labs(title = "ECOSTRESS LST and Modeled Air Temperature by Year", 
       subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
  ylim(min = 0, max = 50) + 
  guides(alpha = FALSE) + 
  theme(legend.title = element_blank())

#DAY TS

#biting rate ts
df <- base::merge(dplyr::select(biting_df, dt, mean), dplyr::select(temperature_df, dt, mean), by = "dt")
coef = (max(df$mean.x)/max(df$mean.y))
ggplot(df, 
       aes(x = hour(dt) + minute(dt)/60)) + 
  geom_point(aes(y = mean.x, alpha = .2, col = "Rate")) +
  geom_smooth(method = "lm", 
              aes(y = mean.x, col = "Rate"),
              formula = y ~ I(sin((2*pi*(x))/24)) + I(cos((2*pi*(x))/24)) + I(sin((2*2*pi*(x))/24)) + I(cos((2*2*pi*(x))/24)), 
              se = FALSE) + 
  geom_smooth(aes(y = mean.y*coef, col = "Temperature"), 
              method = "lm", 
              formula = y ~ I(sin((2*pi*(x))/24)) + I(cos((2*pi*(x))/24)), 
              se = FALSE, 
              linetype = "dashed", 
              alpha = .5) + 
  xlab("Hour of day") +
  labs(title = "Culex Tarsalis Biting Rates", 
       subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
  guides(alpha = FALSE) + 
  scale_y_continuous(
    name = "Biting Rate",
    sec.axis = sec_axis(~./coef, name="Temperature (°C)")) + 
  theme(legend.title = element_blank())

#transmission rate ts
df <- base::merge(dplyr::select(transmission_df, dt, mean), dplyr::select(temperature_df, dt, mean), by = "dt")
coef = (max(df$mean.x)/max(df$mean.y))
ggplot(df, 
       aes(x = hour(dt) + minute(dt)/60)) + 
  geom_point(aes(y = mean.x, alpha = .2, col = "Probability")) +
  geom_smooth(method = "lm", 
              aes(y = mean.x, col = "Probability"),
              formula = y ~ I(sin((2*pi*(x))/24)) + I(cos((2*pi*(x))/24)) + I(sin((2*2*pi*(x))/24)) + I(cos((2*2*pi*(x))/24)), 
              se = FALSE) + 
  geom_smooth(aes(y = mean.y*coef, col = "Temperature"), 
              method = "lm", 
              formula = y ~ I(sin((2*pi*(x))/24)) + I(cos((2*pi*(x))/24)), 
              se = FALSE, 
              linetype = "dashed", 
              alpha = .5) + 
  xlab("Hour of day") +
  labs(title = "Culex Tarsalis Transmission Probability", 
       subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
  guides(alpha = FALSE) + 
  scale_y_continuous(
    name = "Transmission Probability",
    sec.axis = sec_axis(~./coef, name="Temperature (°C)")) + 
  theme(legend.title = element_blank())

make_box_plot <- function(df, cols, newnames, title, ylab, newcols = c()){
  newdf <- pivot_longer(df, cols = all_of(cols), names_transform = newnames, names_to = "Landcover", values_to = "Probability")
  
  for(i in 1:length(newnames)){
    newdf$Landcover <- ifelse(newdf$Landcover == cols[i], newnames[i], newdf$Landcover)
  }
  newdf$Landcover <- factor(newdf$Landcover, levels = newnames)
  
  ggplot(newdf, aes(x = Landcover, y = Probability, col = Landcover)) + 
    geom_boxplot() + 
    xlab("") +
    ylab(ylab) + 
    labs(title = title, 
         subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
    scale_color_manual(values = newcols) + 
    theme(axis.text.x = element_blank())
}

# make_box_plot(biting_df, cols = c("mean", "ag", "urban", "veg", "nveg", "urb_veg", "urb_nveg", "ag_veg", "ag_nveg"), 
#               newnames = c("Whole image", "Agriculture", "Urban", "Vegetated", "Unvegetated", "Vegetated Urban", "Unvegetated Urban", "Vegetated Agriculture", "Unvegetated Agriculture"), 
#               title = "Cx. Tarsalis Biting Rates by Land Cover Type", 
#               ylab = "Biting Rate", 
#               newcols = c("black", "blue", "red", "forestgreen", "grey", "lightblue", "tomato", "darkblue", "darkred"))
# 
# make_box_plot(transmission_df, cols = c("mean", "ag", "urban", "veg", "nveg", "urb_veg", "urb_nveg", "ag_veg", "ag_nveg"), 
#               newnames = c("Whole image", "Agriculture", "Urban", "Vegetated", "Unvegetated", "Vegetated Urban", "Unvegetated Urban", "Vegetated Agriculture", "Unvegetated Agriculture"), 
#               title =  "Cx. Tarsalis WNV Transmission Probability by Land Cover Type", 
#               ylab = "Transmission Probability", 
#               newcols = c("black", "blue", "red", "forestgreen", "grey", "lightblue", "tomato", "darkblue", "darkred"))

make_box_plot(biting_df, cols = c("mean", "ag", "urban", "veg", "nveg"), 
              newnames = c("Whole image", "Agriculture", "Urban", "Vegetated", "Unvegetated"), 
              title = "Cx. Tarsalis Biting Rates by Land Cover Type", 
              ylab = "Biting Rate", 
              newcols = c("black", "blue", "red", "forestgreen", "grey"))

make_box_plot(transmission_df, cols = c("mean", "ag", "urban", "veg", "nveg"), 
              newnames = c("Whole image", "Agriculture", "Urban", "Vegetated", "Unvegetated"), 
              title =  "Cx. Tarsalis WNV Transmission Probability by Land Cover Type", 
              ylab = "Transmission Probability", 
              newcols = c("black", "blue", "red", "forestgreen", "grey"))


#CLand cover comparison time series: 

plot_tar_ts <- function(df, cols, newnames, title, ylab, newcols = c(), ts_type = "day"){
  
  newdf <- pivot_longer(df, cols = all_of(cols), names_transform = newnames, names_to = "Landcover", values_to = "Probability")
  
  for(i in 1:length(newnames)){
    newdf$Landcover <- ifelse(newdf$Landcover == cols[i], newnames[i], newdf$Landcover)
  }
  
  if(ts_type == "day"){
    
    plot <- ggplot(newdf, aes(x = hour(dt) + minute(dt)/60, y = Probability, col = Landcover)) + 
      geom_point(aes(alpha = .2)) + 
      geom_smooth(method = "lm", 
                  formula = y ~ I(sin((2*pi*(x))/24)) + I(cos((2*pi*(x))/24)) + I(sin((2*2*pi*(x))/24)) + I(cos((2*2*pi*(x))/24)), 
                  se = FALSE) + 
      xlab("Hour of day") +
      ylab(ylab) + 
      labs(title = title, 
           subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
      scale_color_manual(values = newcols) + 
      guides(alpha = FALSE)
    
  } else if(ts_type == "month"){
    
    plot <- ggplot(newdf, aes(x = month(dt) + day(dt)/30, y = Probability, col = Landcover)) + 
      geom_point(aes(alpha = .2)) + 
      xlab("Month") +
      ylab(ylab) + 
      labs(title = title, 
           subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
      guides(alpha = FALSE)
    
  }
  plot(plot)
}

plot_tar_ts(df = biting_df, 
            cols = c("ag", "urban"), 
            newnames = c("Agriculture", "Urban area"), 
            title = "Cx. Tarsalis Biting Rates in Agricultural and Urban Areas", 
            ylab = "Biting Rate", 
            newcols = c("blue", "red"), 
            ts_type = "day")

plot_tar_ts(df = biting_df, 
            cols = c("veg", "nveg"), 
            newnames = c("Vegetated", "Not vegetated"), 
            title = "Cx. Tarsalis Biting Rates in Vegetated and Unvegetated Areas", 
            ylab = "Biting Rate",
            newcols = c("grey", "forestgreen"), 
            ts_type = "day")

# plot_tar_ts(df = biting_df, 
#             cols = c("ag_veg", "ag_nveg", "urb_veg", "urb_nveg"), 
#             newnames = c("Vegetated agriculture", "Unvegetated agriculture", "Vegetated urban area", "Unvegetated urban area"), 
#             title = "Cx. Tarsalis Biting Rates in Vegetated and Unvegetated Agricultural and Urban Areas", 
#             ylab = "Biting Rate",
#             newcols = c("lightblue", "tomato", "darkblue", "darkred"),
#             ts_type = "day")

plot_tar_ts(df = transmission_df, 
            cols = c("ag", "urban"), 
            newnames = c("Agriculture", "Urban area"), 
            title = "Cx. Tarsalis WNV Transmission Probability in Agricultural and Urban Areas", 
            ylab = "Transmission Probability", 
            newcols = c("blue", "red"), 
            ts_type = "day")

plot_tar_ts(df = transmission_df, 
            cols = c("veg", "nveg"), 
            newnames = c("Vegetated", "Not vegetated"), 
            title = "Cx. Tarsalis WNV Transmission Probability in Vegetated and Unvegetated Areas", 
            ylab = "Transmission Probability",
            newcols = c("grey", "forestgreen"), 
            ts_type = "day")

# plot_tar_ts(df = transmission_df, 
#             cols = c("ag_veg", "ag_nveg", "urb_veg", "urb_nveg"), 
#             newnames = c("Vegetated agriculture", "Unvegetated agriculture", "Vegetated urban area", "Unvegetated urban area"), 
#             title = "Cx. Tarsalis WNV Transmission Probability in Vegetated and Unvegetated Agricultural and Urban Areas", 
#             ylab = "Transmission Probability",
#             newcols = c("lightblue", "tomato", "darkblue", "darkred"),
#             ts_type = "day")


make_box_plot <- function(df, cols, newnames, title, ylab, newcols = c()){
  newdf <- pivot_longer(df, cols = all_of(cols), names_transform = newnames, names_to = "Landcover", values_to = "Probability")
  
  for(i in 1:length(newnames)){
    newdf$Landcover <- ifelse(newdf$Landcover == cols[i], newnames[i], newdf$Landcover)
  }
  
  ggplot(newdf, aes(x = Landcover, y = Probability, col = Landcover)) + 
    geom_boxplot() + 
    xlab("") +
    ylab(ylab) + 
    labs(title = title, 
         subtitle = "Fresno and surrounding area, June - September, 2018 - 2020") + 
    scale_color_manual(values = newcols) + 
    theme(axis.text.x = element_blank()) + 
    facet_wrap(~tod)
}


make_box_plot(df = biting_df, 
            cols = c("ag", "urban"), 
            newnames = c("Agriculture", "Urban area"), 
            title = "Cx. Tarsalis Biting Rates in Agricultural and Urban Areas", 
            ylab = "Biting Rate", 
            newcols = c("blue", "red"))

make_box_plot(df = biting_df, 
            cols = c("veg", "nveg"), 
            newnames = c("Vegetated", "Not vegetated"), 
            title = "Cx. Tarsalis Biting Rates in Vegetated and Unvegetated Areas", 
            ylab = "Biting Rate",
            newcols = c("grey", "forestgreen"))

# make_box_plot(df = biting_df, 
#             cols = c("ag_veg", "ag_nveg", "urb_veg", "urb_nveg"), 
#             newnames = c("Vegetated agriculture", "Unvegetated agriculture", "Vegetated urban area", "Unvegetated urban area"), 
#             title = "Cx. Tarsalis Biting Rates in Vegetated and Unvegetated Agricultural and Urban Areas", 
#             ylab = "Biting Rate",
#             newcols = c("lightblue", "tomato", "darkblue", "darkred"))

make_box_plot(df = transmission_df, 
            cols = c("ag", "urban"), 
            newnames = c("Agriculture", "Urban area"), 
            title = "Cx. Tarsalis WNV Transmission Probability in Agricultural and Urban Areas", 
            ylab = "Transmission Probability", 
            newcols = c("blue", "red"))

make_box_plot(df = transmission_df, 
            cols = c("veg", "nveg"), 
            newnames = c("Vegetated", "Not vegetated"), 
            title = "Cx. Tarsalis WNV Transmission Probability in Vegetated and Unvegetated Areas", 
            ylab = "Transmission Probability",
            newcols = c("grey", "forestgreen"))

# make_box_plot(df = transmission_df, 
#             cols = c("ag_veg", "ag_nveg", "urb_veg", "urb_nveg"), 
#             newnames = c("Vegetated agriculture", "Unvegetated agriculture", "Vegetated urban area", "Unvegetated urban area"), 
#             title = "Cx. Tarsalis WNV Transmission Probability in Vegetated and Unvegetated Agricultural and Urban Areas", 
#             ylab = "Transmission Probability",
#             newcols = c("lightblue", "tomato", "darkblue", "darkred"))





Comp_temp <- filter(Comp_temp, dt %in% biting_df$dt, year(dt) == 2020)
night <- Comp_temp$dt[Comp_temp$dayness == min(Comp_temp$dayness)][1]
day <- Comp_temp$dt[Comp_temp$dayness == max(Comp_temp$dayness)][1]
dawn <- filter(Comp_temp, morning == "sunrise")$dt[
  abs(filter(Comp_temp, morning == "sunrise")$dayness) == min(abs(filter(Comp_temp, morning == "sunrise")$dayness))][1]
dusk <- filter(Comp_temp, morning == "sunset")$dt[
  abs(filter(Comp_temp, morning == "sunset")$dayness) == min(abs(filter(Comp_temp, morning == "sunset")$dayness))][1]

chosen <- c(night, dawn, day, dusk)

four_panel <- function(type, correction, limits){
  files <- list.files(here::here("data", "June-Sept_2018-2020", type, correction))
  
  date <- ymd(substr(files, 14, 23))
  hhmmss <- str_extract(files, regex('[0-9]{2}:{1}[0-9]{2}:{1}[0-9]{2}'))
  dt <- ymd_hms(paste(date, hhmmss), tz = "America/Los_Angeles")
  
  files <- files[dt %in% chosen]
  
  night <- raster(paste0(here::here("data", "June-Sept_2018-2020", type, correction), "/",
                         files[1]))
  dawn <- raster(paste0(here::here("data", "June-Sept_2018-2020", type, correction), "/",
                        files[2]))
  day <- raster(paste0(here::here("data", "June-Sept_2018-2020", type, correction), "/",
                       files[3]))
  dusk <- raster(paste0(here::here("data", "June-Sept_2018-2020", type, correction), "/",
                        files[4]))
  
  pal <- colorRampPalette(c("blue", "white", "red"))
  
  plot(dawn, 
       col = pal(50), 
       xlim=c(0,55),
       ylim=c(-10,50),
       zlim=limits,
       xlab="Longitude", 
       ylab = "Latitude")
  title(paste0("Dawn (", chosen[2],")"))
  
  plot(day, 
       col = pal(50), 
       xlim=c(0,55),
       ylim=c(-10,50),
       zlim=limits,
       xlab="Longitude", 
       ylab = "Latitude")
  title(paste0("Day (", chosen[3],")"))
  
  plot(dusk, 
       col = pal(50), 
       xlim=c(0,55),
       ylim=c(-10,50),
       zlim=limits,
       xlab="Longitude", 
       ylab = "Latitude")
  title(paste0("Dusk (", chosen[4],")"))
  
  plot(night, 
       col = pal(50), 
       xlim=c(0,55),
       ylim=c(-10,50),
       zlim=limits,
       xlab="Longitude", 
       ylab = "Latitude")
  title(paste0("Night (", chosen[1],")"))
}

four_panel("tarsalis_biting_rate", "corrected", limits = c(0,.3))
four_panel("tarsalis_transmission", "corrected", limits = c(0,20))
four_panel("temperature", "corrected", limits = c(0,70))
four_panel("temperature", "not_corrected", limits = c(0,70))


#modeled vs actual 
final <- lm(Air_Temp ~ vegetated + dayness*ECOSTRESS, 
            data = Comp_temp)
Comp_temp$air_temp_modeled <- final$fitted.values

ggplot(Comp_temp, aes(x = Air_Temp, y = air_temp_modeled)) + 
  geom_point(aes(alpha = .2)) +
  stat_function(fun = function(x){x}, aes(color = "One to one"), show.legend = TRUE) + 
  labs(title = "Modeled vs. Observed Air Temperature", 
       subtitle = "San Joaquin Valley, June - September, 2018 - 2020", 
       caption = "") + 
  ylab("Modeled Air Temperature (C)") + 
  xlab("Observed Air Temperature (C)") + 
  ylim(c(0, 45)) + 
  xlim(c(0, 45)) + 
  geom_smooth(method = "lm", se = FALSE, aes(color = "Linear regression")) +
  guides(alpha = FALSE) +
  theme(legend.title = element_blank())

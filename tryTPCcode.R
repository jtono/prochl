# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

get_model_names()

# load in data
data("chlorella_tpc")

# keep just a single curve
d <- filter(chlorella_tpc, curve_id == 1)

# show the data
ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# choose model
mod = 'sharpschoolhigh_1981'

# get start vals
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')

# get limits
low_lims <- get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')

start_vals
#>     r_tref          e         eh         th
#>  0.7485827  0.8681437  2.4861344 43.0000000
low_lims
#> r_tref      e     eh     th
#>      0      0      0      1
upper_lims
#>    r_tref         e        eh        th
#>  1.616894 10.000000 20.000000 49.000000

# fit model
fit <- nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                     data = d,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y',
                     convergence_count=FALSE)
#'convergence_count=FALSE asks nls_multstart() to try and fit all iterations of model
#'can leave out

fit



# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)
#>   rmax  topt ctmin ctmax    e    eh  q10 thermal_safety_margin
#> 1 1.81 41.65  2.54 45.56 0.58 11.48 2.06                  3.91
#>   thermal_tolerance breadth skewness
#> 1             43.02    5.37    -10.9


# predict new data
new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)

# plot data and model fit
ggplot(d, aes(temp, rate)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')



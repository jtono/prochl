### Mixed Effects Models with Gaussian ####
library(ggplot2)
library(rTPC)
library(nls.multstart)
library(purrr)
library(stargazer)
library(tidyr)
library(dplyr)
library(broom)
library(AICcmodavg)
library(emmeans)

###########function will need###########
getslopeslm <- function(data){
  treatment <- c()
  int <- c()
  sl <- c()
  rep <- c()
  for (i in unique(data$Treatment)){
    sub <- subset(data, Treatment==i)
    for (j in unique(sub$Rep.ID)){
      sub2 <- subset(sub, Rep.ID==j)
      mod <- lm(log(Actual.Cell.count)~day, sub2)
      treatment <- c(treatment, i)
      int <- c(int, coef(mod)[1])
      sl <- c(sl, coef(mod)[2])
      rep <- c(rep, j)
    }
  }
  output <- data.frame(temp=treatment, rep, int, sl)
  return(output)
}

################ load in data ####################
#read in data
b_tpc <- read.csv("Pc_B_TPC.csv", sep=",", as.is=TRUE, header=TRUE)
d_tpc <- read.csv("Pc_D_TPC.csv", sep=",", as.is=TRUE, header=TRUE)

#cut extra rows from b_tpc
b_tpc <- b_tpc[1:408,]

#get rid of ones not numbers in d_tpc and make numeric
d_tpc <- d_tpc[-which(d_tpc$Actual.Cell.count=="#VALUE!"),]
d_tpc$Actual.Cell.count <- as.numeric(d_tpc$Actual.Cell.count)

#extract date information
b_tpc$Date <- NA
for (i in 1:length(b_tpc$File.ID)){
  b_tpc$Date[i] <- strsplit(b_tpc$File.ID[i], "_")[[1]][4]
}

d_tpc$Date <- NA
for (i in 1:length(d_tpc$File.ID)){
  d_tpc$Date[i] <- strsplit(d_tpc$File.ID[i], "_")[[1]][4]
}

b_tpc$Date <- format(as.Date(b_tpc$Date, "%d%m%Y"), "20%y-%m-%d")
d_tpc$Date <- format(as.Date(d_tpc$Date, "%d%m%Y"), "20%y-%m-%d")

#convert to days since beginning
b_tpc$day <- as.numeric(round(difftime(b_tpc$Date, b_tpc$Date[1], units="days"), digits=0))
d_tpc$day <- as.numeric(round(difftime(d_tpc$Date, d_tpc$Date[1], units="days"), digits=0))

#convert treatment temps to numbers
d_tpc$Treatment <- sub("c","",d_tpc$Treatment)
d_tpc$Treatment <- as.numeric(d_tpc$Treatment)

b_tpc$Treatment <- sub("c","",b_tpc$Treatment)
b_tpc$Treatment <- as.numeric(b_tpc$Treatment)

#remove NAs from b_tpc
b_tpc <- b_tpc[-which(is.na(b_tpc$Actual.Cell.count)),]

#add log cell count column
b_tpc$ln_cell_cnt <- log(b_tpc$Actual.Cell.count)
d_tpc$ln_cell_cnt <- log(d_tpc$Actual.Cell.count)

#'get rid of last 3 points, points after day 10
b_tpc13 <- b_tpc[b_tpc$day<14,]
d_tpc13 <- d_tpc[d_tpc$day<14,]
b_tpc10 <- b_tpc[b_tpc$day<10.5,]
d_tpc10 <- d_tpc[d_tpc$day<10.5,]

#do linear regression to get data
b_tpc_sum <- getslopeslm(b_tpc)
d_tpc_sum <- getslopeslm(d_tpc)
b_tpc13_sum <- getslopeslm(b_tpc13)
d_tpc13_sum <- getslopeslm(d_tpc13)
b_tpc10_sum <- getslopeslm(b_tpc10)
d_tpc10_sum <- getslopeslm(d_tpc10)

#put in dataframe together
b_tpc_sum$strain <- "ProB"
b_tpc_sum$rep <- paste(b_tpc_sum$rep, b_tpc_sum$strain)
d_tpc_sum$strain <- "ProD"
d_tpc_sum$rep <- paste(d_tpc_sum$rep, d_tpc_sum$strain)
tpc_sum <- rbind(b_tpc_sum,d_tpc_sum)

b_tpc13_sum$strain <- "ProB"
b_tpc13_sum$rep <- paste(b_tpc13_sum$rep, b_tpc13_sum$strain)
d_tpc13_sum$strain <- "ProD"
d_tpc13_sum$rep <- paste(d_tpc13_sum$rep, d_tpc13_sum$strain)
tpc13_sum <- rbind(b_tpc13_sum,d_tpc13_sum)

b_tpc10_sum$strain <- "ProB"
b_tpc10_sum$rep <- paste(b_tpc10_sum$rep, b_tpc10_sum$strain)
d_tpc10_sum$strain <- "ProD"
d_tpc10_sum$rep <- paste(d_tpc10_sum$rep, d_tpc10_sum$strain)
tpc10_sum <- rbind(b_tpc10_sum,d_tpc10_sum)

# center the value around the rough optimum (may need to find a better value for this)
center_value = 22
tpc_sum <- tpc_sum %>%
  mutate(CenteredTemp = scale(temp, center = center_value , scale = FALSE))
tpc13_sum <- tpc13_sum %>%
  mutate(CenteredTemp = scale(temp, center = center_value , scale = FALSE))
tpc10_sum <- tpc10_sum %>%
  mutate(CenteredTemp = scale(temp, center = center_value , scale = FALSE))

require(nlme)

# Set paramters to correct form for mixed model fitting
tpc_sum$rep <- as.factor(tpc_sum$rep)
tpc_sum$strain <- as.factor(tpc_sum$strain)
tpc_sum$temp <- as.numeric(tpc_sum$temp)
tpc_sum$CenteredTemp <- as.numeric(tpc_sum$CenteredTemp)

tpc13_sum$rep <- as.factor(tpc13_sum$rep)
tpc13_sum$strain <- as.factor(tpc13_sum$strain)
tpc13_sum$temp <- as.numeric(tpc13_sum$temp)
tpc13_sum$CenteredTemp <- as.numeric(tpc13_sum$CenteredTemp)

tpc10_sum$rep <- as.factor(tpc10_sum$rep)
tpc10_sum$strain <- as.factor(tpc10_sum$strain)
tpc10_sum$temp <- as.numeric(tpc10_sum$temp)
tpc10_sum$CenteredTemp <- as.numeric(tpc10_sum$CenteredTemp)

# group the data
df_grp <- groupedData(sl ~ temp | rep, data = tpc_sum)
df_grp13 <- groupedData(sl ~ temp | rep, data = tpc13_sum)
df_grp10 <- groupedData(sl ~ temp | rep, data = tpc10_sum)

#########fit models - gnls - use simplest - gaussian_1987#########
# fit the model - with gnls cuz no random effects

#gaussian_1987
# get_start_vals(df_grp$temp, df_grp$sl, model_name = "gaussian_1987")
# # rmax       topt          a
# #0.5860897 25.0000000 15.0000000
# gau_gnls <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
#                    data = df_grp,
#                    params = list(rmax + topt + a ~ strain),
#                    # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
#                    start = c(0.5860897, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
#                    na.action = na.omit)
#Error in gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a), data = df_grp,  :
#step halving factor reduced below minimum in NLS step


get_start_vals(df_grp13$temp, df_grp13$sl, model_name = "gaussian_1987")
# rmax       topt          a
#0.6228499 25.0000000 15.0000000
gau_gnls13 <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                 data = df_grp13,
                 params = list(rmax + topt + a ~ strain),
                 # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                 start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                 na.action = na.omit)


get_start_vals(df_grp10$temp, df_grp10$sl, model_name = "gaussian_1987")
# rmax       topt          a
#0.6582842 25.0000000 15.0000000
gau_gnls10 <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                   data = df_grp10,
                   params = list(rmax + topt + a ~ strain),
                   # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                   start = c(0.6582842, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                   na.action = na.omit)
#
# #with CenteredTemp
# get_start_vals(df_grp$CenteredTemp, df_grp$sl, model_name = "gaussian_1987")
# # rmax       topt          a
# #0.5860897  3.0000000 15.0000000
# gau_gnlsC <- gnls(sl ~ gaussian_1987(temp = CenteredTemp, rmax, topt, a),
#                  data = df_grp,
#                  params = list(rmax + topt + a ~ strain),
#                  # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
#                  start = c(0.5860897, rep(0, 1), 3, rep(0, 1), 15, rep(0, 1)),
#                  na.action = na.omit)
# #Error in gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a), data = df_grp,  :
# #step halving factor reduced below minimum in NLS step
#

get_start_vals(df_grp13$CenteredTemp, df_grp13$sl, model_name = "gaussian_1987")
# rmax       topt          a
#0.6228499  3.0000000 15.0000000
gau_gnlsC13 <- gnls(sl ~ gaussian_1987(temp = CenteredTemp, rmax, topt, a),
                   data = df_grp13,
                   params = list(rmax + topt + a ~ strain),
                   # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                   start = c(0.6228499, rep(0, 1), 3, rep(0, 1), 15, rep(0, 1)),
                   na.action = na.omit)


get_start_vals(df_grp10$CenteredTemp, df_grp10$sl, model_name = "gaussian_1987")
# rmax       topt          a
#0.6582842  3.0000000 15.0000000
gau_gnlsC10 <- gnls(sl ~ gaussian_1987(temp = CenteredTemp, rmax, topt, a),
                   data = df_grp10,
                   params = list(rmax + topt + a ~ strain),
                   # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                   start = c(0.6582842, rep(0, 1), 3, rep(0, 1), 15, rep(0, 1)),
                   na.action = na.omit)

anova(gau_gnls13, gau_gnls10, gau_gnlsC13, gau_gnlsC10)
#centered temp same as not, use not

##############drop each parameter singly###########

# remove treatment effect on rmax (max rate at optimum)
gau_gnls13_2<- update(gau_gnls13, params = list(rmax ~ 1, topt + a ~ strain), start = c(0.6228499, 25, rep(0, 1), 15, rep(0, 1)))
gau_gnls10_2<- update(gau_gnls10, params = list(rmax ~ 1, topt + a ~ strain), start = c(0.6582842, 25, rep(0, 1), 15, rep(0, 1)))

# remove treatment effect on topt
gau_gnls13_3<- update(gau_gnls13, params = list(topt ~ 1, rmax + a ~ strain), start = c(25, 0.6228499, rep(0, 1), 15, rep(0, 1)))
gau_gnls10_3<- update(gau_gnls10, params = list(topt ~ 1, rmax + a ~ strain), start = c(25, 0.6582842, rep(0, 1), 15, rep(0, 1)))

# remove treatment effect on a (width)
gau_gnls13_4<- update(gau_gnls13, params = list(rmax + topt ~ strain, a ~ 1), start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 15))
gau_gnls10_4<- update(gau_gnls10, params = list(rmax + topt ~ strain, a ~ 1), start = c(0.6582842, rep(0, 1), 25, rep(0, 1), 15))

anova(gau_gnls13, gau_gnls13_2)
#dropping rmax sig worse
anova(gau_gnls13, gau_gnls13_3)
#dropping topt sig worse
anova(gau_gnls13, gau_gnls13_4)
#dropping a about the same, slightly better AIC but worse lgoLib

anova(gau_gnls10, gau_gnls10_2)
#dropping rmax sig worse
anova(gau_gnls10, gau_gnls10_3)
#dropping topt sig worse
anova(gau_gnls10, gau_gnls10_4)
#dropping a about the same, slightly better AIC but worse lgoLib

##############drop pairs of parameters###########
# remove treatment effect on rmax (max rate at optimum) and topt
gau_gnls13_5<- update(gau_gnls13, params = list(rmax + topt ~ 1, a ~ strain), start = c(0.6228499, 25, 15, rep(0, 1)))
gau_gnls10_5<- update(gau_gnls10, params = list(rmax + topt ~ 1, a ~ strain), start = c(0.6582842, 25, 15, rep(0, 1)))

# remove treatment effect on topt and a
gau_gnls13_6<- update(gau_gnls13, params = list(topt + a ~ 1, rmax ~ strain), start = c(25, 15, 0.6228499, rep(0, 1)))
gau_gnls10_6<- update(gau_gnls10, params = list(topt + a ~ 1, rmax ~ strain), start = c(25, 15, 0.6582842, rep(0, 1)))

# remove treatment effect on a (width) and rmax
gau_gnls13_7<- update(gau_gnls13, params = list(topt ~ strain, rmax + a ~ 1), start = c(25, rep(0, 1), 0.6228499, 15))
gau_gnls10_7<- update(gau_gnls10, params = list(topt ~ strain, rmax + a ~ 1), start = c(25, rep(0, 1), 0.6582842, 15))

#remove all
gau_gnls13_8<- update(gau_gnls13, params = list(topt + rmax + a ~ 1), start = c(25, 0.6228499, 15))
gau_gnls10_8<- update(gau_gnls10, params = list(topt + rmax + a ~ 1), start = c(25, 0.6582842, 15))

anova(gau_gnls13, gau_gnls13_4, gau_gnls13_5, gau_gnls13_6, gau_gnls13_7, gau_gnls13_8)
#overall either full or dropping just a best
anova(gau_gnls10, gau_gnls10_4, gau_gnls10_5, gau_gnls10_6, gau_gnls10_7, gau_gnls10_8)
#overall either full or dropping just a best

anova(gau_gnls13, gau_gnls13_5)
anova(gau_gnls13, gau_gnls13_6)
anova(gau_gnls13, gau_gnls13_7)
anova(gau_gnls13, gau_gnls13_8)
#all double drops sig worse

anova(gau_gnls10, gau_gnls10_5)
anova(gau_gnls10, gau_gnls10_6)
anova(gau_gnls10, gau_gnls10_7)
anova(gau_gnls10, gau_gnls10_8)
#all double drops sig worse

#use most reduced? - gau_gnls13_4 and gau_gnls10_4
gau_gnls_final13 <- update(gau_gnls13_4)
gau_gnls_final10 <- update(gau_gnls10_4)

summary(gau_gnls_final13)
summary(gau_gnls_final10)

#################get the diffs!!#########
comparisons_topt13 <- emmeans(gau_gnls_final13, ~ strain, param = "topt")
pairs(comparisons_topt13)
comparisons_topt10 <- emmeans(gau_gnls_final10, ~ strain, param = "topt")
pairs(comparisons_topt10)

comparisons_rmax13 <- emmeans(gau_gnls_final13, ~ strain, param = "rmax")
pairs(comparisons_rmax13)
comparisons_rmax10 <- emmeans(gau_gnls_final10, ~ strain, param = "rmax")
pairs(comparisons_rmax10)

#this is confusing??
comparisons_width_full13 <- emmeans(gau_gnls13, ~ strain, param = "a")
pairs(comparisons_width_full13)
comparisons_width_full10 <- emmeans(gau_gnls10, ~ strain, param = "a")
pairs(comparisons_width_full10)

# sig diff between all param estimated, even though model better without fitting a separately

#####################plot##########
# Create a data frame with the desired temperature levels
newX <- expand.grid(
  temp = seq(0, 40,length = 100),
  strain = unique(df_grp13$strain))

#add in topt
topt13 <- as.data.frame(comparisons_topt13)
newX$topt <- 0
newX[newX$strain=="ProB",]$topt <- topt13[topt13$strain=="ProB",]$emmean
newX[newX$strain=="ProD",]$topt <- topt13[topt13$strain=="ProD",]$emmean

# Obtain the predicted values for the new temperature levels
newX$sl <- predict(gau_gnls_final13, newdata = newX)

#plot
tpc_plot13 <-
  ggplot(df_grp13, aes(temp, sl, col=strain)) +
  geom_point() +
  geom_line(data=newX, size = 2, aes(group=strain, col=strain)) +
  geom_vline(data = newX, aes(xintercept=topt, col=strain), linetype="dashed")
tpc_plot13

################old############
# fit the model
#fit with ML so could compare to gnls
get_start_vals(df_grp$temp, df_grp$sl, model_name = "gaussian_1987")
#   rmax       topt          a
#0.6228499 25.0000000 15.0000000
gaus_nlme <- nlme(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                  data = df_grp,
                  fixed = list(rmax + topt + a ~ 1+strain),
                  random = pdDiag(rmax + topt + a ~ 1),
                  # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                  start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                  na.action = na.omit,
                  method = 'ML')
summary(gaus_nlme)


# check the best model fit and check each combination of random effects (update updates the model above without having to write it all out)
gaus_nlme2 <- update(gaus_nlme, random = pdDiag(rmax + topt ~ 1))
gaus_nlme3 <- update(gaus_nlme, random = pdDiag(rmax + a ~ 1))
gaus_nlme4 <- update(gaus_nlme, random = pdDiag(topt + a ~ 1))
gaus_nlme5 <- update(gaus_nlme, random = pdDiag(topt ~ 1))
gaus_nlme6 <- update(gaus_nlme, random = pdDiag(rmax ~ 1))
gaus_nlme7 <- update(gaus_nlme, random = pdDiag(a ~ 1))

gaus_gnls <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                  data = df_grp,
                  params = list(rmax + topt + a ~ 1 + strain),
                  #fixed = list(rmax + topt + a ~ 1),
                  #random = pdDiag(rmax + topt + a ~ 1),
                  # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                  start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                  na.action = na.omit)
summary(gaus_gnls)


# Do ANOVA style analysis to get AIC and p values on which model fits best
anova(gaus_nlme,gaus_nlme2,gaus_nlme3,gaus_nlme4,gaus_nlme5,gaus_nlme6,gaus_nlme7,gaus_gnls)
# all pretty much the same. No random effects slightly better AIC.


# remove treatment effect on rmax (max rate at optimum)
#gaus_nlme2<- update(gaus_nlme, fixed = list(rmax ~ 1, topt + a ~ 1 + strain), start = c(0.6228499, 25, rep(0, 1), 15, rep(0, 1)))
gaus_gnls2<- update(gaus_gnls, params = list(rmax ~ 1, topt + a ~ 1 + strain), start = c(0.6228499, 25, rep(0, 1), 15, rep(0, 1)))

# remove treatment effect on topt
gaus_nlme3<- update(gaus_nlme, fixed = list(topt ~ 1, rmax + a ~ 1 + strain), start = c(25, 0.6228499, rep(0,1), 15, rep(0, 1)))
gaus_gnls3<- update(gaus_gnls, params = list(topt ~ 1, rmax + a ~ 1 + strain), start = c(25, 0.6228499,rep(0,1), 15, rep(0, 1)))

# remove treatment effect on a (width)
gaus_nlme4<- update(gaus_nlme, fixed = list(rmax + topt ~ 1+strain, a ~ 1), start = c(0.6228499, rep(0,1), 25, rep(0, 1), 15))
gaus_gnls4<- update(gaus_gnls, params = list(rmax + topt ~ 1+strain, a ~ 1), start = c(0.6228499, rep(0,1), 25, rep(0, 1), 15))

anova(gaus_nlme, gaus_nlme3, gaus_nlme4, gaus_gnls, gaus_gnls2, gaus_gnls3, gaus_gnls4)
anova(gaus_nlme, gaus_nlme4, gaus_gnls, gaus_gnls4)
#quite similar with or without treatment effect on a and with or without random effects. no random and remove treatment effect on a lowest AIC but not significant
#notably, removing treatment effect on topt (model3) sig worse than full model
anova(gaus_nlme, gaus_nlme3, gaus_gnls, gaus_gnls3)
#notably, removing treatment effect on rmax (model2) sig worse than full model
anova(gaus_gnls, gaus_gnls2)


# remove treatment effect on rmax (max rate at optimum) and topt
#gaus_nlme5<- update(gaus_nlme, fixed = list(rmax +topt ~ 1, a ~ 1 + strain), start = c(0.6228499, 25, 15, rep(0, 1)))
gaus_gnls5<- update(gaus_gnls, params = list(rmax +topt ~ 1, a ~ 1 + strain), start = c(0.6228499, 25, 15, rep(0, 1)))

# remove treatment effect on topt and a
#gaus_nlme6<- update(gaus_nlme, fixed = list(topt +a ~ 1, rmax ~ 1 + strain), start = c(25, 15, 0.6228499, rep(0,1)))
gaus_gnls6<- update(gaus_gnls, params = list(topt +a~ 1, rmax ~ 1 + strain), start = c(25, 15, 0.6228499,rep(0,1)))

# remove treatment effect on a (width) and rmax
#gaus_nlme7<- update(gaus_nlme, fixed = list(topt ~ 1+strain, rmax + a ~ 1), start = c(25, rep(0, 1), 0.6228499, 15))
gaus_gnls7<- update(gaus_gnls, params = list(topt ~ 1+strain, rmax + a ~ 1), start = c(25, rep(0, 1), 0.6228499, 15))

gaus_gnls8<- update(gaus_gnls, params = list(topt + rmax + a ~ 1), start = c(25, 0.6228499, 15))

anova(gaus_gnls, gaus_gnls2, gaus_gnls3, gaus_gnls4, gaus_gnls5, gaus_gnls6, gaus_gnls7, gaus_gnls8)

anova(gaus_nlme, gaus_nlme4, gaus_gnls, gaus_gnls4, gaus_gnls5, gaus_gnls6, gaus_gnls7)
#double drops all worse
anova(gaus_nlme, gaus_gnls,gaus_gnls5, gaus_gnls6, gaus_gnls7)
anova(gaus_nlme, gaus_gnls,gaus_gnls6, gaus_gnls7)
anova(gaus_nlme, gaus_gnls,gaus_gnls7)
anova(gaus_nlme4, gaus_gnls4,gaus_gnls5, gaus_gnls6, gaus_gnls7)
anova(gaus_nlme4, gaus_gnls4,gaus_gnls6, gaus_gnls7)
anova(gaus_nlme4, gaus_gnls4,gaus_gnls7)

#use most reduced? - gaus_gnls4
gaus_gnls_final <- update(gaus_gnls4)

summary(gaus_gnls_final)

#########fit full models - gnls - only those with topt#########
# fit the model - with gnls cuz no random effects
#'for each data source, fit all models in: deutsch_2008, gaussian_1987, joehnk_2008, johnsonlewin_1946, lrf_1991, modifiedgaussian_2006, oneill_1972, pawar_2018, thomas_2012, weibull_1995
#'track AICc, AIC
fits_aicc <- c()
fits_aic <- c()
model <- c()

#deutsch_2008
get_start_vals(df_grp$temp, df_grp$sl, model_name = "deutsch_2008")
# rmax       topt      ctmax          a
#0.5860897 25.0000000 30.0000000  3.0000000
deu08_gnls <- gnls(sl ~ deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                   data = df_grp,
                   params = list(rmax + topt + ctmax + a ~ 1 + strain),
                   # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                   start = c(0.5860897, rep(0, 1), 25, rep(0, 1), 30, rep(0, 1), 3, rep(0, 1)),
                   na.action = na.omit)
#summary(deu08_gnls)
fits_aicc <- c(fits_aicc, AICc(deu08_gnls))
fits_aic <- c(fits_aic, AIC(deu08_gnls))
model <- c(model, "deu08_gnls")

#cutoff at 13
get_start_vals(df_grp13$temp, df_grp13$sl, model_name = "deutsch_2008")
# rmax       topt      ctmax          a
#0.6228499 25.0000000 30.0000000  3.0000000
deu08_gnls13 <- gnls(sl ~ deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                     data = df_grp13,
                     params = list(rmax + topt + ctmax + a ~ 1 + strain),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 30, rep(0, 1), 3, rep(0, 1)),
                     na.action = na.omit)
#summary(deu08_gnls13)
fits_aicc <- c(fits_aicc, AICc(deu08_gnls13))
fits_aic <- c(fits_aic, AIC(deu08_gnls13))
model <- c(model, "deu08_gnls13")

#cutoff at 10
get_start_vals(df_grp10$temp, df_grp10$sl, model_name = "deutsch_2008")
# rmax       topt      ctmax          a
#0.6582842 25.0000000 30.0000000  3.0000000
deu08_gnls10 <- gnls(sl ~ deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                     data = df_grp10,
                     params = list(rmax + topt + ctmax + a ~ 1 + strain),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.6582842, rep(0, 1), 25, rep(0, 1), 30, rep(0, 1), 3, rep(0, 1)),
                     na.action = na.omit)
#summary(deu08_gnls13)
fits_aicc <- c(fits_aicc, AICc(deu08_gnls10))
fits_aic <- c(fits_aic, AIC(deu08_gnls10))
model <- c(model, "deu08_gnls10")

#gaussian_1987
# get_start_vals(df_grp$temp, df_grp$sl, model_name = "gaussian_1987")
# #  rmax       topt          a
# #0.5860897 25.0000000 15.0000000
# gau87_gnls <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
#                    data = df_grp,
#                    params = list(rmax + topt + a ~ 1 + strain),
#                    # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
#                    start = c(0.5860897, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
#                    na.action = na.omit)
# fits_aicc <- c(fits_aicc, AICc(gau87_gnls))
# fits_aic <- c(fits_aic, AIC(gau87_gnls))
# model <- c(model, "gau87_gnls")
#Error in gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a), data = df_grp,  :
#step halving factor reduced below minimum in NLS step

#gaussian_1987 - cutoff at 13
get_start_vals(df_grp13$temp, df_grp13$sl, model_name = "gaussian_1987")
#  rmax       topt          a
#0.6228499 25.0000000 15.0000000
gau87_gnls13 <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                     data = df_grp13,
                     params = list(rmax + topt + a ~ 1 + strain),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                     na.action = na.omit)
fits_aicc <- c(fits_aicc, AICc(gau87_gnls13))
fits_aic <- c(fits_aic, AIC(gau87_gnls13))
model <- c(model, "gau87_gnls13")

#gaussian_1987 - cutoff at 10
get_start_vals(df_grp10$temp, df_grp10$sl, model_name = "gaussian_1987")
#  rmax       topt          a
#0.6582842 25.0000000 15.0000000
gau87_gnls10 <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                     data = df_grp10,
                     params = list(rmax + topt + a ~ 1 + strain),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.6582842, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                     na.action = na.omit)
fits_aicc <- c(fits_aicc, AICc(gau87_gnls10))
fits_aic <- c(fits_aic, AIC(gau87_gnls10))
model <- c(model, "gau87_gnls10")

# #joehnk_2008
# get_start_vals(df_grp$temp, df_grp$sl, model_name = "joehnk_2008")
# # rmax       topt          a          b          c
# #0.5860897 25.0000000  9.6866667  1.1133333  1.1866667
# joe08_gnls <- gnls(sl ~ joehnk_2008(temp = temp, rmax, topt, a, b, c),
#                      data = df_grp,
#                      params = list(rmax + topt + a + b + c ~ 1 + strain),
#                      # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
#                      start = c(0.5860897, rep(0, 1), 25, rep(0, 1), 9.6866667, rep(0, 1), 1.1133333, rep(0, 1), 1.1866667, rep(0, 1)),
#                      na.action = na.omit)
# fits_aicc <- c(fits_aicc, AICc(joe08_gnls))
# fits_aic <- c(fits_aic, AIC(joe08_gnls))
# model <- c(model, "joe08_gnls")
#Error in gnls(sl ~ joehnk_2008(temp = temp, rmax, topt, a, b, c), data = df_grp,  :
#step halving factor reduced below minimum in NLS step
#In addition: There were 14 warnings (use warnings() to see them)

# #joehnk_2008 - cutoff 13
# get_start_vals(df_grp13$temp, df_grp13$sl, model_name = "joehnk_2008")
# # rmax       topt          a          b          c
# #0.6228499 25.0000000  9.6866667  1.1133333  1.1866667
# joe08_gnls13 <- gnls(sl ~ joehnk_2008(temp = temp, rmax, topt, a, b, c),
#                    data = df_grp13,
#                    params = list(rmax + topt + a + b + c ~ 1 + strain),
#                    # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
#                    start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 9.6866667, rep(0, 1), 1.1133333, rep(0, 1), 1.1866667, rep(0, 1)),
#                    na.action = na.omit)
# fits_aicc <- c(fits_aicc, AICc(joe08_gnls13))
# fits_aic <- c(fits_aic, AIC(joe08_gnls13))
# model <- c(model, "joe08_gnls13")
# #[1] "approximate covariance matrix for parameter estimates not of full rank"

# #joehnk_2008 - cutoff 10
# get_start_vals(df_grp10$temp, df_grp10$sl, model_name = "joehnk_2008")
# # rmax       topt          a          b          c
# #0.6582842 25.0000000  9.6866667  1.1133333  1.1866667
# joe08_gnls10 <- gnls(sl ~ joehnk_2008(temp = temp, rmax, topt, a, b, c),
#                      data = df_grp10,
#                      params = list(rmax + topt + a + b + c ~ 1 + strain),
#                      # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
#                      start = c(0.6582842, rep(0, 1), 25, rep(0, 1), 9.6866667, rep(0, 1), 1.1133333, rep(0, 1), 1.1866667, rep(0, 1)),
#                      na.action = na.omit)
# fits_aicc <- c(fits_aicc, AICc(joe08_gnls10))
# fits_aic <- c(fits_aic, AIC(joe08_gnls10))
# model <- c(model, "joe08_gnls10")
# #1] "approximate covariance matrix for parameter estimates not of full rank"

#joehnk_2008 didn't run at all

#johnsonlewin_1946
get_start_vals(df_grp$temp, df_grp$sl, model_name = "johnsonlewin_1946")
#        r0           e          eh        topt
#-0.01809075  1.51895649  5.10506243 25.00000000
joh46_gnls <- gnls(sl ~ johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                   data = df_grp13,
                   params = list(r0 + e + eh + topt ~ 1 + strain),
                   # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                   start = c(-0.01809075, rep(0, 1), 1.51895649, rep(0, 1), 5.10506243, rep(0, 1), 25, rep(0, 1)),
                   na.action = na.omit)
#summary(deu08_gnls)
fits_aicc <- c(fits_aicc, AICc(joh46_gnls))
fits_aic <- c(fits_aic, AIC(joh46_gnls))
model <- c(model, "joh46_gnls")
#Error in gnls(sl ~ johnsonlewin_1946(temp = temp, r0, e, eh, topt), data = df_grp,  :
#step halving factor reduced below minimum in NLS step

#####lots of problems fitting these models. Just go with best fit individually.

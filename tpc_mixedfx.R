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

#'get rid of last 3 points
b_tpc <- b_tpc[b_tpc$day<14,]
d_tpc <- d_tpc[d_tpc$day<14,]

#do linear regression to get data
b_tpc_sum <- getslopeslm(b_tpc)
d_tpc_sum <- getslopeslm(d_tpc)

#put in dataframe together
b_tpc_sum$strain <- "ProB"
b_tpc_sum$rep <- paste(b_tpc_sum$rep, b_tpc_sum$strain)
d_tpc_sum$strain <- "ProD"
d_tpc_sum$rep <- paste(d_tpc_sum$rep, d_tpc_sum$strain)
tpc_sum <- rbind(b_tpc_sum,d_tpc_sum)


# center the value around the rough optimum (may need to find a better value for this)
center_value_b = 23
center_value_d = 21
b_tpc_sum <- b_tpc_sum %>%
  mutate(CenteredTemp = scale(temp, center = center_value_b , scale = FALSE))
d_tpc_sum <- d_tpc_sum %>%
  mutate(CenteredTemp = scale(temp, center = center_value_d , scale = FALSE))
center_value = 22
tpc_sum <- tpc_sum %>%
  mutate(CenteredTemp = scale(temp, center = center_value , scale = FALSE))

require(nlme)

# Set paramters to correct form for mixed model fitting
b_tpc_sum$rep <- as.factor(b_tpc_sum$rep)
b_tpc_sum$temp <- as.numeric(b_tpc_sum$temp)
b_tpc_sum$CenteredTemp <- as.numeric(b_tpc_sum$CenteredTemp)
d_tpc_sum$rep <- as.factor(d_tpc_sum$rep)
d_tpc_sum$temp <- as.numeric(d_tpc_sum$temp)
d_tpc_sum$CenteredTemp <- as.numeric(d_tpc_sum$CenteredTemp)

tpc_sum$rep <- as.factor(tpc_sum$rep)
tpc_sum$strain <- as.factor(tpc_sum$strain)
tpc_sum$temp <- as.numeric(tpc_sum$temp)
tpc_sum$CenteredTemp <- as.numeric(tpc_sum$CenteredTemp)

# group the data
df_grp_b <- groupedData(sl ~ temp | rep, data = b_tpc_sum)
df_grp_b <- groupedData(sl ~ CenteredTemp | rep, data = b_tpc_sum)
df_grp_d <- groupedData(sl ~ temp | rep, data = d_tpc_sum)
df_grp_d <- groupedData(sl ~ CenteredTemp | rep, data = d_tpc_sum)
df_grp <- groupedData(sl ~ temp | rep, data = tpc_sum)
df_grp <- groupedData(sl ~ CenteredTemp | rep, data = tpc_sum)

# fit the model
#fit with ML so could compare to gnls
get_start_vals(tpc_sum$temp, tpc_sum$sl, model_name = "gaussian_1987")
#   rmax       topt          a
#0.6228499 25.0000000 15.0000000
gaus_nlme <- nlme(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                  data = df_grp,
                  fixed = list(rmax + topt + a ~ 1+strain),
                  random = pdDiag(rmax + topt + a ~ 1),
                  # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                  start = c(0.6228499, rep(0, 1), 25, rep(0, 1), 15, rep(0, 1)),
                  na.action = na.omit,
                  method = 'REML')
summary(gaus_nlme)

#########left off here###########


# check the best model fit and check each combination of random effects (update updates the model above without having to write it all out)
gaus_nlme_b2 <- update(gaus_nlme_b, random = pdDiag(rmax + topt ~ 1))
gaus_nlme_b3 <- update(gaus_nlme_b, random = pdDiag(rmax + a ~ 1))
gaus_nlme_b4 <- update(gaus_nlme_b, random = pdDiag(topt + a ~ 1))
gaus_nlme_b5 <- update(gaus_nlme_b, random = pdDiag(topt ~ 1))
gaus_nlme_b6 <- update(gaus_nlme_b, random = pdDiag(rmax ~ 1))
gaus_nlme_b7 <- update(gaus_nlme_b, random = pdDiag(a ~ 1))

gaus_gnls_b8 <- gnls(sl ~ gaussian_1987(temp = temp, rmax, topt, a),
                    data = df_grp_b,
                    #fixed = list(rmax + topt + a ~ 1),
                    #random = pdDiag(rmax + topt + a ~ 1),
                    # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                    start = get_start_vals(b_tpc_sum$temp, b_tpc_sum$sl, model_name = "gaussian_1987"),
                    na.action = na.omit)
summary(gaus_gnls_b8)


# Do ANOVA style analysis to get AIC and p values on which model fits best
anova(gaus_nlme_b,gaus_nlme_b2,gaus_nlme_b3,gaus_nlme_b4,gaus_nlme_b5,gaus_nlme_b6,gaus_nlme_b7,gaus_gnls_b8)
# all pretty much the same. No random effects slightly better AIC.
########left off here - add more comparisons from Hebe's code with strain##########


#################

comparisons_height <- emmeans(quad_temp_nlme_final, ~ Stress, param = "a")

pairs(comparisons_height)

comparisons_breadth <- emmeans(quad_temp_nlme_final, ~ Stress, param = "c")

pairs(comparisons_breadth)

# we see that there are significant differences between the height of all curves other than between pH and salinity
# we see that there is a signifcaint difference between the breadth of the control and the combined stressor trearment
# but no differences between other treatments.

#####################
# Create a data frame with the desired temperature levels
newX <- expand.grid(
  CenteredTemp= seq(-12, 20,length = 100),
  Stress= unique(d_temp$Stress),
  Id= unique(d_temp$Id))

# Obtain the predicted values for the new temperature levels
newX$estimate <- predict(quad_temp_nlme_final, newdata = newX, level = 0)

Temp_nlme_Plot <- ggplot(d_temp, aes(CenteredTemp, estimate)) +
  geom_point() +
  geom_line(data=newX, size = 2) +
  facet_grid(~Stress)
Temp_nlme_Plot

################ pH mixed effects model ####################
#### pH ####
d_pH <-  read.csv("/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/pH_Aug22/LogisticGrowth_pH_Aug22/LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv")
d_pH  <- d_pH  %>% mutate(Stress = case_when((SalinityLevel == 0 & TempLevel == 20) ~ "Control",
                                             (SalinityLevel == 0  & TempLevel == 38) ~ "High Temperature",
                                             (SalinityLevel == 20 & TempLevel == 20) ~ "High Salinity",
                                             (SalinityLevel == 20 & TempLevel == 38) ~ "High Salinity x High Temperature"))

# Set paramters to correct form for mixed model fitting
d_pH$Rep <- as.factor(d_pH$Rep)
d_pH$Id <- as.factor(d_pH$Id)
d_pH$Stress <- as.factor(d_pH$Stress)
d_pH$pH_Level <- as.numeric(d_pH$pH_Level)

# group the data
df_grp <- groupedData(estimate ~ pH_Level | Id, data = d_pH)

# fit the model
gaus_pH_nlme <- nlme(estimate ~ gaussian_1987(temp = pH_Level, rmax, topt, a),
                     data = df_grp,
                     fixed = list(rmax + topt + a ~ 1 + Stress),
                     random = pdDiag(rmax + topt + a ~ 1),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.004, rep(0, 3), 7.2, rep(0, 3), 2,rep(0, 3)),
                     na.action = na.omit,
                     method = 'REML',
                     control = nlmeControl(pnlsTol = 0.02))

# fit the model
quad_pH_nlme <- nlme(estimate ~ quadratic_2008(temp = pH_Level, a, b, c),
                     data = df_grp,
                     fixed = list(a + b+ c ~ 1 + Stress),
                     random = pdDiag(a + b+ c ~ 1),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.001, rep(0, 3), 1, rep(0, 3), 1,rep(0, 3)),
                     na.action = na.omit,
                     method = 'REML',
                     control = nlmeControl(pnlsTol = 0.02))
summary(quad_pH_nlme)

# Fit the model using glmmTMB
gaus_pH_glmmTMB <- glmmTMB(estimate ~ gaussian_1987(temp = pH_Level, rmax, topt, a),
                           data = df_grp,  # Specify the correct dataframe
                           family = Gamma(link = "log"),
                           start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 8, rep(0, 3)),
                           control = glmmTMBControl(pnlsTol = 0.02))
summary(gaus_pH_nlme)

gaus_pH_nlme_A <- update(gaus_pH_nlme, start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 2.2,rep(0, 3)))
gaus_pH_nlme_B <- update(gaus_pH_nlme, start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 3,rep(0, 3)))
gaus_pH_nlme_C <- update(gaus_pH_nlme, start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 2.5,rep(0, 3)))
gaus_pH_nlme_D <- update(gaus_pH_nlme, start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 2,rep(0, 3)))
gaus_pH_nlme_E <- update(gaus_pH_nlme, start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 8,rep(0, 3)))

anova(gaus_pH_nlme, gaus_pH_nlme_B, gaus_pH_nlme_A, gaus_pH_nlme_C, gaus_pH_nlme_D, gaus_pH_nlme_E)


# check the best model fit and check each combination of random effects (update updates the model above without having to write it all out)
gaus_pH_nlme_2 <- update(gaus_pH_nlme, random = pdDiag(rmax + topt ~ 1)) # failed
gaus_pH_nlme_3 <- update(gaus_pH_nlme, random = pdDiag(rmax + a ~ 1)) # failed
gaus_pH_nlme_4 <- update(gaus_pH_nlme, random = pdDiag(a + topt ~ 1))# failed
gaus_pH_nlme_5 <- update(gaus_pH_nlme, random = pdDiag(rmax ~ 1))
gaus_pH_nlme_6 <- update(gaus_pH_nlme, random = pdDiag(topt ~ 1))
gaus_pH_nlme_7 <- update(gaus_pH_nlme, random = pdDiag(a ~ 1))


# Do ANOVA style analysis to get AIC and p values on which model fits best
anova(gaus_pH_nlme, gaus_pH_nlme_6, gaus_pH_nlme_5, gaus_pH_nlme_7) # full random effects model is best


# fit the model with ML to look at fixed effects
gaus_pH_nlme <- nlme(estimate ~ gaussian_1987(temp = pH_Level, rmax, topt, a),
                     data = df_grp,
                     fixed = list(rmax + topt + a ~ 1 + Stress),
                     random = pdDiag(rmax + topt + a ~ 1),
                     # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                     start = c(0.006, rep(0, 3),7.2, rep(0, 3), 8, rep(0, 3)),
                     na.action = na.omit,
                     method = 'ML',
                     control = nlmeControl(pnlsTol = 0.02))

# remove treatment effect on rmax
gaus_pH_nlme1<- update(gaus_pH_nlme, fixed = list(rmax ~ 1, topt + a ~ 1 + Stress), start = c(0.006, 7.2, rep(0, 3), 8, rep(0, 3)))

# remove treatment effect on topt
gaus_pH_nlme2<- update(gaus_pH_nlme, fixed = list(topt ~ 1, rmax + a ~ 1 + Stress), start = c(7.2, 0.006, rep(0, 3), 8, rep(0, 3))) #failed

# remove treatment effect on breadth (a)
gaus_pH_nlme3<- update(gaus_pH_nlme, fixed = list(rmax + topt ~ 1 + Stress, a ~ 1), start = c(0.006, rep(0, 3), 7.2, rep(0, 3), 8)) #failed


anova(gaus_pH_nlme, gaus_pH_nlme1) # full model best

# now remove more and see the effects and which is best model
# remove treatment effect on rmax	 and topt
gaus_pH_nlme4<- update(gaus_pH_nlme, fixed = list(rmax + topt ~ 1, a ~ 1 + Stress), start = c(0.006, 7.2, 8, rep(0, 3))) #failed

# remove treatment effect on topt and a
gaus_pH_nlme5<- update(gaus_pH_nlme, fixed = list(topt ~ 1 + Stress, rmax + a  ~ 1), start = c(0.006, 7.2, rep(0, 3), 8 )) # failed

# remove treatment effect on breadth and topt
gaus_pH_nlme6<- update(gaus_pH_nlme, fixed = list(rmax ~ 1 + Stress, topt + a ~ 1 ), start = c(0.006, rep(0, 3), 7.2, 8))

anova(gaus_pH_nlme, gaus_pH_nlme6) # full model best here

# full model best
gaus_pH_nlme_final <- update(gaus_pH_nlme, method = 'REML')

#################


#####################
# Create a data frame with the desired pHerature levels
newX <- expand.grid(
  pH_Level= seq(4, 11,length = 100),
  Stress= unique(d_pH$Stress),
  Id= unique(d_pH$Id))

# Obtain the predicted values for the new pHerature levels
newX$estimate <- predict(quad_pH_nlme, newdata = newX, level = 0)

pH_nlme_Plot <- ggplot(d_pH, aes(pH_Level, estimate)) +
  geom_point() +
  geom_line(data=newX, size = 2) +
  facet_grid(~Stress)
pH_nlme_Plot

################ Sal mixed effects model ####################
#### Sal ####
d_Sal <-  read.csv("/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/Sal_Aug22/LogisticGrowth_Sal_Aug22/LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv")

d_Sal  <- d_Sal  %>% mutate(Stress = case_when((pH_Level == 7.2 & TempLevel == 20) ~ "Control",
                                               (pH_Level == 7.2 & TempLevel == 38) ~ "High Temperature",
                                               (pH_Level == 5.5 & TempLevel == 20) ~ "Low pH",
                                               (pH_Level == 5.5 & TempLevel == 38) ~ "Low pH x High Temperature"))

# Set paramters to correct form for mixed model fitting
d_Sal$Rep <- as.factor(d_Sal$Rep)
d_Sal$Id <- as.factor(d_Sal$Id)
d_Sal$Stress <- as.factor(d_Sal$Stress)
d_Sal$SalinityLevel <- as.numeric(d_Sal$SalinityLevel)

# group the data
df_grp <- groupedData(estimate ~ SalinityLevel | Id, data = d_Sal)

# fit the model
gaus_Sal_nlme <- nlme(estimate ~ gaussian_1987(temp = SalinityLevel, rmax, topt, a),
                      data = df_grp,
                      fixed = list(rmax + topt + a ~ 1 + Stress),
                      random = pdDiag(rmax + topt + a ~ 1),
                      # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                      start = c(0.004, rep(0, 3), 6, rep(0, 3), 20,rep(0, 3)),
                      na.action = na.omit,
                      method = 'REML',
                      control = nlmeControl(pnlsTol = 0.02))
summary(gaus_Sal_nlme)

# fit the model
quad_Sal_nlme <- nlme(estimate ~ quadratic_2008(temp = SalinityLevel, a, b, c),
                      data = df_grp,
                      fixed = list(a + b+ c ~ 1 + Stress),
                      random = pdDiag(a + b+ c ~ 1),
                      # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                      start = c(0.0025, rep(0, 3), 1, rep(0, 3), 1,rep(0, 3)),
                      na.action = na.omit,
                      method = 'REML',
                      control = nlmeControl(pnlsTol = 0.02))

gaus_Sal_nlme_A <- update(gaus_Sal_nlme, start = c(0.004, rep(0, 3), -0.5, rep(0, 3), 20,rep(0, 3)))
gaus_Sal_nlme_B <- update(gaus_Sal_nlme, start = c(0.004, rep(0, 3), -1, rep(0, 3), 20,rep(0, 3)))
gaus_Sal_nlme_C <- update(gaus_Sal_nlme, start = c(0.004, rep(0, 3), -2, rep(0, 3), 20,rep(0, 3)))
gaus_Sal_nlme_D <- update(gaus_Sal_nlme, start = c(0.004, rep(0, 3), 6, rep(0, 3), 20,rep(0, 3)))
gaus_Sal_nlme_E <- update(gaus_Sal_nlme, start = c(0.004, rep(0, 3), 4, rep(0, 3), 20,rep(0, 3)))

anova(gaus_Sal_nlme, gaus_Sal_nlme_A, gaus_Sal_nlme_B, gaus_Sal_nlme_C, gaus_Sal_nlme_D, gaus_Sal_nlme_E)

# check the best model fit and check each combination of random effects (update updates the model above without having to write it all out)
gaus_Sal_nlme_2 <- update(gaus_Sal_nlme, random = pdDiag(rmax + topt ~ 1))
gaus_Sal_nlme_3 <- update(gaus_Sal_nlme, random = pdDiag(rmax + a ~ 1))
gaus_Sal_nlme_4 <- update(gaus_Sal_nlme, random = pdDiag(a + topt ~ 1))
gaus_Sal_nlme_5 <- update(gaus_Sal_nlme, random = pdDiag(rmax ~ 1))
gaus_Sal_nlme_6 <- update(gaus_Sal_nlme, random = pdDiag(topt ~ 1))
gaus_Sal_nlme_7 <- update(gaus_Sal_nlme, random = pdDiag(a ~ 1)) # failed


# Do ANOVA style analysis to get AIC and p values on which model fits best
anova(gaus_Sal_nlme, gaus_Sal_nlme_2, gaus_Sal_nlme_3, gaus_Sal_nlme_4, gaus_Sal_nlme_6, gaus_Sal_nlme_5) # full random effects model is best


# fit the model with ML to look at fixed effects
gaus_Sal_nlme <- nlme(estimate ~ gaussian_1987(temp = SalinityLevel, rmax, topt, a),
                      data = df_grp,
                      fixed = list(rmax + topt + a ~ 1 + Stress),
                      random = pdDiag(rmax + topt + a ~ 1),
                      # start tells the model roughly where to start based on values in the data and the 3 represents the number of variables it will fit the model to not including the control (three other stressor combinations)
                      start = c(0.004, rep(0, 3), 6, rep(0, 3), 20,rep(0, 3)),
                      na.action = na.omit,
                      method = 'ML',
                      control = nlmeControl(pnlsTol = 0.02))

# remove treatment effect on rmax
gaus_Sal_nlme1<- update(gaus_Sal_nlme, fixed = list(rmax ~ 1, topt + a ~ 1 + Stress), start = c(0.004, 6, rep(0, 3), 20, rep(0, 3)))

# remove treatment effect on topt
gaus_Sal_nlme2<- update(gaus_Sal_nlme, fixed = list(topt ~ 1, rmax + a ~ 1 + Stress), start = c(6, 0.004, rep(0, 3), 20, rep(0, 3)))

# remove treatment effect on breadth (a)
gaus_Sal_nlme3<- update(gaus_Sal_nlme, fixed = list(rmax + topt ~ 1 + Stress, a ~ 1), start = c(0.004, rep(0, 3), 6, rep(0, 3), 20)) #failed


anova(gaus_Sal_nlme, gaus_Sal_nlme1, gaus_Sal_nlme2) # full model best

# now remove more and see the effects and which is best model
# remove treatment effect on rmax	 and topt
gaus_Sal_nlme4<- update(gaus_Sal_nlme, fixed = list(rmax + topt ~ 1, a ~ 1 + Stress), start = c(0.004, 6, 20, rep(0, 3)))

# remove treatment effect on topt and a
gaus_Sal_nlme5<- update(gaus_Sal_nlme, fixed = list(topt ~ 1 + Stress, rmax + a  ~ 1), start = c(0.004, 6, rep(0, 3), 20 )) # failed

# remove treatment effect on breadth and topt
gaus_Sal_nlme6<- update(gaus_Sal_nlme, fixed = list(rmax ~ 1 + Stress, topt + a ~ 1 ), start = c(0.004, rep(0, 3), 6, 20))

anova(gaus_Sal_nlme, gaus_Sal_nlme1, gaus_Sal_nlme2, gaus_Sal_nlme4, gaus_Sal_nlme6) # full model best here

# full model best
gaus_Sal_nlme_final <- update(gaus_Sal_nlme, method = 'REML')

#################


#####################
# Create a data frame with the desired Salerature levels
newX <- expand.grid(
  SalinityLevel= seq(0, 35,length = 100),
  Stress= unique(d_Sal$Stress),
  Id= unique(d_Sal$Id))

# Obtain the predicted values for the new Salerature levels
newX$estimate <- predict(quad_Sal_nlme, newdata = newX, level = 0)

Sal_nlme_Plot <- ggplot(d_Sal, aes(SalinityLevel, estimate)) +
  geom_point() +
  geom_line(data=newX, size = 2) +
  facet_grid(~Stress)
Sal_nlme_Plot

summary(quad_Sal_nlme)
fixef(gaus_Sal_nlme)



#################################################
######## Individual OTU Gaussian fits ###########
#################################################

################ load in data ####################
#### Temperature ####
d_temp <-  read.csv("/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/Temp_Aug22/LogisticGrowth_Temp_Aug22/LogisticGrowth_AllOTUs_gompertz_Temp_zeros.csv")

#### Gaussian Temperature #####
Gaussian_Temp_fits <- group_by(d_temp, Id, Stress) %>%
  nest() %>%
  dplyr::mutate(.,Gaussian = purrr::map(data, ~nls_multstart(estimate~gaussian_1987(temp = TempLevel, rmax,topt,a),
                                                             data = .x,
                                                             iter = c(3,3,3),
                                                             start_lower = get_start_vals(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987') - 10,
                                                             start_upper = get_start_vals(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987') + 10,
                                                             lower = get_lower_lims(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             upper = get_upper_lims(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             supp_errors = 'Y',
                                                             convergence_count = FALSE)))

Gaussian_Temp_fits <- group_by(d_temp, Id, Stress) %>%
  nest() %>%
  dplyr::mutate(.,Gaussian = purrr::map(data, ~nls_multstart(estimate~quadratic_2008(temp = TempLevel, a, b, c),
                                                             data = .x,
                                                             iter = c(3,3,3),
                                                             start_lower = get_start_vals(.x$TempLevel, .x$estimate, model_name = 'quadratic_2008') - 10,
                                                             start_upper = get_start_vals(.x$TempLevel, .x$estimate, model_name = 'quadratic_2008') + 10,
                                                             lower = get_lower_lims(.x$TempLevel, .x$estimate, model_name = 'quadratic_2008'),
                                                             upper = get_upper_lims(.x$TempLevel, .x$estimate, model_name = 'quadratic_2008'),
                                                             supp_errors = 'Y',
                                                             convergence_count = FALSE)))

# create new list column of for high resolution data
Temp_d_predsG <- mutate(Gaussian_Temp_fits, new_data = map(data, ~tibble(TempLevel = seq(min(.x$TempLevel), max(.x$TempLevel), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(Gaussian)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # unlist the preds list column
  unnest(preds)

Temp_params_Gaussian <- Gaussian_Temp_fits %>%
  mutate(., params = map(Gaussian, tidy)) %>%
  unnest(params, keep_empty = T) %>%
  dplyr::select(-data, -Gaussian)

Params_temp <- Temp_params_Gaussian %>%
  pivot_wider(id_cols = c(Id, Stress),
              names_from = (term),
              values_from = (estimate))

#write.csv(Temp_params_Gaussian, "/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/TempParamsGaussian.csv")
Temp_d_predsG$Stress = factor(Temp_d_predsG$Stress, levels = c("Control",  "Sal", "pH", "pH_Sal"))
d_temp$Stress = factor(d_temp$Stress, levels = c("Control",  "Sal", "pH", "pH_Sal"))

# plot
All_OTU_Temp_Plot <- ggplot(Temp_d_predsG) +
  geom_line(aes(TempLevel, .fitted, colour = Stress), size = 1) +
  geom_point(aes(TempLevel, estimate, colour = Stress), d_temp) +
  facet_grid(Id ~ Stress, scales = 'free_y') +
  theme_few() +
  scale_colour_manual(values = c( "#666666", "#D95F02", "#1B9E77", "#7570B3"))+
  labs(x = 'Temperature (ºC)',
       y = 'Growth Rate')

All_OTU_Temp_Plot

### All together ####
Gaussian_Temp_fits_all <- group_by(d_temp, Stress) %>%
  nest() %>%
  dplyr::mutate(.,Gaussian = purrr::map(data, ~nls_multstart(estimate~gaussian_1987(temp = TempLevel, rmax,topt,a),
                                                             data = .x,
                                                             iter = c(3,3,3),
                                                             start_lower = get_start_vals(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987') - 10,
                                                             start_upper = get_start_vals(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987') + 10,
                                                             lower = get_lower_lims(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             upper = get_upper_lims(.x$TempLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             supp_errors = 'Y',
                                                             convergence_count = FALSE)))

# create new list column of for high resolution data
Temp_d_predsG <- mutate(Gaussian_Temp_fits_all, new_data = map(data, ~tibble(TempLevel = seq(min(.x$TempLevel), max(.x$TempLevel), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(Gaussian)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # unlist the preds list column
  unnest(preds)

Temp_params_Gaussian <- Gaussian_Temp_fits_all %>%
  mutate(., params = map(Gaussian, tidy)) %>%
  unnest(params, keep_empty = T) %>%
  dplyr::select(-data, -Gaussian)

Params_temp <- Temp_params_Gaussian %>%
  pivot_wider(id_cols = c( Stress),
              names_from = (term),
              values_from = (estimate))

Temp_d_predsG$Stress = factor(Temp_d_predsG$Stress, levels = c("Control",  "Sal", "pH", "pH_Sal"))
d_temp$Stress = factor(d_temp$Stress, levels = c("Control",  "Sal", "pH", "pH_Sal"))

# plot
Overall_plotTemp <- ggplot(d_temp) +
  geom_line(aes(TempLevel, .fitted), Temp_d_predsG, size  = 2) +
  geom_point(aes(TempLevel, estimate, alpha = 2)) +
  scale_alpha(guide = 'none') +
  facet_grid( ~ Stress, scales = 'free_y') +
  theme_few() +
  scale_colour_manual(values = c( "#666666", "#D95F02", "#1B9E77", "#7570B3"))+
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate')

Overall_plotTemp

############ pH ################
################ load in data ####################

d_pH <-  read.csv("/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/pH_Aug22/LogisticGrowth_pH_Aug22/LogisticGrowth_AllOTUs_gompertz_pH_zeros_edited.csv")
d_pH <- d_pH %>% mutate(Stress = case_when((SalinityLevel == 0 & TempLevel == 20) ~ "Control",
                                           (SalinityLevel == 0  & TempLevel == 38) ~ "HighTemp",
                                           (SalinityLevel == 20 & TempLevel == 20) ~ "HighSal",
                                           (SalinityLevel == 20 & TempLevel == 38) ~ "HighSal_HighTemp"))
d_pH <- d_pH %>%
  dplyr::select(-TempLevel, -SalinityLevel, -X, -X.1)

d1_pH <- d_pH %>%
  filter(Id != "D11" | Stress != "HighSal_HighTemp") %>%
  filter(Id != "I22" | Stress != "HighTemp") %>%
  filter(Id != "I23" | Stress != "HighSal_HighTemp") %>%
  filter(Id != "I23" | Stress != "HighTemp")

# fit five chosen model formulation in rTPC
gaussian_pH_fits <- group_by(d1_pH, Id, Stress) %>%
  nest() %>%
  dplyr::mutate(., gaussianMod = purrr::map(data, ~nls_multstart(estimate~gaussian_1987(temp = pH_Level, rmax,topt,a),
                                                                 data = .x,
                                                                 iter = c(3,3,3),
                                                                 start_lower = get_start_vals(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987') - 10,
                                                                 start_upper = get_start_vals(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987') + 10,
                                                                 lower = get_lower_lims(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987'),
                                                                 upper = get_upper_lims(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987'),
                                                                 supp_errors = 'Y',
                                                                 convergence_count = FALSE)))

gaussian_pH_fits <- group_by(d1_pH, Id, Stress) %>%
  nest() %>%
  dplyr::mutate(., gaussianMod = purrr::map(data, ~nls_multstart(estimate~quadratic_2008(temp = pH_Level, a, b, c),
                                                                 data = .x,
                                                                 iter = c(3,3,3),
                                                                 start_lower = get_start_vals(.x$pH_Level, .x$estimate, model_name = 'quadratic_2008') - 10,
                                                                 start_upper = get_start_vals(.x$pH_Level, .x$estimate, model_name = 'quadratic_2008') + 10,
                                                                 lower = get_lower_lims(.x$pH_Level, .x$estimate, model_name = 'quadratic_2008'),
                                                                 upper = get_upper_lims(.x$pH_Level, .x$estimate, model_name = 'quadratic_2008'),
                                                                 supp_errors = 'Y',
                                                                 convergence_count = FALSE)))

# create new list column of for high resolution data
gaussian_preds <- dplyr::mutate(gaussian_pH_fits, new_data = map(data, ~tibble(pH_Level = seq(min(.x$pH_Level), max(.x$pH_Level), length.out = 100)))) %>%
  # get rid of original data column
  dplyr::select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(gaussianMod)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  dplyr::mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  #select(Id, Temp, process, flux, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

pH_params_Gaussian <- gaussian_pH_fits %>%
  mutate(., params = map(gaussianMod, tidy)) %>%
  unnest(params, keep_empty = T) %>%
  dplyr::select(-data, -gaussianMod)

Params_pH <- pH_params_Gaussian %>%
  pivot_wider(id_cols = c(Id, Stress),
              names_from = (term),
              values_from = (estimate))

#write.csv(pH_params_Gaussian, "/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/pHParamsGaussian.csv")

gaussian_preds$Stress = factor(gaussian_preds$Stress, levels = c("Control",  "HighSal", "HighTemp", "HighSal_HighTemp"))
d_pH$Stress = factor(d_pH$Stress, levels = c("Control",  "HighSal", "HighTemp", "HighSal_HighTemp"))

# plot
PanelPlot <- ggplot(gaussian_preds) +
  geom_line(aes(pH_Level, .fitted, col = Stress)) +
  geom_point(aes(pH_Level, estimate, col = Stress), d_pH) +
  facet_grid( Id ~ Stress, scales = 'free_y') +
  theme_few() +
  scale_colour_manual(values = c( "#666666", "#D95F02", "#E7298A", "#7570B3"))+
  labs(x = 'pH Level',
       y = 'Growth rate')
PanelPlot

OverLap_Plot <- ggplot(gaussian_preds) +
  geom_line(aes(pH_Level, .fitted, col = Stress)) +
  geom_point(aes(pH_Level, estimate, col = Stress), d) +
  facet_wrap(~ Id, scales = 'free_y', ncol = 3) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'pH',
       y = 'estimate rate') + ggtitle("gaussian")
OverLap_Plot



#### Get overall fit to all OTUs - pH ####
# fit five chosen model formulation in rTPC
gaussian_fits_all <- group_by(d1_pH, Stress) %>%
  nest() %>%
  dplyr::mutate(., gaussianMod = purrr::map(data, ~nls_multstart(estimate~gaussian_1987(temp = pH_Level, rmax,topt,a),
                                                                 data = .x,
                                                                 iter = c(3,3,3),
                                                                 start_lower = get_start_vals(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987') - 10,
                                                                 start_upper = get_start_vals(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987') + 10,
                                                                 lower = get_lower_lims(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987'),
                                                                 upper = get_upper_lims(.x$pH_Level, .x$estimate, model_name = 'gaussian_1987'),
                                                                 supp_errors = 'Y',
                                                                 convergence_count = FALSE)))









# create new list column of for high resolution data
gaussian_preds <- dplyr::mutate(gaussian_fits_all, new_data = map(data, ~tibble(pH_Level = seq(min(.x$pH_Level), max(.x$pH_Level), length.out = 100)))) %>%
  # get rid of original data column
  dplyr::select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(gaussianMod)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  dplyr::mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  #select(Id, Temp, process, flux, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

gaussian_preds$Stress = factor(gaussian_preds$Stress, levels = c("Control",  "HighSal", "HighTemp", "HighSal_HighTemp"))
d_pH$Stress = factor(d_pH$Stress, levels = c("Control",  "HighSal", "HighTemp", "HighSal_HighTemp"))

# plot
PanelPlot_pH_all<- ggplot(gaussian_preds) +
  geom_line(aes(pH_Level, .fitted), size = 2) +
  geom_point(aes(pH_Level, estimate, alpha = 2), d_pH) +
  scale_alpha(guide = 'none')  +
  facet_grid( ~Stress, scales = 'free_y') +
  theme_few() +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'pH Level',
       y = 'Growth Rate')
PanelPlot_pH_all

####### Salinity #######

Sal_d <- read.csv("/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/Sal_Aug22/LogisticGrowth_Sal_Aug22/LogisticGrowth_AllOTUs_gompertz_Sal_zeros.csv")
Sal_d <- Sal_d %>% mutate(Stress = case_when((pH_Level == 7.2 & TempLevel == 20) ~ "Control",
                                             (pH_Level == 7.2 & TempLevel == 38) ~ "HighTemp",
                                             (pH_Level == 5.5 & TempLevel == 20) ~ "Low_pH",
                                             (pH_Level == 5.5 & TempLevel == 38) ~ "LowpH_HighTemp"))

Sal_d1 <- Sal_d %>% filter(Id != "I18" | Stress != "LowpH_HighTemp") %>%
  dplyr::select(-TempLevel, -pH_Level, -X)

#### Gaussian Salinity #####
Gaussian_Sal_fits <- group_by(Sal_d1, Id, Stress) %>%
  nest() %>%
  dplyr::mutate(.,Gaussian = purrr::map(data, ~nls_multstart(estimate~gaussian_1987(temp = SalinityLevel, rmax,topt,a),
                                                             data = .x,
                                                             iter = c(3,3,3),
                                                             start_lower = get_start_vals(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987') - 10,
                                                             start_upper = get_start_vals(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987') + 10,
                                                             lower = get_lower_lims(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             upper = get_upper_lims(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             supp_errors = 'Y',
                                                             convergence_count = FALSE)))

Gaussian_Sal_fits <- group_by(Sal_d1, Id, Stress) %>%
  nest() %>%
  dplyr::mutate(.,Gaussian = purrr::map(data, ~nls_multstart(estimate~quadratic_2008(temp = SalinityLevel, a, b, c),
                                                             data = .x,
                                                             iter = c(3,3,3),
                                                             start_lower = get_start_vals(.x$SalinityLevel, .x$estimate, model_name = 'quadratic_2008') - 10,
                                                             start_upper = get_start_vals(.x$SalinityLevel, .x$estimate, model_name = 'quadratic_2008') + 10,
                                                             lower = get_lower_lims(.x$SalinityLevel, .x$estimate, model_name = 'quadratic_2008'),
                                                             upper = get_upper_lims(.x$SalinityLevel, .x$estimate, model_name = 'quadratic_2008'),
                                                             supp_errors = 'Y',
                                                             convergence_count = FALSE)))

# create new list column of for high resolution data
Sal_d_predsG <- mutate(Gaussian_Sal_fits, new_data = map(data, ~tibble(SalinityLevel = seq(min(.x$SalinityLevel), max(.x$SalinityLevel), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(Gaussian)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # unlist the preds list column
  unnest(preds)

Sal_params_Gaussian <- Gaussian_Sal_fits %>%
  mutate(., params = map(Gaussian, tidy)) %>%
  unnest(params, keep_empty = T) %>%
  dplyr::select(-data, -Gaussian)

#write.csv(Sal_params_Gaussian, "/Users/HCarmichael/Documents/PhD_ExeterUniversity/Exeter_PhD/Lab_work/Aug22_MultiStress_Individual/MultiStress_Aug22/SalParamsGaussian.csv")


Params_Sal <- Sal_params_Gaussian %>%
  pivot_wider(id_cols = c(Id, Stress),
              names_from = (term),
              values_from = (estimate))

Sal_d_predsG$Stress = factor(Sal_d_predsG$Stress, levels = c("Control", "Low_pH", "HighTemp", "LowpH_HighTemp"))
Sal_d$Stress = factor(Sal_d$Stress, levels = c("Control", "Low_pH", "HighTemp", "LowpH_HighTemp"))


# plot
AllOTU_Sal_plot <- ggplot(Sal_d_predsG) +
  geom_line(aes(SalinityLevel, .fitted, col = Stress), size = 1) +
  geom_point(aes(SalinityLevel, estimate, col = Stress), Sal_d) +
  facet_grid(Id ~ Stress, scales = 'free_y') +
  theme_few() +
  scale_colour_manual(values = c( "#666666", "#1B9E77", "#E7298A", "#7570B3"))+
  labs(x = 'Salinity (g/L)',
       y = 'Growth Rate')

AllOTU_Sal_plot

#### All OTUs ####
Gaussian_Sal_fits_all <- group_by(Sal_d, Stress) %>%
  nest() %>%
  dplyr::mutate(.,Gaussian = purrr::map(data, ~nls_multstart(estimate~gaussian_1987(temp = SalinityLevel, rmax,topt,a),
                                                             data = .x,
                                                             iter = c(3,3,3),
                                                             start_lower = get_start_vals(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987') - 10,
                                                             start_upper = get_start_vals(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987') + 10,
                                                             lower = get_lower_lims(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             upper = get_upper_lims(.x$SalinityLevel, .x$estimate, model_name = 'gaussian_1987'),
                                                             supp_errors = 'Y',
                                                             convergence_count = FALSE)))

# create new list column of for high resolution data
Sal_d_predsG <- mutate(Gaussian_Sal_fits_all, new_data = map(data, ~tibble(SalinityLevel = seq(min(.x$SalinityLevel), max(.x$SalinityLevel), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(Gaussian)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # unlist the preds list column
  unnest(preds)

Sal_params_Gaussian <- Gaussian_Sal_fits_all %>%
  mutate(., params = map(Gaussian, tidy)) %>%
  unnest(params, keep_empty = T) %>%
  dplyr::select(-data, -Gaussian)

Params_Sal <- Sal_params_Gaussian %>%
  pivot_wider(id_cols = c( Stress),
              names_from = (term),
              values_from = (estimate))

Sal_d_predsG$Stress = factor(Sal_d_predsG$Stress, levels = c("Control", "Low_pH", "HighTemp", "LowpH_HighTemp"))
Sal_d$Stress = factor(Sal_d$Stress, levels = c("Control", "Low_pH", "HighTemp", "LowpH_HighTemp"))

# plot
Overall_plotSal <- ggplot(Sal_d_predsG) +
  geom_line(aes(SalinityLevel, .fitted), size = 2) +
  geom_point(aes(SalinityLevel, estimate, alpha = 2), Sal_d) +
  scale_alpha(guide = 'none')  +
  facet_grid( ~ Stress, scales = 'free_y') +
  theme_few() +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Salinity (g/L)',
       y = 'Growth rate')

Overall_plotSal


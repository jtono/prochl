library(ggplot2)
library(dplyr)
library(lme4)
library(tidyverse) #install.packages(tidyverse)
library(zoo) #install.packages(zoo)
library(broom) #install.packages(broom)
library(nls.multstart) # install.packages(nls.multstart)
library(rTPC)
library("growthrates")
library(MuMIn)
########functions########
getslopes <- function(data){
  treatment <- c()
  int <- c()
  sl <- c()
  rep <- c()
  for (i in unique(data$Treatment)){
    sub <- subset(data, Treatment==i)
    mod <- lmer(log(Actual.Cell.count)~day+(day|Rep.ID), sub)
    n <- length(unique(sub$Rep.ID))
    treatment <- c(treatment, rep(i,n))
    int <- c(int, coef(mod)$Rep.ID[,1])
    sl <- c(sl, coef(mod)$Rep.ID[,2])
    rep <- c(rep, rownames(coef(mod)$Rep.ID))
  }
  output <- data.frame(temp=treatment, rep, int, sl)
  return(output)
}

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



#####getting data ready#####
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

#####plot growth curves######
ggplot(data=b_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

ggplot(data=d_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data - all data
b_tpc_sum <- getslopes(b_tpc)
d_tpc_sum <- getslopes(d_tpc)
b_tpc_sumlm <- getslopeslm(b_tpc)
d_tpc_sumlm <- getslopeslm(d_tpc)

#####do rolling regression - 6#####
# create the rolling regression function
roll_regress <- function(x){
  temp <- data.frame(x)
  mod <- lm(temp)
  temp <- data.frame(slope = coef(mod)[[2]],
                     slope_lwr = confint(mod)[2, ][[1]],
                     slope_upr = confint(mod)[2, ][[2]],
                     intercept = coef(mod)[[1]],
                     rsq = summary(mod)$r.squared, stringsAsFactors = FALSE)
  return(temp)
}
# define window - here every 3 points
num_points = 3

# run rolling regression on ln_cell_cnt ~ day
#b
models_b <- b_tpc %>%
  group_by(Rep.ID, Treatment) %>%
  do(cbind(model = dplyr::select(., ln_cell_cnt, day) %>%
             zoo::rollapplyr(width = num_points, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           time = dplyr::select(., day),
           ln_od = dplyr::select(., ln_cell_cnt))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')
# calculate growth rate for each one
b_gr_roll6 <- models_b %>%
  filter(slope == max(slope, na.rm = TRUE)) %>%
  ungroup()

names(b_gr_roll6)[1] <- "rep"
names(b_gr_roll6)[2] <- "temp"
names(b_gr_roll6)[3] <- "sl"

#d
models_d <- d_tpc %>%
  group_by(Rep.ID, Treatment) %>%
  do(cbind(model = dplyr::select(., ln_cell_cnt, day) %>%
             zoo::rollapplyr(width = num_points, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           time = dplyr::select(., day),
           ln_od = dplyr::select(., ln_cell_cnt))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')
# calculate growth rate for each one
d_gr_roll6 <- models_d %>%
  filter(slope == max(slope, na.rm = TRUE)) %>%
  ungroup()

names(d_gr_roll6)[1] <- "rep"
names(d_gr_roll6)[2] <- "temp"
names(d_gr_roll6)[3] <- "sl"

#####fit with growthrates package - 7#####
#a - with easy linear fit
L_b <- all_easylinear(Actual.Cell.count ~ day | Rep.ID + Treatment, data=b_tpc)
b_tpc7a_sum <- results(L_b)

names(b_tpc7a_sum)[1] <- "rep"
names(b_tpc7a_sum)[2] <- "temp"
names(b_tpc7a_sum)[5] <- "sl"

L_d <- all_easylinear(Actual.Cell.count ~ day | Rep.ID + Treatment, data=d_tpc)
d_tpc7a_sum <- results(L_d)

names(d_tpc7a_sum)[1] <- "rep"
names(d_tpc7a_sum)[2] <- "temp"
names(d_tpc7a_sum)[5] <- "sl"

#b - with nonparametric smoothing splines
many_spline_fits_b <- all_splines(Actual.Cell.count ~ day | Rep.ID + Treatment, data = b_tpc)
b_tpc7b_sum <- results(many_spline_fits_b)

names(b_tpc7b_sum)[1] <- "rep"
names(b_tpc7b_sum)[2] <- "temp"
names(b_tpc7b_sum)[4] <- "sl"

many_spline_fits_d <- all_splines(Actual.Cell.count ~ day | Rep.ID + Treatment, data = d_tpc)
d_tpc7b_sum <- results(many_spline_fits_d)

names(d_tpc7b_sum)[1] <- "rep"
names(d_tpc7b_sum)[2] <- "temp"
names(d_tpc7b_sum)[4] <- "sl"

#again with spar set to 0.5 (moderate value)
#https://tpetzoldt.github.io/growthrates/doc/Introduction.html#nonparametric-smoothing-splines
many_spline_fits_b <- all_splines(Actual.Cell.count ~ day | Rep.ID + Treatment, data = b_tpc, spar = 0.5)
b_tpc7b05_sum <- results(many_spline_fits_b)

names(b_tpc7b05_sum)[1] <- "rep"
names(b_tpc7b05_sum)[2] <- "temp"
names(b_tpc7b05_sum)[4] <- "sl"

many_spline_fits_d <- all_splines(Actual.Cell.count ~ day | Rep.ID + Treatment, data = d_tpc, spar = 0.5)
d_tpc7b05_sum <- results(many_spline_fits_d)

names(d_tpc7b05_sum)[1] <- "rep"
names(d_tpc7b05_sum)[2] <- "temp"
names(d_tpc7b05_sum)[4] <- "sl"


#####cut to exponential phase - 1) cut off points clearly decreased at end####
#'get rid of:
#'b - 25 - 16 - all reps
#'b - 27 - 15+16 - all reps
#'d - 20 - 16 - all
#'d - 23 - 15+16 - all
#'d - 25 - 15+16 - all
b_tpc1 <- b_tpc[-which(b_tpc$Treatment==25&b_tpc$day==16|b_tpc$Treatment==27&b_tpc$day%in%c(15,16)),]
d_tpc1 <- d_tpc[-which(d_tpc$Treatment==20&d_tpc$day==16|d_tpc$Treatment==23&d_tpc$day%in%c(15,16)|d_tpc$Treatment==25&d_tpc$day%in%c(15,16)),]

#check growth curves - edited to show fit lines
b_tpc1$temp_rep <- paste(b_tpc1$Treatment,b_tpc1$Rep.ID)
b_tpc1_sumlm$temp_rep <- paste(b_tpc1_sumlm$temp,b_tpc1_sumlm$rep)
p <- ggplot(data=b_tpc1, aes(x=day, y=log(Actual.Cell.count), col=c(Treatment)))+
  geom_point()+
  geom_line(aes(group=temp_rep))+
  facet_wrap(~Treatment)
p+geom_abline(slope=b_tpc1_sumlm$sl, intercept=b_tpc1_sumlm$int, col=c(b_tpc1_sumlm$temp))+facet_wrap(~b_tpc1_sumlm$temp)


ggplot(data=d_tpc1, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data
b_tpc1_sum <- getslopes(b_tpc1)
d_tpc1_sum <- getslopes(d_tpc1)

b_tpc1_sumlm <- getslopeslm(b_tpc1)
d_tpc1_sumlm <- getslopeslm(d_tpc1)

#####cut to exponential phase - 1b) cut off points clearly decreased at end####
#'get rid of last 3 points
b_tpc1b <- b_tpc[b_tpc$day<14,]
d_tpc1b <- d_tpc[d_tpc$day<14,]

#check growth curves
ggplot(data=b_tpc1b, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

ggplot(data=d_tpc1b, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data
b_tpc1b_sumlm <- getslopeslm(b_tpc1b)
d_tpc1b_sumlm <- getslopeslm(d_tpc1b)

#if cut to first 10 is it different?
b_tpc1b10 <- b_tpc[b_tpc$day<10,]
d_tpc1b10 <- d_tpc[d_tpc$day<10,]

b_tpc1b10_sumlm <- getslopeslm(b_tpc1b10)
d_tpc1b10_sumlm <- getslopeslm(d_tpc1b10)

#compare to each other and full
ggplot(data=b_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  geom_vline(xintercept=9.5)+
  geom_vline(xintercept=13.5)+
  facet_wrap(~Treatment)

ggplot(data=d_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  geom_vline(xintercept=9.5)+
  geom_vline(xintercept=13.5)+
  facet_wrap(~Treatment)

b_compare <- merge(b_tpc1b_sumlm, b_tpc1b10_sumlm, by=c("temp","rep"))
b_comp <- merge(b_compare, b_tpc_sumlm, by=c("temp","rep"))
names(b_comp) = c("temp","rep","int1b","sl1b","int1b10","sl1b10","int","sl")
plot(b_comp$sl~jitter(b_comp$temp,2), pch=16, col=1, ylim=c(-0.2,0.8))
points(b_comp$sl1b~jitter(b_comp$temp,2), pch=16, col=2)
points(b_comp$sl1b10~jitter(b_comp$temp,2), pch=16, col=3)

d_compare <- merge(d_tpc1b_sumlm, d_tpc1b10_sumlm, by=c("temp","rep"))
d_comp <- merge(d_compare, d_tpc_sumlm, by=c("temp","rep"))
names(d_comp) = c("temp","rep","int1b","sl1b","int1b10","sl1b10","int","sl")
plot(d_comp$sl~jitter(d_comp$temp,2), pch=16, col=1, ylim=c(-0.2,0.8))
points(d_comp$sl1b~jitter(d_comp$temp,2), pch=16, col=2)
points(d_comp$sl1b10~jitter(d_comp$temp,2), pch=16, col=3)



#####cut to exponential phase - 2) also cut off weird measurements####
#'get rid of:
#'b - 18 - 11 - R3
#'b - 27 - 11 - R3
#'b-29-15-all
#'b-30-5-R1
#'d-15-7-all
#'d-20-7-R1
#'d-30-7-R1

b_tpc2 <- b_tpc1[-which(b_tpc1$Treatment==18&b_tpc1$day==11&b_tpc1$Rep.ID=="R3"|b_tpc1$Treatment==27&b_tpc1$day==11&b_tpc1$Rep.ID=="R3"|b_tpc1$Treatment==29&b_tpc1$day==15|b_tpc1$Treatment==30&b_tpc1$day==5&b_tpc1$Rep.ID=="R1"),]
d_tpc2 <- d_tpc1[-which(d_tpc1$Treatment==15&d_tpc1$day==7|d_tpc1$Treatment==20&d_tpc1$day==7&d_tpc1$Rep.ID=="R1"|d_tpc1$Treatment==30&d_tpc1$day==7&d_tpc1$Rep.ID=="R1"),]

#check growth curves
ggplot(data=b_tpc2, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

ggplot(data=d_tpc2, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data
b_tpc2_sum <- getslopes(b_tpc2)
d_tpc2_sum <- getslopes(d_tpc2)

b_tpc2_sumlm <- getslopeslm(b_tpc2)
d_tpc2_sumlm <- getslopeslm(d_tpc2)

#####cut to exponential phase - 3) trim ends - only days 5-11 a) with other mods####
b_tpc3a <- b_tpc2[-which(b_tpc2$day<5|b_tpc2$day>11),]
d_tpc3a <- d_tpc2[-which(d_tpc2$day<5|d_tpc2$day>11),]

#check growth curves
ggplot(data=b_tpc3a, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

ggplot(data=d_tpc3a, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data
b_tpc3a_sum <- getslopes(b_tpc3a)
d_tpc3a_sum <- getslopes(d_tpc3a)

b_tpc3a_sumlm <- getslopeslm(b_tpc3a)
d_tpc3a_sumlm <- getslopeslm(d_tpc3a)

#####cut to exponential phase - 3) trim ends - only days 5-11 b) without other mods####
b_tpc3b <- b_tpc[-which(b_tpc$day<5|b_tpc$day>11),]
d_tpc3b <- d_tpc[-which(d_tpc$day<5|d_tpc$day>11),]

#check growth curves
ggplot(data=b_tpc3b, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

ggplot(data=d_tpc3b, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data
b_tpc3b_sum <- getslopes(b_tpc3b)
d_tpc3b_sum <- getslopes(d_tpc3b)

b_tpc3b_sumlm <- getslopeslm(b_tpc3b)
d_tpc3b_sumlm <- getslopeslm(d_tpc3b)

#####!!DOESNT WORK>???? cut to exponential phase - 4) try to automatically detect linear portion####
par(mfrow=c(3,3), mar=c(1,1,1,1))
temp <- c()
rep <- c()
sl <- c()
int <- c()
for (i in unique(b_tpc$Treatment)){
  for (j in unique(b_tpc$Rep.ID)){
    sub <- subset(b_tpc, Treatment==i&Rep.ID==j)
    dat <- data.frame(x=sub$day, y=log(sub$Actual.Cell.count))
    mod <- lm(y~x, dat)
    while(summary(mod)$r.squared < 0.95){
      dat <- dat[-order(abs(mod$residuals), decreasing=TRUE)[1],]
      mod <- lm(y~x, dat)
    }
    plot(y~x, data=dat, ylim=c(6,20), xlim=c(0,18))
    temp <- c(temp, i)
    rep <- c(rep, j)
    sl <- c(sl, mod$coefficients[2])
    int <- c(int, mod$coefficients[1])
  }
}

#put data in dataframe
b_tpc4_sum <- data.frame(temp, rep, sl, int)

#d
par(mfrow=c(3,3), mar=c(1,1,1,1))
temp <- c()
rep <- c()
sl <- c()
int <- c()
for (i in unique(d_tpc$Treatment)){
  for (j in unique(d_tpc$Rep.ID)){
    sub <- subset(d_tpc, Treatment==i&Rep.ID==j)
    dat <- data.frame(x=sub$day, y=log(sub$Actual.Cell.count))
    mod <- lm(y~x, dat)
    while(summary(mod)$r.squared < 0.95){
      dat <- dat[-order(abs(mod$residuals), decreasing=TRUE)[1],]
      mod <- lm(y~x, dat)
    }
    plot(y~x, data=dat, ylim=c(6,20), xlim=c(0,18))
    temp <- c(temp, i)
    rep <- c(rep, j)
    sl <- c(sl, mod$coefficients[2])
    int <- c(int, mod$coefficients[1])
  }
}

#put data in dataframe
d_tpc4_sum <- data.frame(temp, rep, sl, int)

#####cut to exponential phase - 5) try max gr a) use code from PhD####
#functions needed
nderiv <- function(fit, x, eps=1e-5)
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)

spline.slope <- function(x, y, n, eps=1e-5)
  max(nderiv(loess(log(y) ~ x, degree=1, span=0.2),
             seq(min(x), max(x), length=n)), na.rm=TRUE)
#d
temp <- c()
rep <- c()
sl <- c()
for (i in unique(d_tpc$Treatment)){
  for (j in unique(d_tpc$Rep.ID)){
    sub <- subset(d_tpc, Treatment==i&Rep.ID==j)
    k <- spline.slope(sub$day, sub$Actual.Cell.count, length(sub$File.ID))
    temp <- c(temp, i)
    rep <- c(rep, j)
    sl <- c(sl, k)
  }
}

#put data in dataframe
d_tpc5a_sum <- data.frame(temp, rep, sl)

#b
temp <- c()
rep <- c()
sl <- c()
for (i in unique(b_tpc$Treatment)){
  for (j in unique(b_tpc$Rep.ID)){
    sub <- subset(b_tpc, Treatment==i&Rep.ID==j)
    k <- spline.slope(sub$day, sub$Actual.Cell.count, length(sub$File.ID))
    temp <- c(temp, i)
    rep <- c(rep, j)
    sl <- c(sl, k)
  }
}

#put data in dataframe
b_tpc5a_sum <- data.frame(temp, rep, sl)

#!!!!many warnings
#!!!didn't work for d at all

#####cut to exponential phase - 5) try max gr b) lm 3 pts with R2 > 0.95####
#d
temp <- c()
rep <- c()
maxsl <- c()
for (i in unique(d_tpc$Treatment)){
  for (j in unique(d_tpc$Rep.ID)){
    sub <- subset(d_tpc, Treatment==i&Rep.ID==j)
    dat <- data.frame(x=sub$day, y=log(sub$Actual.Cell.count))
    sl<-c()
    for (k in 1:(length(dat$x)-2)){
      mod <- lm(y[k:(k+2)]~x[k:(k+2)], dat)
      if(summary(mod)$r.squared > 0.95){
        sl <- c(sl, mod$coefficients[2])}
    }
    maxsl <- c(maxsl,max(sl))
    temp <- c(temp, i)
    rep <- c(rep, j)
    }
}

#put data in dataframe
d_tpc5b_sum <- data.frame(temp, rep, sl=maxsl)

#b
temp <- c()
rep <- c()
maxsl <- c()
for (i in unique(b_tpc$Treatment)){
  for (j in unique(b_tpc$Rep.ID)){
    sub <- subset(b_tpc, Treatment==i&Rep.ID==j)
    dat <- data.frame(x=sub$day, y=log(sub$Actual.Cell.count))
    sl<-c()
    for (k in 1:(length(dat$x)-2)){
      mod <- lm(y[k:(k+2)]~x[k:(k+2)], dat)
      if(summary(mod)$r.squared > 0.95){
        sl <- c(sl, mod$coefficients[2])}
    }
    maxsl <- c(maxsl,max(sl))
    temp <- c(temp, i)
    rep <- c(rep, j)
  }
}

#put data in dataframe
b_tpc5b_sum <- data.frame(temp, rep, sl=maxsl)

#####compare slopes#######
#compare slopes
plot(b_tpc_sum$sl~jitter(b_tpc_sum$temp,2), ylim=c(0,2))
points(b_tpc_sumlm$sl~jitter(b_tpc_sumlm$temp,2), col="grey")
points(b_tpc1_sum$sl~jitter(b_tpc1_sum$temp,2), col="blue")
points(b_tpc1_sumlm$sl~jitter(b_tpc1_sumlm$temp,2), col="darkblue")
points(b_tpc2_sum$sl~jitter(b_tpc2_sum$temp,2), col="red")
points(b_tpc2_sumlm$sl~jitter(b_tpc2_sumlm$temp,2), col="darkred")
points(b_tpc3a_sum$sl~jitter(b_tpc3a_sum$temp,2), col="violet")
points(b_tpc3a_sumlm$sl~jitter(b_tpc3a_sumlm$temp,2), col="purple")
points(b_tpc3b_sum$sl~jitter(b_tpc3b_sum$temp,2), col="green")
points(b_tpc3b_sumlm$sl~jitter(b_tpc3b_sumlm$temp,2), col="darkgreen")
points(b_tpc4_sum$sl~jitter(b_tpc4_sum$temp,2), col="orange", pch=16)
points(b_tpc5a_sum$sl~jitter(b_tpc5a_sum$temp,2), col="lightblue", pch=16)
points(b_tpc5b_sum$sl~jitter(b_tpc5b_sum$temp,2), col="yellow", pch=16)
points(b_gr_roll6$sl~jitter(b_gr_roll6$temp,2), col="pink", pch=18)
points(b_tpc7a_sum$sl~jitter(b_tpc7a_sum$temp,2), col="cyan", pch=18)
points(b_tpc7b_sum$sl~jitter(b_tpc7b_sum$temp,2), col="magenta", pch=18)
points(b_tpc7b05_sum$sl~jitter(b_tpc7b05_sum$temp,2), col="brown", pch=18)
#'similar by eye except light blue
#'some points for pink, magenta quite different (because didn't trim data?)

plot(d_tpc_sum$sl~jitter(d_tpc_sum$temp,2), ylim=c(0,1))
points(d_tpc_sumlm$sl~jitter(d_tpc_sumlm$temp,2), col="grey")
points(d_tpc1_sum$sl~jitter(d_tpc1_sum$temp,2), col="blue")
points(d_tpc1_sumlm$sl~jitter(d_tpc1_sumlm$temp,2), col="darkblue")
points(d_tpc2_sum$sl~jitter(d_tpc2_sum$temp,2), col="red")
points(d_tpc2_sumlm$sl~jitter(d_tpc2_sumlm$temp,2), col="darkred")
points(d_tpc3a_sum$sl~jitter(d_tpc3a_sum$temp,2), col="violet")
points(d_tpc3a_sumlm$sl~jitter(d_tpc3a_sumlm$temp,2), col="purple")
points(d_tpc3b_sum$sl~jitter(d_tpc3b_sum$temp,2), col="green")
points(d_tpc3b_sumlm$sl~jitter(d_tpc3b_sumlm$temp,2), col="darkgreen")
points(d_tpc4_sum$sl~jitter(d_tpc4_sum$temp,2), col="orange", pch=16)
points(d_tpc5a_sum$sl~jitter(d_tpc5a_sum$temp,2), col="lightblue", pch=16)
points(d_tpc5b_sum$sl~jitter(d_tpc5b_sum$temp,2), col="yellow", pch=16)
points(d_gr_roll6$sl~jitter(d_gr_roll6$temp,2), col="pink", pch=18)
points(d_tpc7a_sum$sl~jitter(d_tpc7a_sum$temp,2), col="cyan", pch=18)
points(d_tpc7b_sum$sl~jitter(d_tpc7b_sum$temp,2), col="magenta", pch=18)
points(d_tpc7b05_sum$sl~jitter(d_tpc7b05_sum$temp,2), col="brown", pch=18)
#pretty similar by eye except light blue and yellow
#'light blue very wrong
#'yellow, pink, magenta (although magenta even higher)  much higher for higher ones
#'black and grey lower for higher growing ones

#######fit tpcs - b#####

#'all data sources for b:
data_b <- list("lmer" = b_tpc_sum,"lm" = b_tpc_sumlm, "lmer1"=b_tpc1_sum, "lm1"=b_tpc1_sumlm, "lmer2"=b_tpc2_sum,"lm2"=b_tpc2_sumlm, "lmer3a"=b_tpc3a_sum,"lm3a"=b_tpc3a_sumlm,"lmer3b"=b_tpc3b_sum,"lm3b"=b_tpc3b_sumlm,"gr4"=b_tpc4_sum,"gr5a"=b_tpc5a_sum, "gr5b"=b_tpc5b_sum, "gr6"=b_gr_roll6, "gr7a"=b_tpc7a_sum, "gr7b"=b_tpc7b_sum, "gr7b05"=b_tpc7b05_sum)

#'for each data source, fit all models in: -	Boatman_2017, sharpeschoolfull_1981, modifiedgaussian_2006, oneill_1972, Thomas_2012, briere2_1999, quadratic_2008, johnsonlewin_1946 - cut cuz errors, Hinshelwood_1947, lactin2_1995 added cuz good with -ve values
#'extract convergence tolerance and AIC for each
#'keep track of: dataset, rep, model, AICc
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()

#Boatman_2017
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "boatman_2017")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "boatman_2017")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "boatman_2017")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "boatman_2017")
    fit <- nls_multstart(sl~boatman_2017(temp = temp, rmax, tmin, tmax, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

b_fit_results <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

#sharpeschoolfull_1981, tref=23
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "sharpeschoolfull_1981")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "sharpeschoolfull_1981")
    start_vals[which(is.na(start_vals))]<-1
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "sharpeschoolfull_1981")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "sharpeschoolfull_1981")
    fit <- nls_multstart(sl~sharpeschoolfull_1981(temp = temp, r_tref, e, e1, t1, eh, th, tref=23),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res2 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res2)

#modifiedgaussian_2006
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "modifiedgaussian_2006")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    fit <- nls_multstart(sl~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res3 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res3)

#oneill_1972
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "oneill_1972")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "oneill_1972")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    fit <- nls_multstart(sl~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res4 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res4)

#Thomas_2012
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "thomas_2012")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "thomas_2012")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    fit <- nls_multstart(sl~thomas_2012(temp = temp, a,b,c,tref),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 2,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res5 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res5)

#briere2_1999
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "briere2_1999")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "briere2_1999")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "briere2_1999")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "briere2_1999")
    fit <- nls_multstart(sl~briere2_1999(temp = temp, tmin, tmax, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res6 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res6)

#quadratic_2008
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "quadratic_2008")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "quadratic_2008")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "quadratic_2008")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "quadratic_2008")
    fit <- nls_multstart(sl~quadratic_2008(temp = temp, a, b, c),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 0.5,
                         start_upper = start_vals + 0.5,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res7 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res7)

#lactin2_1995
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "lactin2_1995")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "lactin2_1995")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "lactin2_1995")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "lactin2_1995")
    fit <- nls_multstart(sl~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 1,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res8 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res8)

#Hinshelwood_1947
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b)){
  for (j in unique(data_b[[i]]$rep)){
    sub <- subset(data_b[[i]], rep==j)
    d_name <- c(d_name,names(data_b[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "hinshelwood_1947")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "hinshelwood_1947")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "hinshelwood_1947")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "hinshelwood_1947")
    fit <- nls_multstart(sl~hinshelwood_1947(temp = temp, a, e, b, eh),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 1,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res9 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res9)

write.csv(b_fit_results, "b_fits.csv")
#'(before added 6, 7) best model according to AICc is quadratic_2008, best according to AIC is modifiedgaussian_2006 (done by hand)
b_fits <- read.csv("b_fits.csv", header=TRUE)

mean(b_fits[b_fits$mod=="quadratic_2008",]$topt)
#23.58784
mean(b_fits[b_fits$mod=="modifiedgaussian_2006",]$topt)
#23.87157

b_fits[which(b_fits$aicc==min(b_fits$aicc)),]
#  X d_name rep_name            mod       aic    aicc  topt
#349 349   gr7a       R1 quadratic_2008 -21.29823 -7.9649 24.57
b_fits[which(b_fits$aic==min(b_fits$aic)),]
# X d_name rep_name          mod       aic     aicc  topt
#18 18    lm2       R3 boatman_2017 -38.80806 45.19194 24.67

#find overall best fitting gr data - by median - using aicc
b_fits_gr <- aggregate(b_fits$aicc, list(b_fits$d_name), FUN=median)
b_fits_gr[which(b_fits_gr$x==min(b_fits_gr$x)),]
#lm1 19.55381
b_fits_mod <- aggregate(b_fits$aicc, list(b_fits$mod), FUN=median)
b_fits_mod[which(b_fits_mod$x==min(b_fits_mod$x)),]
#quadratic_2008 5.244838



#######fit tpcs - d#####
#'all data sources for d:
#'first get rid of -Inf
d_tpc5b_sum <- d_tpc5b_sum[-which(d_tpc5b_sum$sl==-Inf),]
#'now put together
data_d <- list("lmer" = d_tpc_sum,"lm" = d_tpc_sumlm, "lmer1"=d_tpc1_sum, "lm1"=d_tpc1_sumlm, "lmer2"=d_tpc2_sum,"lm2"=d_tpc2_sumlm, "lmer3a"=d_tpc3a_sum,"lm3a"=d_tpc3a_sumlm,"lmer3b"=d_tpc3b_sum,"lm3b"=d_tpc3b_sumlm,"gr4"=d_tpc4_sum,"gr5a"=d_tpc5a_sum, "gr5b"=d_tpc5b_sum,"gr6"=d_gr_roll6, "gr7a"=d_tpc7a_sum, "gr7b"=d_tpc7b_sum, "gr7b05"=d_tpc7b05_sum)


#'for each data source, fit all models in: -	Boatman_2017, sharpeschoolfull_1981, modifiedgaussian_2006, oneill_1972, Thomas_2012, briere2_1999, quadratic_2008, johnsonlewin_1946 - cut cuz errors, Hinshelwood_1947, lactin2_1995 added cuz good with -ve values
#'extract convergence tolerance and AIC for each
#'keep track of: dataset, rep, model, AICc
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()

#Boatman_2017
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "boatman_2017")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "boatman_2017")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "boatman_2017")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "boatman_2017")
    fit <- nls_multstart(sl~boatman_2017(temp = temp, rmax, tmin, tmax, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

d_fit_results <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

#sharpeschoolfull_1981, tref=23
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "sharpeschoolfull_1981")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "sharpeschoolfull_1981")
    start_vals[which(is.na(start_vals))]<-1
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "sharpeschoolfull_1981")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "sharpeschoolfull_1981")
    fit <- nls_multstart(sl~sharpeschoolfull_1981(temp = temp, r_tref, e, e1, t1, eh, th, tref=23),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res2 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res2)

#modifiedgaussian_2006
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "modifiedgaussian_2006")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    fit <- nls_multstart(sl~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res3 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res3)

#oneill_1972
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "oneill_1972")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "oneill_1972")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    fit <- nls_multstart(sl~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res4 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res4)

#Thomas_2012
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "thomas_2012")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "thomas_2012")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    fit <- nls_multstart(sl~thomas_2012(temp = temp, a,b,c,tref),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 2,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res5 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res5)

#briere2_1999
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "briere2_1999")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "briere2_1999")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "briere2_1999")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "briere2_1999")
    fit <- nls_multstart(sl~briere2_1999(temp = temp, tmin, tmax, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res6 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res6)

#quadratic_2008
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "quadratic_2008")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "quadratic_2008")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "quadratic_2008")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "quadratic_2008")
    fit <- nls_multstart(sl~quadratic_2008(temp = temp, a, b, c),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 0.5,
                         start_upper = start_vals + 0.5,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res7 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res7)

#lactin2_1995
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "lactin2_1995")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "lactin2_1995")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "lactin2_1995")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "lactin2_1995")
    fit <- nls_multstart(sl~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 1,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res8 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res8)

#Hinshelwood_1947
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d)){
  for (j in unique(data_d[[i]]$rep)){
    sub <- subset(data_d[[i]], rep==j)
    d_name <- c(d_name,names(data_d[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "hinshelwood_1947")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "hinshelwood_1947")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "hinshelwood_1947")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "hinshelwood_1947")
    fit <- nls_multstart(sl~hinshelwood_1947(temp = temp, a, e, b, eh),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 1,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res9 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res9)

write.csv(d_fit_results, "d_fits.csv")


#'(before added 6, 7) best model according to AICc is quadratic_2008, best according to AIC is sharpeschoolfull_1981 (which doesn't work with AICc), next best with AIC was gaussian (done by hand)
d_fits <- read.csv("d_fits.csv", header=TRUE)

mean(d_fits[d_fits$mod=="quadratic_2008",]$topt)
#21.3375
mean(d_fits[d_fits$mod=="modifiedgaussian_2006",]$topt)
#21.23562
mean(d_fits[d_fits$mod=="sharpeschoolfull_1981",]$topt)
#23.46187

d_fits[which(d_fits$aicc==min(d_fits$aicc)),]
#X d_name rep_name         mod       aic      aicc topt
#226 226   gr5a       R1 thomas_2012 -121.9215 -91.92155   30
d_fits[which(d_fits$aic==min(d_fits$aic)),]
#   X d_name rep_name         mod       aic      aicc topt
#226 226   gr5a       R1 thomas_2012 -121.9215 -91.92155   30

#find overall best fitting gr data - by median - using aicc
d_fits_gr <- aggregate(b_fits$aicc, list(b_fits$d_name), FUN=median)
d_fits_gr[which(b_fits_gr$x==min(b_fits_gr$x)),]
#Group.1        x
#9     lm1 19.55381
d_fits_gr[which(b_fits_gr$x<20),]


d_fits_mod <- aggregate(d_fits$aicc, list(d_fits$mod), FUN=median)
d_fits_mod[which(d_fits_mod$x==min(d_fits_mod$x)),]
#        Group.1          x
#7 quadratic_2008 -0.5454946
d_fits_mod[which(d_fits_mod$x<10),]

######plot data with gr's####
#overall lm1 was best, and quadratic fit
#plot these to data
b_gr<- b_tpc1_sumlm
names(b_gr) <- c("Treatment","Rep.ID","int","sl")

#plotted to original data
p <- ggplot(data=b_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)
p
p + geom_abline(data=b_gr, aes(slope=sl, intercept=int, col=Rep.ID))

d_gr<- d_tpc1_sumlm
names(d_gr) <- c("Treatment","Rep.ID","int","sl")

#plotted to original data
p <- ggplot(data=d_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)
p
p + geom_abline(data=d_gr, aes(slope=sl, intercept=int, col=Rep.ID))


######plot tpcs#####
library(lattice)
#use quadratic_2008 fits for lm1
xyplot(sl ~ temp|rep, data=b_tpc1_sumlm)
abline(h=25)
xyplot(sl ~ temp|rep, data=d_tpc1_sumlm)

ggplot(b_tpc1_sumlm, aes(temp, sl)) +
  geom_point(aes(temp, sl), b_tpc1_sumlm) +
  #geom_line(aes(temp, .fitted), col = 'blue') +
  #facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth Rate',
       title='Prochlorococcus B') +
  geom_hline(aes(yintercept = 0), linetype = 2)

ggplot(d_tpc1_sumlm, aes(temp, sl)) +
  geom_point(aes(temp, sl), d_tpc1_sumlm) +
  #geom_line(aes(temp, .fitted), col = 'blue') +
  #facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth Rate',
       title='Prochlorococcus D') +
  geom_hline(aes(yintercept = 0), linetype = 2)

#######fit tpcs - b - only those with topt and using only lm - full, cut@13, cut@9#####

#'all data sources for b:
data_b_lm <- list("lm" = b_tpc_sumlm, "lm1b"=b_tpc1b_sumlm, "lm1b10"=b_tpc1b10_sumlm)

#'for each data source, fit all models in: deutsch_2008, gaussian_1987, joehnk_2008, johnsonlewin_1946, lrf_1991, modifiedgaussian_2006, oneill_1972, pawar_2018, thomas_2012, weibull_1995
#'keep track of: dataset, rep, model, AIC, AICc, topt
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()

#modifiedgaussian_2006
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "modifiedgaussian_2006")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    fit <- nls_multstart(sl~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

b_fit_results <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

#oneill_1972
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "oneill_1972")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "oneill_1972")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    fit <- nls_multstart(sl~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res2 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res2)

#Thomas_2012
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "thomas_2012")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "thomas_2012")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    fit <- nls_multstart(sl~thomas_2012(temp = temp, a,b,c,tref),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 2,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res3 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res3)

# #deutsch_2008
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "deutsch_2008")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "deutsch_2008")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "deutsch_2008")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "deutsch_2008")
    fit <- nls_multstart(sl~deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res4 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res4)

#gaussian_1987
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "gaussian_1987")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "gaussian_1987")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "gaussian_1987")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "gaussian_1987")
    fit <- nls_multstart(sl~gaussian_1987(temp = temp,rmax, topt, a),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}


res5 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res5)

#joehnk_2008
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "joehnk_2008")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "joehnk_2008")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "joehnk_2008")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "joehnk_2008")
    fit <- nls_multstart(sl~joehnk_2008(temp = temp,rmax, topt, a, b, c),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}
res6 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res6)

#johnsonlewin_1946
# d_name <- c()
# rep_name <- c()
# mod <- c()
# aic <- c()
# aicc <- c()
# topt <- c()
# for (i in 1:length(data_b_lm)){
#   for (j in unique(data_b_lm[[i]]$rep)){
#     sub <- subset(data_b_lm[[i]], rep==j)
#     d_name <- c(d_name,names(data_b_lm[i]))
#     rep_name <- c(rep_name, j)
#     mod <- c(mod, "johnsonlewin_1946")
#
#     # get start vals
#     start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "johnsonlewin_1946")
#     # get limits
#     low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "johnsonlewin_1946")
#     upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "johnsonlewin_1946")
#     fit <- nls_multstart(sl~johnsonlewin_1946(temp = temp,r0, e, eh, topt),
#                          data = sub,
#                          iter = 500,
#                          start_lower = start_vals - 1,
#                          start_upper = start_vals + 1,
#                          lower = low_lims,
#                          upper = upper_lims,
#                          supp_errors = 'Y',
#                          convergence_count=FALSE)
#     aic <- c(aic, AIC(fit))
#     aicc <- c(aicc, AICc(fit))
#     param <- calc_params(fit) %>%mutate_all(round, 2)
#     topt <- c(topt, param$topt)
#   }
# }
#
# res7 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)
#
# b_fit_results <- rbind(b_fit_results, res7)

#lrf_1991
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "lrf_1991")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "lrf_1991")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "lrf_1991")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "lrf_1991")
    fit <- nls_multstart(sl~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res8 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res8)


#pawar_2018
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "pawar_2018")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "pawar_2018")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "pawar_2018")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "pawar_2018")
    fit <- nls_multstart(sl~pawar_2018(temp = temp, r_tref, e, eh, topt, tref=23),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res9 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res9)

#weibull_1995
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_b_lm)){
  for (j in unique(data_b_lm[[i]]$rep)){
    sub <- subset(data_b_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_b_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "weibull_1995")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "weibull_1995")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "weibull_1995")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "weibull_1995")
    fit <- nls_multstart(sl~weibull_1995(temp = temp, a, topt, b, c),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res10 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

b_fit_results <- rbind(b_fit_results, res10)

#write model fit output to file
write.csv(b_fit_results, "b_fits_lm.csv")
#load model fit outputs
b_fits <- read.csv("b_fits_lm.csv", header=TRUE)



#######fit tpcs - d - only those with topt and using only lm - full, cut@13, cut@9#####

#'all data sources for d:
data_d_lm <- list("lm" = d_tpc_sumlm, "lm1b"=d_tpc1b_sumlm, "lm1b10"=d_tpc1b10_sumlm)

#'for each data source, fit all models in: deutsch_2008, gaussian_1987, joehnk_2008, johnsonlewin_1946, lrf_1991, modifiedgaussian_2006, oneill_1972, pawar_2018, thomas_2012, weibull_1995
#'keep track of: dataset, rep, model, AIC, AICc, topt
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()

#modifiedgaussian_2006
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "modifiedgaussian_2006")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "modifiedgaussian_2006")
    fit <- nls_multstart(sl~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

d_fit_results <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

#oneill_1972
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "oneill_1972")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "oneill_1972")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "oneill_1972")
    fit <- nls_multstart(sl~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res2 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res2)

#Thomas_2012
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "thomas_2012")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "thomas_2012")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "thomas_2012")
    fit <- nls_multstart(sl~thomas_2012(temp = temp, a,b,c,tref),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 2,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res3 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res3)

# #deutsch_2008
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "deutsch_2008")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "deutsch_2008")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "deutsch_2008")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "deutsch_2008")
    fit <- nls_multstart(sl~deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res4 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res4)

#gaussian_1987
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "gaussian_1987")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "gaussian_1987")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "gaussian_1987")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "gaussian_1987")
    fit <- nls_multstart(sl~gaussian_1987(temp = temp,rmax, topt, a),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}


res5 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res5)

#joehnk_2008
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "joehnk_2008")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "joehnk_2008")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "joehnk_2008")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "joehnk_2008")
    fit <- nls_multstart(sl~joehnk_2008(temp = temp,rmax, topt, a, b, c),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}
res6 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res6)

#johnsonlewin_1946
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "johnsonlewin_1946")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "johnsonlewin_1946")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "johnsonlewin_1946")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "johnsonlewin_1946")
    fit <- nls_multstart(sl~johnsonlewin_1946(temp = temp,r0, e, eh, topt),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 1,
                         start_upper = start_vals + 1,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res7 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res7)

#lrf_1991
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "lrf_1991")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "lrf_1991")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "lrf_1991")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "lrf_1991")
    fit <- nls_multstart(sl~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res8 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res8)


#pawar_2018
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "pawar_2018")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "pawar_2018")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "pawar_2018")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "pawar_2018")
    fit <- nls_multstart(sl~pawar_2018(temp = temp, r_tref, e, eh, topt, tref=23),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res9 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res9)

#weibull_1995
d_name <- c()
rep_name <- c()
mod <- c()
aic <- c()
aicc <- c()
topt <- c()
for (i in 1:length(data_d_lm)){
  for (j in unique(data_d_lm[[i]]$rep)){
    sub <- subset(data_d_lm[[i]], rep==j)
    d_name <- c(d_name,names(data_d_lm[i]))
    rep_name <- c(rep_name, j)
    mod <- c(mod, "weibull_1995")

    # get start vals
    start_vals <- get_start_vals(sub$temp, sub$sl, model_name = "weibull_1995")
    # get limits
    low_lims <- get_lower_lims(sub$temp, sub$sl, model_name = "weibull_1995")
    upper_lims <- get_upper_lims(sub$temp, sub$sl, model_name = "weibull_1995")
    fit <- nls_multstart(sl~weibull_1995(temp = temp, a, topt, b, c),
                         data = sub,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y',
                         convergence_count=FALSE)
    aic <- c(aic, AIC(fit))
    aicc <- c(aicc, AICc(fit))
    param <- calc_params(fit) %>%mutate_all(round, 2)
    topt <- c(topt, param$topt)
  }
}

res10 <- data.frame(d_name, rep_name, mod, aic, aicc, topt)

d_fit_results <- rbind(d_fit_results, res10)
#write output to a file
write.csv(d_fit_results, "d_fits_lm.csv")


####find topts for models of interst####

#load model fits output
b_fits <- read.csv("b_fits_lm.csv", header=TRUE)
d_fits <- read.csv("d_fits_lm.csv", header=TRUE)

#gaussian_1987
gau_b <- subset(b_fits, mod=="gaussian_1987"&d_name=="lm1b")
mean(gau_b$topt)
#23.27333
gau_d <- subset(d_fits, mod=="gaussian_1987"&d_name=="lm1b")
mean(gau_d$topt)
#21.08

#modifiedgaussian_2006
gau_b <- subset(b_fits, mod=="modifiedgaussian_2006")
mean(gau_b$topt)
#23.47222
gau_d <- subset(d_fits, mod=="modifiedgaussian_2006")
mean(gau_d$topt)
#20.68111

#lrf_1991
gau_b <- subset(b_fits, mod=="lrf_1991")
mean(gau_b$topt)
#28.81222
gau_d <- subset(d_fits, mod=="lrf_1991")
mean(gau_d$topt)
#21.11111

#pawar_2018
gau_b <- subset(b_fits, mod=="pawar_2018")
mean(gau_b$topt)
#25.39889
gau_d <- subset(d_fits, mod=="pawar_2018")
mean(gau_d$topt)
#23.83333



#####plot lm tpcs#####
###right now this is average curve
data_b_lm <- list("lm" = b_tpc_sumlm, "lm1b"=b_tpc1b_sumlm, "lm1b10"=b_tpc1b10_sumlm)
#all data lm
d <- subset(d_tpc1b_sumlm, rep=="R1")
# show the data
ggplot(d, aes(temp, sl)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'ProB lm13')

# choose model
mod = 'gaussian_1987'

# get start vals
start_vals <- get_start_vals(d$temp, d$sl, model_name = 'gaussian_1987')

# get limits
low_lims <- get_lower_lims(d$temp, d$sl, model_name = 'gaussian_1987')
upper_lims <- get_upper_lims(d$temp, d$sl, model_name = 'gaussian_1987')

start_vals
#rmax       topt          a
#0.5860897 25.0000000 15.0000000
low_lims
#rmax        topt           a
#-0.01809075 15.00000000  0.00000000

upper_lims
#rmax       topt          a
#5.860897  30.000000 150.000000

# fit model
fit <- nls_multstart(sl~gaussian_1987(temp = temp, rmax, topt, a),
                     data = d,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y',
                     convergence_count=FALSE)

fit
# Nonlinear regression model
# model: sl ~ gaussian_1987(temp = temp, rmax, topt, a)
# data: data
# rmax    topt       a
# 0.6068 22.9852  4.2478
# residual sum-of-squares: 0.2067
#
# Number of iterations to convergence: 19
# Achieved convergence tolerance: 1.49e-08

# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)
#rmax  topt ctmin ctmax    e   eh   q10 thermal_safety_margin thermal_tolerance breadth skewness
# 0.61 22.98 12.59 33.38 2.44 1.11 28.02                  10.4              20.8    5.67     1.33

# predict new data
new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)

# plot data and model fit
ggplot(d, aes(temp, sl)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'ProB lm')


####many curves####
b_tpc_sumlm$method <- "b.lm"
b_tpc1b_sumlm$method <- "b.lm13"
b_tpc1b10_sumlm$method <- "b.lm9"
d_tpc_sumlm$method <- "d.lm"
d_tpc1b_sumlm$method <- "d.lm13"
d_tpc1b10_sumlm$method <- "d.lm9"

all_tpc_data <- rbind(b_tpc_sumlm, b_tpc1b_sumlm, b_tpc1b10_sumlm, d_tpc_sumlm, d_tpc1b_sumlm, d_tpc1b10_sumlm)

d <- all_tpc_data

# when scaling up our code to fit hundreds of models, its nice to have a progress bar
# edit nls_multstart to allow for a progress bar
nls_multstart_progress <- function(formula, data = parent.frame(), iter, start_lower,
                                   start_upper, supp_errors = c("Y", "N"), convergence_count = 100,
                                   control, modelweights, ...){
  if(!is.null(pb)){
    pb$tick()
  }
  nls_multstart(formula = formula, data = data, iter = iter, start_lower = start_lower,
                start_upper = start_upper, supp_errors = supp_errors, convergence_count = convergence_count,
                control = control, modelweights = modelweights, ...)
}

# start progress bar and estimate time it will take
number_of_models <- 2
number_of_curves <- length(unique(d$rep))*length(unique(d$method))

# setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models,
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")

# fit two chosen model formulation in rTPC
d_fits <- nest(d, data = c(temp, sl, int)) %>%
  mutate(gaussian = map(data, ~nls_multstart(sl~gaussian_1987(temp = temp, rmax,topt,a),
                                                      data = .x,
                                                      iter = c(3,3,3),
                                                      start_lower = get_start_vals(.x$temp, .x$sl, model_name = 'gaussian_1987') - 10,
                                                      start_upper = get_start_vals(.x$temp, .x$sl, model_name = 'gaussian_1987') + 10,
                                                      lower = get_lower_lims(.x$temp, .x$sl, model_name = 'gaussian_1987'),
                                                      upper = get_upper_lims(.x$temp, .x$sl, model_name = 'gaussian_1987'),
                                                      supp_errors = 'Y',
                                                      convergence_count = FALSE)),
         lrf = map(data, ~nls_multstart(sl~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                                 data = .x,
                                                 iter = 500,
                                                 start_lower = get_start_vals(.x$temp, .x$sl, model_name = 'lrf_1991') - 10,
                                                 start_upper = get_start_vals(.x$temp, .x$sl, model_name = 'lrf_1991') + 10,
                                                 lower = get_lower_lims(.x$temp, .x$sl, model_name = 'lrf_1991'),
                                                 upper = get_upper_lims(.x$temp, .x$sl, model_name = 'lrf_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)))


d_preds <- mutate(d_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(gaussian, lrf)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(rep, method, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(d_preds)

# plot
ggplot(d_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, sl), d) +
  facet_wrap(~rep*method, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'All fitted thermal performance curves')






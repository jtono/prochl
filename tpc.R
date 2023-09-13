library(ggplot2)
library(dplyr)
library(lme4)
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
  output <- data.frame(treatment, rep, int, sl)
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
  output <- data.frame(treatment, rep, int, sl)
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

#####cut to exponential phase - 1) cut off points clearly decreased at end####
#'get rid of:
#'b - 25 - 16 - all reps
#'b - 27 - 15+16 - all reps
#'d - 20 - 16 - all
#'d - 23 - 15+16 - all
#'d - 25 - 15+16 - all
b_tpc1 <- b_tpc[-which(b_tpc$Treatment==25&b_tpc$day==16|b_tpc$Treatment==27&b_tpc$day%in%c(15,16)),]
d_tpc1 <- d_tpc[-which(d_tpc$Treatment==20&d_tpc$day==16|d_tpc$Treatment==23&d_tpc$day%in%c(15,16)|d_tpc$Treatment==25&d_tpc$day%in%c(15,16)),]

#check growth curves
ggplot(data=b_tpc1, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

ggplot(data=d_tpc1, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)

#do linear regression to get data
b_tpc1_sum <- getslopes(b_tpc1)
d_tpc1_sum <- getslopes(d_tpc1)

b_tpc1_sumlm <- getslopeslm(b_tpc1)
d_tpc1_sumlm <- getslopeslm(d_tpc1)

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

#####!!DOESNT WORK!!! cut to exponential phase - 4) try to automatically detect linear portion####
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
d_tpc5b_sum <- data.frame(temp, rep, maxsl)

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
b_tpc5b_sum <- data.frame(temp, rep, maxsl)

#####compare slopes#######
#compare slopes
plot(b_tpc_sum$sl~jitter(b_tpc_sum$treatment,2), ylim=c(0,2))
points(b_tpc1_sumlm$sl~jitter(b_tpc1_sumlm$treatment,2), col="grey")
points(b_tpc1_sum$sl~jitter(b_tpc1_sum$treatment,2), col="blue")
points(b_tpc1_sumlm$sl~jitter(b_tpc1_sumlm$treatment,2), col="darkblue")
points(b_tpc2_sum$sl~jitter(b_tpc2_sum$treatment,2), col="red")
points(b_tpc2_sumlm$sl~jitter(b_tpc2_sumlm$treatment,2), col="darkred")
points(b_tpc3a_sum$sl~jitter(b_tpc3a_sum$treatment,2), col="violet")
points(b_tpc3a_sumlm$sl~jitter(b_tpc3a_sumlm$treatment,2), col="purple")
points(b_tpc3b_sum$sl~jitter(b_tpc3b_sum$treatment,2), col="green")
points(b_tpc3b_sumlm$sl~jitter(b_tpc3b_sumlm$treatment,2), col="darkgreen")
points(b_tpc4_sum$sl~jitter(b_tpc4_sum$temp,2), col="orange", pch=16)
points(b_tpc5a_sum$sl~jitter(b_tpc5a_sum$temp,2), col="lightblue", pch=16)
points(b_tpc5b_sum$maxsl~jitter(b_tpc5b_sum$temp,2), col="yellow", pch=16)
#'similar by eye except light blue

plot(d_tpc_sum$sl~jitter(d_tpc_sum$treatment,2), ylim=c(0,1))
points(d_tpc1_sumlm$sl~jitter(d_tpc1_sumlm$treatment,2), col="grey")
points(d_tpc1_sum$sl~jitter(d_tpc1_sum$treatment,2), col="blue")
points(d_tpc1_sumlm$sl~jitter(d_tpc1_sumlm$treatment,2), col="darkblue")
points(d_tpc2_sum$sl~jitter(d_tpc2_sum$treatment,2), col="red")
points(d_tpc2_sumlm$sl~jitter(d_tpc2_sumlm$treatment,2), col="darkred")
points(d_tpc3a_sum$sl~jitter(d_tpc3a_sum$treatment,2), col="violet")
points(d_tpc3a_sumlm$sl~jitter(d_tpc3a_sumlm$treatment,2), col="purple")
points(d_tpc3b_sum$sl~jitter(d_tpc3b_sum$treatment,2), col="green")
points(d_tpc3b_sumlm$sl~jitter(d_tpc3b_sumlm$treatment,2), col="darkgreen")
points(d_tpc4_sum$sl~jitter(d_tpc4_sum$temp,2), col="orange", pch=16)
points(d_tpc5a_sum$sl~jitter(d_tpc5a_sum$temp,2), col="lightblue", pch=16)
points(d_tpc5b_sum$maxsl~jitter(d_tpc5b_sum$temp,2), col="yellow", pch=16)
#pretty similar by eye except light blue and yellow
#'black lower for higher growing ones





######plot tpcs#####
#b
#find means and variance
b_gr_summary <- b_gr %>% group_by(Treatment) %>%summarise_at(vars(sl), list(avg=mean, sd=sd))
#plot
ggplot(b_gr_summary, aes(x=Treatment, y=avg)) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.1) +
  geom_line() +
  geom_point() +
  geom_point(aes(x=23, y=0.525), color="#2BCC41") +
  geom_point(aes(x=29, y=0.362), color="red")


#d
#find means and variance
d_gr_summary <- d_gr %>% group_by(Treatment) %>%summarise_at(vars(sl), list(avg=mean, sd=sd))
#plot
ggplot(d_gr_summary, aes(x=Treatment, y=avg)) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.1) +
  geom_line() +
  geom_point() +
  geom_point(aes(x=23, y=0.356), color="#2BCC41") +
  geom_point(aes(x=27, y=0.176), color="red")

#d with points removed
#find means and variance
d_gr_cut_summary <- d_gr_cut %>% group_by(Treatment) %>%summarise_at(vars(sl), list(avg=mean, sd=sd))
#plot
ggplot(d_gr_cut_summary, aes(x=Treatment, y=avg)) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.1) +
  geom_line() +
  geom_point() +
  geom_point(aes(x=23, y=0.470), color="#2BCC41") +
  geom_point(aes(x=27, y=0.176), color="red")



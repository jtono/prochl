library(ggplot2)
library(dplyr)

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

#convert to days since begining
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

#####fit lines####
#basic regression line - by replicate
#for b
rep <- c()
temp <- c()
int <- c()
sl <- c()
for (i in unique(b_tpc$Treatment)){
  sub <- subset(b_tpc, Treatment==i)
  for (j in unique(sub$Rep.ID)){
    sub2 <- subset(sub, Rep.ID==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, j)
    temp <- c(temp, i)
  }
}
  
b_gr <- data.frame(temp, rep, int, sl)  
names(b_gr) <- c("Treatment","Rep.ID","int","sl")
  
p <- ggplot(data=b_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)
p
p + geom_abline(data=b_gr, aes(slope=sl, intercept=int, col=Rep.ID))

#for d
rep <- c()
temp <- c()
int <- c()
sl <- c()
for (i in unique(d_tpc$Treatment)){
  sub <- subset(d_tpc, Treatment==i)
  for (j in unique(sub$Rep.ID)){
    sub2 <- subset(sub, Rep.ID==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, j)
    temp <- c(temp, i)
  }
}

d_gr <- data.frame(temp, rep, int, sl)  
names(d_gr) <- c("Treatment","Rep.ID","int","sl")

p <- ggplot(data=d_tpc, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)
p
p + geom_abline(data=d_gr, aes(slope=sl, intercept=int, col=Rep.ID))

#cut last 2 days for 20, 23 and 25 
##FIND REASON FOR THIS!!!!
d_tpc_cut <- d_tpc[!(d_tpc$Treatment%in%c(20,23,25)&d_tpc$day%in%c(15,16)),]

rep <- c()
temp <- c()
int <- c()
sl <- c()
for (i in unique(d_tpc_cut$Treatment)){
  sub <- subset(d_tpc_cut, Treatment==i)
  for (j in unique(sub$Rep.ID)){
    sub2 <- subset(sub, Rep.ID==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, j)
    temp <- c(temp, i)
  }
}

d_gr_cut <- data.frame(temp, rep, int, sl)  
names(d_gr_cut) <- c("Treatment","Rep.ID","int","sl")

p <- ggplot(data=d_tpc_cut, aes(x=day, y=log(Actual.Cell.count), col=Rep.ID))+
  geom_point()+
  geom_line(aes(group=Rep.ID))+
  facet_wrap(~Treatment)
p
p + geom_abline(data=d_gr_cut, aes(slope=sl, intercept=int, col=Rep.ID))

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



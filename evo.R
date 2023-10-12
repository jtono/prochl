library(ggplot2)

#####getting data ready#####
#read in data
b_evo <- read.csv("Pc_B_Evo.csv", sep=",", as.is=TRUE, header=TRUE)
d_evo <- read.csv("Pc_D_Evo.csv", sep=",", as.is=TRUE, header=TRUE)

#strip extra rows
b_evo<- b_evo[1:288,]
d_evo <- d_evo[1:276,]

#strip extra columns
b_evo <- b_evo[,1:7]
d_evo <- d_evo[,c(1,2,3,4,5,7,10)]

#format date info
b_evo$Date <- as.Date(b_evo$Date, "%d/%m/%Y")
d_evo$Date <- as.Date(d_evo$Date, "%d/%m/%Y")

#convert to days since begining
b_evo$day <- as.numeric(round(difftime(b_evo$Date, b_evo$Date[1], units="days"), digits=0))
d_evo$day <- as.numeric(round(difftime(d_evo$Date, d_evo$Date[1], units="days"), digits=0))

#make column of cell counts that are numeric
b_evo$Cell.count <- sub(",","",b_evo$Cell.count)
b_evo$Cell.count <- as.numeric(b_evo$Cell.count)
#for b, NAs introduced
#should be 14 - checks out - these were originally nas

d_evo$Cell.count <- sub(",","",d_evo$Count)
d_evo$Cell.count <- as.numeric(d_evo$Cell.count)
#for d, NAs introduced
#should be 6 - checks out - these were originally nas

#####find growth rates####
#?I think Ruth told me that the cells were diluted to 1000 cells/10ul for each transfer
#?not sure what time 0 means then...
#?I'm assuming these cell counts are /10ul for this analysis
#to calculate growth rate: (log(Cell.count) - log(1000))/#days
#?does rep 1 in 23 correspond to rep 1 in 29 in some way? for now, assume no and make unique
#?need to check when dilutions not done/done - NAs right now stopping calculations

b_evo$Rep <- paste(b_evo$Rep.ID,b_evo$Treatment)
d_evo$Rep <- paste(d_evo$Rep.ID,d_evo$Treatment)

#calculate growth rates for b
b_evo$gr <- NA
for (i in unique(b_evo$Rep)){
  sub <- subset(b_evo, Rep==i)
  for (j in 2:length(sub$File.ID)){
    b_evo[b_evo$Rep==i,]$gr[j] <- (log(sub$Cell.count[j])-log(1000))/(sub$day[j]-sub$day[j-1])
  }
}

#if remove rows with NAs, assuming dilution not done on any of those days
b_evo_nona <- b_evo[-which(is.na(b_evo$Cell.count)),]

#calculate growth rates for b without NAs
b_evo_nona$gr <- NA
for (i in unique(b_evo_nona$Rep)){
  sub <- subset(b_evo_nona, Rep==i)
  for (j in 2:length(sub$File.ID)){
    b_evo_nona[b_evo_nona$Rep==i,]$gr[j] <- (log(sub$Cell.count[j])-log(1000))/(sub$day[j]-sub$day[j-1])
  }
}


#calculate growth rates for d
d_evo$gr <- NA
for (i in unique(d_evo$Rep)){
  sub <- subset(d_evo, Rep==i)
  for (j in 2:length(sub$File.ID)){
    d_evo[d_evo$Rep==i,]$gr[j] <- (log(sub$Cell.count[j])-log(1000))/(sub$day[j]-sub$day[j-1])
  }
}

#if remove rows with NAs, assuming dilution not done on any of those days
d_evo_nona <- d_evo[-which(is.na(d_evo$Cell.count)),]

#calculate growth rates for b without NAs
d_evo_nona$gr <- NA
for (i in unique(d_evo_nona$Rep)){
  sub <- subset(d_evo_nona, Rep==i)
  for (j in 2:length(sub$File.ID)){
    d_evo_nona[d_evo_nona$Rep==i,]$gr[j] <- (log(sub$Cell.count[j])-log(1000))/(sub$day[j]-sub$day[j-1])
  }
}

#####plotting#####
#plot b without removing nas
ggplot(b_evo, aes(x=day, y=gr, col=as.factor(Treatment))) +
  geom_line(aes(group=Rep))+
  geom_point()

#plot b with removing nas
ggplot(b_evo_nona, aes(x=day, y=gr, col=as.factor(Treatment))) +
  geom_line(aes(group=Rep))+
  geom_point()
#mostly ok but huge dip at one timepoint - day 125 = 2022-03-31

#plot d without removing nas
ggplot(d_evo, aes(x=day, y=gr, col=as.factor(Treatment))) +
  geom_line(aes(group=Rep))+
  geom_point()

#plot d with removing nas
ggplot(d_evo_nona, aes(x=day, y=gr, col=as.factor(Treatment))) +
  geom_line(aes(group=Rep))+
  geom_point()
#they have a big dip as well - day 118 = 2022-03-31 - same date
#seems to affect next points as well - related to covid gap? maybe stationary phase and took awhile to recover?

###try to find overall lines of best fit
#b - b_evo_nona
#subset by Treatment - 23
b_subt23 <- subset(b_evo_nona, Treatment==23)
#make version where remove weird points
b_subt23.rm <- b_subt23[-which(b_subt23$day==125),]
b_subt23.rm <- b_subt23.rm[-which(b_subt23.rm$day==131),]
b_subt23.rm <- b_subt23.rm[-which(b_subt23.rm$day==138),]
#fit a linear model
fit1 <- lm(gr~day, data=b_subt23)
fit1.rm <- lm(gr~day, data=b_subt23.rm)
#look at the fits
plot(b_subt23$day, b_subt23$gr, pch=19, col="red")
points(b_subt23.rm$day, b_subt23.rm$gr, pch=19, col="black")
x_axis <- seq(0,175,length=100)
lines(x_axis, predict(fit1, data.frame(day=x_axis)), col="red")
lines(x_axis, predict(fit1.rm, data.frame(day=x_axis)), col="black")

#find values
summary(fit1)$adj.r.squared
fit1$coefficients
summary(fit1.rm)$adj.r.squared
fit1.rm$coefficients

#subset by Treatment - 29
b_subt29 <- subset(b_evo_nona, Treatment==29)
#make version where remove weird points
b_subt29.rm <- b_subt29[-which(b_subt29$day==125),]
b_subt29.rm <- b_subt29.rm[-which(b_subt29.rm$day==131),]
b_subt29.rm <- b_subt29.rm[-which(b_subt29.rm$day==138),]
#fit a linear model
fit1 <- lm(gr~day, data=b_subt29)
fit1.rm <- lm(gr~day, data=b_subt29.rm)
#look at the fits
plot(b_subt29$day, b_subt29$gr, pch=19, col="red")
points(b_subt29.rm$day, b_subt29.rm$gr, pch=19, col="black")
x_axis <- seq(0,175,length=100)
lines(x_axis, predict(fit1, data.frame(day=x_axis)), col="red")
lines(x_axis, predict(fit1.rm, data.frame(day=x_axis)), col="black")

#find values
summary(fit1)$adj.r.squared
fit1$coefficients
summary(fit1.rm)$adj.r.squared
fit1.rm$coefficients

#d - d_evo_nona
#subset by Treatment - 23
d_subt23 <- subset(d_evo_nona, Treatment==23)
#make version where remove weird points
d_subt23.rm <- d_subt23[-which(d_subt23$day==118),]
d_subt23.rm <- d_subt23.rm[-which(d_subt23.rm$day==124),]
d_subt23.rm <- d_subt23.rm[-which(d_subt23.rm$day==131),]
#fit a linear model
fit1 <- lm(gr~day, data=d_subt23)
fit1.rm <- lm(gr~day, data=d_subt23.rm)
#look at the fits
plot(d_subt23$day, d_subt23$gr, pch=19, col="red")
points(d_subt23.rm$day, d_subt23.rm$gr, pch=19, col="black")
x_axis <- seq(0,175,length=100)
lines(x_axis, predict(fit1, data.frame(day=x_axis)), col="red")
lines(x_axis, predict(fit1.rm, data.frame(day=x_axis)), col="black")

#find values
summary(fit1)$adj.r.squared
fit1$coefficients
summary(fit1.rm)$adj.r.squared
fit1.rm$coefficients

#subset by Treatment - 27
d_subt27 <- subset(d_evo_nona, Treatment==27)
#make version where remove weird points
d_subt27.rm <- d_subt27[-which(d_subt27$day==118),]
d_subt27.rm <- d_subt27.rm[-which(d_subt27.rm$day==124),]
d_subt27.rm <- d_subt27.rm[-which(d_subt27.rm$day==131),]
#fit a linear model
fit1 <- lm(gr~day, data=d_subt27)
fit1.rm <- lm(gr~day, data=d_subt27.rm)
#look at the fits
plot(d_subt27$day, d_subt27$gr, pch=19, col="red")
points(d_subt27.rm$day, d_subt27.rm$gr, pch=19, col="black")
x_axis <- seq(0,175,length=100)
lines(x_axis, predict(fit1, data.frame(day=x_axis)), col="red")
lines(x_axis, predict(fit1.rm, data.frame(day=x_axis)), col="black")

#find values
summary(fit1)$adj.r.squared
fit1$coefficients
summary(fit1.rm)$adj.r.squared
fit1.rm$coefficients


#######compare beginning and end######

mean(d_subt27.rm[d_subt27.rm$day%in%unique(d_subt27.rm$day)[2:6],]$gr)
#0.2445429
mean(d_subt27.rm[d_subt27.rm$day%in%unique(d_subt27.rm$day)[15:19],]$gr)
#0.5135323

d_subt27_beg <- d_subt27[d_subt27$day%in%unique(d_subt27$day)[2:6],]
d_subt27_end <- d_subt27[d_subt27$day%in%unique(d_subt27$day)[(length(unique(d_subt27$day))-4):length(unique(d_subt27$day))],]
t.test(d_subt27_beg$gr, d_subt27_end$gr)
# Welch Two Sample t-test
#
# data:  d_subt27.rm_beg$gr and d_subt27.rm_end$gr
# t = -7.2127, df = 40.038, p-value = 9.448e-09
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3443608 -0.1936180
# sample estimates:
#   mean of x mean of y
# 0.2445429 0.5135323


d_subt23_beg <- d_subt23[d_subt23$day%in%unique(d_subt23$day)[2:6],]
d_subt23_end <- d_subt23[d_subt23$day%in%unique(d_subt23$day)[(length(unique(d_subt23$day))-4):length(unique(d_subt23$day))],]
t.test(d_subt23_beg$gr, d_subt23_end$gr)
#different in wrong direction.
#maybe do by sample?















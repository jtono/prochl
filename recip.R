library(ggplot2)
library(dplyr)

#####getting data ready#####
#read in data
b_recip <- read.csv("Pc_B_Recip.csv", sep=",", as.is=TRUE, header=TRUE)
d_recip <- read.csv("Pc_D_Recip.csv", sep=",", as.is=TRUE, header=TRUE)

#cut extra rows
b_recip <- b_recip[1:384,]
d_recip <- d_recip[1:384,]

#get rid of ones not numbers in d_recip and make numeric
d_recip <- d_recip[-which(d_recip$Actual.Cell.count=="#VALUE!"),]
d_recip$Actual.Cell.count <- as.numeric(d_recip$Actual.Cell.count)
#get rid of 0s in d_recip
d_recip<-d_recip[-which(d_recip$Actual.Cell.count==0),]


#extract date information
b_recip$Date <- NA
for (i in 1:length(b_recip$File.ID)){
  b_recip$Date[i] <- strsplit(b_recip$File.ID[i], "_")[[1]][4]
}

d_recip$Date <- NA
for (i in 1:length(d_recip$File.ID)){
  d_recip$Date[i] <- strsplit(d_recip$File.ID[i], "_")[[1]][4]
}

b_recip$Date <- format(as.Date(b_recip$Date, "%d%m%Y"), "20%y-%m-%d")
d_recip$Date <- format(as.Date(d_recip$Date, "%d%m%Y"), "20%y-%m-%d")

#convert to days since beginning
b_recip$day <- as.numeric(round(difftime(b_recip$Date, b_recip$Date[1], units="days"), digits=0))
d_recip$day <- as.numeric(round(difftime(d_recip$Date, d_recip$Date[1], units="days"), digits=0))

#figure out which temp evolved in
#in col Evo.Temp
#add to column Evo.Temp when control (currently blank)
for (i in 1:length(d_recip$File.ID)){
  if(is.na(d_recip$Evo.Temp[i])){
    d_recip$Evo.Temp[i] <- as.numeric(gsub("c","",d_recip$Assay.Temp[i]))
  }
}

#add to column Evo.Temp when control (currently blank)
for (i in 1:length(b_recip$File.ID)){
  if(is.na(b_recip$Evo.Temp[i])){
    b_recip$Evo.Temp[i] <- as.numeric(gsub("c","",b_recip$Assay.Temp[i]))
  }
}

#make new replicate ids
b_recip$rep <- paste(b_recip$Rep.ID,b_recip$Evo.Temp)
d_recip$rep <- paste(d_recip$Rep.ID,d_recip$Evo.Temp)

#####plot growth curves######
ggplot(data=b_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)

ggplot(data=d_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)

#####fit lines####
#basic regression line - by replicate
#for b
rep <- c()
Assay.Temp <- c()
int <- c()
sl <- c()
Evo.Temp <- c()
for (i in unique(b_recip$rep)){
  sub <- subset(b_recip, rep==i)
  for (j in unique(sub$Assay.Temp)){
    sub2 <- subset(sub, Assay.Temp==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, i)
    Assay.Temp <- c(Assay.Temp, j)
    Evo.Temp <- c(Evo.Temp, sub2$Evo.Temp[1])
  }
}

b_gr <- data.frame(Assay.Temp, rep, int, sl, Evo.Temp)
names(b_gr) <- c("Assay.Temp","Rep","int","sl", "Evo.Temp")

p <- ggplot(data=b_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)
p
p + geom_abline(data=b_gr, aes(slope=sl, intercept=int, col=Evo.Temp))

#for d
rep <- c()
Assay.Temp <- c()
int <- c()
sl <- c()
Evo.Temp <- c()
for (i in unique(d_recip$rep)){
  sub <- subset(d_recip, rep==i)
  for (j in unique(sub$Assay.Temp)){
    sub2 <- subset(sub, Assay.Temp==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, i)
    Assay.Temp <- c(Assay.Temp, j)
    Evo.Temp <- c(Evo.Temp, sub2$Evo.Temp[1])
  }
}

d_gr <- data.frame(Assay.Temp, rep, int, sl, Evo.Temp)
names(d_gr) <- c("Assay.Temp","Rep","int","sl", "Evo.Temp")

p <- ggplot(data=d_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)
p
p + geom_abline(data=d_gr, aes(slope=sl, intercept=int, col=Evo.Temp))


######plot grs#####
#b
b_gr$Treatment <- paste(b_gr$Evo.Temp, "@",b_gr$Assay.Temp)
p <- ggplot(b_gr, aes(x=Treatment, y=sl, col=Evo.Temp)) +
  geom_boxplot()
p

#d
d_gr$Treatment <- paste(d_gr$Evo.Temp, "@",d_gr$Assay.Temp)
p <- ggplot(d_gr, aes(x=Treatment, y=sl, col=Evo.Temp)) +
  geom_boxplot()
p

###########stats###############
mod1 <- lm(b_gr$sl~d_gr$Assay.Temp*d_gr$Evo.Temp)
summary(mod1)
anova(mod1)
mod2 <- lm(b_gr$sl~d_gr$Treatment)
summary(mod2)
anova(mod2)
mod3 <- lm(d_gr$sl~d_gr$Assay.Temp*d_gr$Evo.Temp)
summary(mod3)
anova(mod3)
mod4 <- lm(d_gr$sl~d_gr$Treatment)
summary(mod4)
anova(mod4)

t.test(b_gr$)

########check output######
b_o23 <- subset(b_recip, Orig.Temp=="23c")
b_o23_con <- subset(b_o23, Treatment=="Control")
plot(log(b_o23_con$Actual.Cell.count)~b_o23_con$day, type="b")
b_o23_recip <- subset(b_o23, Treatment=="recip")
plot(log(b_o23_recip$Actual.Cell.count)~b_o23_recip$day, type="b")
b_o29 <- subset(b_recip, Orig.Temp=="29c")
b_o29_con <- subset(b_o29, Treatment=="Control")
plot(log(b_o29_con$Actual.Cell.count)~b_o29_con$day, type="b")
b_o29_recip <- subset(b_o29, Treatment=="recip")
plot(log(b_o29_recip$Actual.Cell.count)~b_o29_recip$day, type="b")

d_o23 <- subset(d_recip, Orig.Temp=="23c")
d_o23_con <- subset(d_o23, Treatment=="Control")
plot(log(d_o23_con$Actual.Cell.count)~d_o23_con$day, type="b")
d_o23_recip <- subset(d_o23, Treatment=="recip")
plot(log(d_o23_recip$Actual.Cell.count)~d_o23_recip$day, type="b")
d_o27 <- subset(d_recip, Orig.Temp=="27c")
d_o27_con <- subset(d_o27, Treatment=="Control")
plot(log(d_o27_con$Actual.Cell.count)~d_o27_con$day)
d_o27_recip <- subset(d_o27, Treatment=="recip")
plot(log(d_o27_recip$Actual.Cell.count)~d_o27_recip$day, type="b")
#same sample from D low, so maybe mixup in which is which???


b_15 <- subset(b_recip, day==15)
p <- ggplot(b_15, aes(x=Treatment, y=Actual.Cell.count, col=Orig.Temp)) +
  geom_boxplot()
p

b_15$og.gr <- paste(b_15$Orig.Temp,b_15$grtemp)
p <- ggplot(b_15, aes(x=og.gr, y=Actual.Cell.count, col=grtemp)) +
  geom_boxplot()
p

b_0 <- subset(b_recip, day==0)
p <- ggplot(b_0, aes(x=Treatment, y=Actual.Cell.count, col=Orig.Temp)) +
  geom_boxplot()
p

b_0$og.gr <- paste(b_0$Orig.Temp,b_0$grtemp)
p <- ggplot(b_0, aes(x=og.gr, y=Actual.Cell.count, col=grtemp)) +
  geom_boxplot()
p

d_15 <- subset(d_recip, day==15)
p <- ggplot(d_15, aes(x=Treatment, y=Actual.Cell.count, col=Orig.Temp)) +
  geom_boxplot()
p

d_15$og.gr <- paste(d_15$Orig.Temp,d_15$grtemp)
p <- ggplot(d_15, aes(x=og.gr, y=Actual.Cell.count, col=grtemp)) +
  geom_boxplot()
p

d_1 <- subset(d_recip, day==1)
p <- ggplot(d_1, aes(x=Treatment, y=Actual.Cell.count, col=Orig.Temp)) +
  geom_boxplot()
p

d_1$og.gr <- paste(d_1$Orig.Temp,d_1$grtemp)
p <- ggplot(d_1, aes(x=og.gr, y=Actual.Cell.count, col=grtemp)) +
  geom_boxplot()
p




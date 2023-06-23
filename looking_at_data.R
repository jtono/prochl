
library(ggplot2)

#read in data
b_evo <- read.csv("Pc_B_Evo.csv", sep=",", as.is=TRUE, header=TRUE)
b_recip <- read.csv("Pc_B_Recip.csv", sep=",", as.is=TRUE, header=TRUE)
b_tpc <- read.csv("Pc_B_TPC.csv", sep=",", as.is=TRUE, header=TRUE)

d_evo <- read.csv("Pc_D_Evo.csv", sep=",", as.is=TRUE, header=TRUE)
d_recip <- read.csv("Pc_D_Recip.csv", sep=",", as.is=TRUE, header=TRUE)
d_tpc <- read.csv("Pc_D_TPC.csv", sep=",", as.is=TRUE, header=TRUE)

#cut extra rows from b_tpc
b_tpc <- b_tpc[1:408,]

#quick plot of data
ggplot(b_tpc, aes(x=Treatment, y=Actual.Cell.count, color=Timepoint)) + geom_jitter()

#get rid of ones not numbers in d
which(d_tpc$Actual.Cell.count=="#VALUE!")
d_tpc <- d_tpc[-which(d_tpc$Actual.Cell.count=="#VALUE!"),]
d_tpc$Actual.Cell.count <- as.numeric(d_tpc$Actual.Cell.count)

d_tpc$Treatment <- sub("c","",d_tpc$Treatment)
d_tpc$Treatment <- as.numeric(d_tpc$Treatment)

b_tpc$Treatment <- sub("c","",b_tpc$Treatment)
b_tpc$Treatment <- as.numeric(b_tpc$Treatment)


ggplot(d_tpc, aes(x=Treatment, y=Actual.Cell.count, color=Timepoint)) + geom_jitter()

#make factors of timepoints and plot
d_tpc$Timepoint <- as.factor(d_tpc$Timepoint)
d_tpc$Timepoint <- factor(d_tpc$Timepoint, c("T0","T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16"))

ggplot(data=d_tpc, aes(x=Timepoint, y=log(Actual.Cell.count)))+
  geom_jitter()+
  facet_wrap(~Treatment)


b_tpc$Timepoint <- as.factor(b_tpc$Timepoint)
b_tpc$Timepoint <- factor(b_tpc$Timepoint, c("T0","T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16"))

ggplot(data=b_tpc, aes(x=Timepoint, y=log(Actual.Cell.count)))+
  geom_jitter()+
  facet_wrap(~Treatment)


#####look at evo
View(b_evo)

b_evo<- b_evo[1:288,]
b_evo$Date <- as.Date(b_evo$Date, "%d/%m/%Y")
b_evo$Cell.count <- sub(",","",b_evo$Cell.count)
b_evo$Cell.count <- as.numeric(b_evo$Cell.count)

plot(b_evo$Cell.count~b_evo$Date, xaxt="n", col=b_evo$Treatment)
axis(1, b_evo$Date, format(b_evo$Date, "%Y%b%d"), cex.axis = .7, las=2)



d_evo$Date <- as.Date(d_evo$Date, "%d/%m/%Y")
d_evo$Count <- sub(",","",d_evo$Count)
d_evo$Count <- as.numeric(d_evo$Count)

plot(d_evo$Count~d_evo$Date, xaxt="n", col=d_evo$Treatment)
axis(1, d_evo$Date, format(d_evo$Date, "%Y%b%d"), cex.axis = .7, las=2)

###left off here####
ggplot(data=b_evo, aes(x=Date, y=Cell.count))+
  geom_jitter()

as.numeric(b_evo$Cell.count)



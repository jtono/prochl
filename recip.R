library(ggplot2)
library(dplyr)
library("multcompView")
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
    d_recip$Evo.Temp[i] <- gsub("c","",d_recip$Assay.Temp[i])
  }
}

#add to column Evo.Temp when control (currently blank)
for (i in 1:length(b_recip$File.ID)){
  if(is.na(b_recip$Evo.Temp[i])){
    b_recip$Evo.Temp[i] <- gsub("c","",b_recip$Assay.Temp[i])
  }
}

#make both Evo.Temp and Assay.Temp factors
b_recip$Evo.Temp <- as.factor(b_recip$Evo.Temp)
b_recip$Assay.Temp <- as.factor(b_recip$Assay.Temp)
d_recip$Evo.Temp <- as.factor(d_recip$Evo.Temp)
d_recip$Assay.Temp <- as.factor(d_recip$Assay.Temp)

#make new replicate ids
b_recip$rep <- paste(b_recip$Rep.ID,b_recip$Evo.Temp)
d_recip$rep <- paste(d_recip$Rep.ID,d_recip$Evo.Temp)

#####plot growth curves######
ggplot(data=b_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  geom_vline(xintercept=10.5, linetype="dotted")+
  facet_wrap(~Assay.Temp)+
  labs(x = "Day",
       y = "log(Actual Cell Count)",
       title="eMIT9312")

ggplot(data=d_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  geom_vline(xintercept=10.5, linetype="dotted")+
  facet_wrap(~Assay.Temp)+
  labs(x = "Day",
       y = "log(Actual Cell Count)",
       title="eMED4")

#####fit lines####
#basic regression line - by replicate
#for b
rep <- c()
A.Temp <- c()
int <- c()
sl <- c()
E.Temp <- c()
for (i in unique(b_recip$rep)){
  sub <- subset(b_recip, rep==i)
  for (j in unique(sub$Assay.Temp)){
    sub2 <- subset(sub, Assay.Temp==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, i)
    A.Temp <- c(A.Temp, j)
    E.Temp <- c(E.Temp, as.character(sub2$Evo.Temp[1]))
  }
}

b_gr <- data.frame(A.Temp, rep, int, sl, E.Temp)
names(b_gr) <- c("Assay.Temp","Rep","int","sl", "Evo.Temp")
b_gr$Evo.Temp <- as.factor(b_gr$Evo.Temp)

p <- ggplot(data=b_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)
p
p + geom_abline(data=b_gr, aes(slope=sl, intercept=int, col=Evo.Temp))

#for d
rep <- c()
A.Temp <- c()
int <- c()
sl <- c()
E.Temp <- c()
for (i in unique(d_recip$rep)){
  sub <- subset(d_recip, rep==i)
  for (j in unique(sub$Assay.Temp)){
    sub2 <- subset(sub, Assay.Temp==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, i)
    A.Temp <- c(A.Temp, j)
    E.Temp <- c(E.Temp, as.character(sub2$Evo.Temp[1]))
  }
}

d_gr <- data.frame(A.Temp, rep, int, sl, E.Temp)
names(d_gr) <- c("Assay.Temp","Rep","int","sl", "Evo.Temp")
d_gr$Evo.Temp <- as.factor(d_gr$Evo.Temp)

p <- ggplot(data=d_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)
p
p + geom_abline(data=d_gr, aes(slope=sl, intercept=int, col=Evo.Temp))

#basic regression line - by replicate - cutoff at 10
b_recip10 <- subset(b_recip, day<10.5)
#for b
rep <- c()
A.Temp <- c()
int <- c()
sl <- c()
E.Temp <- c()
for (i in unique(b_recip10$rep)){
  sub <- subset(b_recip10, rep==i)
  for (j in unique(sub$Assay.Temp)){
    sub2 <- subset(sub, Assay.Temp==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, i)
    A.Temp <- c(A.Temp, j)
    E.Temp <- c(E.Temp, as.character(sub2$Evo.Temp[1]))
  }
}

b_gr10 <- data.frame(A.Temp, rep, int, sl, E.Temp)
names(b_gr10) <- c("Assay.Temp","Rep","int","sl", "Evo.Temp")
b_gr10$Evo.Temp <- as.factor(b_gr10$Evo.Temp)

p <- ggplot(data=b_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)
p
p + geom_abline(data=b_gr10, aes(slope=sl, intercept=int, col=Evo.Temp))
p + geom_abline(data=b_gr, aes(slope=sl, intercept=int, col=Evo.Temp))

d_recip10 <- subset(d_recip, day<10.5)
#for d
rep <- c()
A.Temp <- c()
int <- c()
sl <- c()
E.Temp <- c()
for (i in unique(d_recip10$rep)){
  sub <- subset(d_recip10, rep==i)
  for (j in unique(sub$Assay.Temp)){
    sub2 <- subset(sub, Assay.Temp==j)
    mod <- lm(log(Actual.Cell.count)~day, sub2)
    int <- c(int, mod$coefficients[1])
    sl <- c(sl, mod$coefficients[2])
    rep <- c(rep, i)
    A.Temp <- c(A.Temp, j)
    E.Temp <- c(E.Temp, as.character(sub2$Evo.Temp[1]))
  }
}

d_gr10 <- data.frame(A.Temp, rep, int, sl, E.Temp)
names(d_gr10) <- c("Assay.Temp","Rep","int","sl", "Evo.Temp")
d_gr10$Evo.Temp <- as.factor(d_gr10$Evo.Temp)

p <- ggplot(data=d_recip, aes(x=day, y=log(Actual.Cell.count), col=Evo.Temp))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~Assay.Temp)
p
p + geom_abline(data=d_gr10, aes(slope=sl, intercept=int, col=Evo.Temp))
p + geom_abline(data=d_gr, aes(slope=sl, intercept=int, col=Evo.Temp))




###########stats###############
b_gr$Treatment <- paste(b_gr$Assay.Temp, ":",b_gr$Evo.Temp, sep="")
b_gr10$Treatment <- paste(b_gr10$Assay.Temp, ":",b_gr10$Evo.Temp, sep="")
d_gr$Treatment <- paste(d_gr$Assay.Temp, ":",d_gr$Evo.Temp, sep="")
d_gr10$Treatment <- paste(d_gr10$Assay.Temp, ":",d_gr10$Evo.Temp, sep="")

aov.model1 <- aov(sl~Assay.Temp*Evo.Temp, data=b_gr)
summary(aov.model1)
aov.model2 <- aov(sl~Treatment, data=b_gr)
summary(aov.model2)
b_treat_tk <- TukeyHSD(aov.model2, conf.level=.95)
plot(b_treat_tk)


aov.model3 <- aov(sl~Assay.Temp*Evo.Temp, data=d_gr)
summary(aov.model3)
aov.model4 <- aov(sl~Treatment, data=d_gr)
summary(aov.model4)
d_treat_tk <- TukeyHSD(aov.model4, conf.level=.95)
plot(d_treat_tk)


aov.modelb10 <- aov(sl~Assay.Temp*Evo.Temp, data=b_gr10)
summary(aov.modelb10)
aov.model210 <- aov(sl~Treatment, data=b_gr10)
summary(aov.model210)
b_treat_tk10 <- TukeyHSD(aov.modelb10, conf.level=.95)
plot(b_treat_tk10)
#write.csv(b_treat_tk10$`Assay.Temp:Evo.Temp`, "b_treat_tk10.csv")


aov.modeld10 <- aov(sl~Assay.Temp*Evo.Temp, data=d_gr10)
summary(aov.modeld10)
aov.model410 <- aov(sl~Treatment, data=d_gr10)
summary(aov.model410)
d_treat_tk10 <- TukeyHSD(aov.modeld10, conf.level=.95)
plot(d_treat_tk10)
#write.csv(d_treat_tk10$`Assay.Temp:Evo.Temp`, "d_treat_tk10.csv")





######plot grs - boxplots with tukey letters#####

#b
p <- ggplot(b_gr, aes(x=Assay.Temp, y=sl, col=Evo.Temp)) +
  geom_boxplot()+
  labs(x = 'Assay Temperature',
       y = 'Growth rate',
       title = 'eMIT9312')
p

letters.df <- data.frame(multcompLetters(TukeyHSD(aov.modelb10, conf.level=.95)$`Assay.Temp:Evo.Temp`[,4])$Letters)
colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Treatment <- rownames(letters.df) #Create column based on rownames
placement <- b_gr10 %>% #We want to create a dataframe to assign the letter position.
  group_by(Assay.Temp,Evo.Temp) %>%
  summarise(quantile(sl)[4])
placement$Treatment <- paste(placement$Assay.Temp, ":", placement$Evo.Temp, sep="")
colnames(placement)[3] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) #Merge dataframes

b_gr10 <- left_join(b_gr10, letters.df, by="Treatment")


b_box <- ggplot(b_gr10, aes(x=Assay.Temp.x, y=sl, group=interaction(Assay.Temp.x,Evo.Temp.x))) +
  geom_boxplot(aes(col=Evo.Temp.x))+
  labs(x = 'Assay Temperature',
       y = 'Growth rate',
       title = 'eMIT9312')+
  theme_bw() +
  scale_color_manual(values=c("royalblue1","red4"), name = "Treatment", labels = c("Control", "Warmed")) +
  geom_text(data = b_gr10, aes(x = Assay.Temp.x, y = Placement.Value, label = Letter), size = 4, color = "black" , hjust = -0.5, vjust = -0.5, fontface = "bold", position=position_dodge(0.75))+
  scale_x_discrete(breaks=c("23c","29c"),labels=c("23","29"))+
  ylim(0,0.7)
b_box


#d
p <- ggplot(d_gr, aes(x=Assay.Temp, y=sl, col=Evo.Temp)) +
  geom_boxplot()+
  labs(x = 'Assay Temperature',
       y = 'Growth rate',
       title = 'eMED4')
p

letters.df <- data.frame(multcompLetters(TukeyHSD(aov.modeld10, conf.level=.95)$`Assay.Temp:Evo.Temp`[,4])$Letters)
colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Treatment <- rownames(letters.df) #Create column based on rownames
placement <- d_gr10 %>% #We want to create a dataframe to assign the letter position.
  group_by(Assay.Temp,Evo.Temp) %>%
  summarise(quantile(sl)[4])
placement$Treatment <- paste(placement$Assay.Temp, ":", placement$Evo.Temp, sep="")
colnames(placement)[3] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) #Merge dataframes

d_gr10 <- left_join(d_gr10, letters.df, by="Treatment")


d_box <- ggplot(d_gr10, aes(x=Assay.Temp.x, y=sl, group=interaction(Assay.Temp.x,Evo.Temp.x))) +
  geom_boxplot(aes(col=Evo.Temp.x))+
  labs(x = 'Assay Temperature',
       y = 'Growth rate',
       title = 'eMED4')+
  theme_bw() +
  scale_color_manual(values=c("royalblue1","red4"), name = "Treatment", labels = c("Control", "Warmed")) +
  geom_text(data = d_gr10, aes(x = Assay.Temp.x, y = Placement.Value, label = Letter), size = 4, color = "black" , hjust = -0.5, vjust = -0.5, fontface = "bold", position=position_dodge(0.75))+
  scale_x_discrete(breaks=c("23c","27c"),labels=c("23","27"))+
  ylim(0,0.7)
d_box

#combine two figs
comb_box <- ggarrange(d_box + rremove("ylab"), b_box + rremove("ylab"), ncol=2, nrow=1, legend = "right", common.legend=TRUE)
annotate_figure(comb_box, left = text_grob("Growth rate", rot=90))

########find one ProD replicate with low growth - rerun analysis with omitted#########
d_gr10[which(d_gr10$sl<0.2),]
#Assay.Temp   Rep      int        sl Evo.Temp Treatment
#2         27c R1 23 7.311519 0.1960166       23  23 @ 27c
#23        27c R6 27 6.882175 0.1583683       27  27 @ 27c

d_gr10[d_gr10$Rep=="R6 27",]
# Assay.Temp   Rep      int        sl Evo.Temp Treatment
#23        27c R6 27 6.882175 0.1583683       27  27 @ 27c
#24        23c R6 27 6.657928 0.3936455       27  27 @ 23c
#this replicate lower in both

d_gr10_noout <- d_gr10[-which(d_gr10$Rep=="R6 27"),]

p <- ggplot(d_gr10_noout, aes(x=Treatment, y=sl, col=Evo.Temp)) +
  geom_boxplot()+
  labs(x = 'Treatment',
       y = 'Growth rate',
       title = 'ProD')
p

aov.model310.noout <- aov(sl~Assay.Temp*Evo.Temp, data=d_gr10_noout)
summary(aov.model310.noout)
aov.model410.noout <- aov(sl~Treatment, data=d_gr10_noout)
summary(aov.model410.noout)
d_treat_tk10_noout <- TukeyHSD(aov.model410.noout, conf.level=.95)
plot(d_treat_tk10_noout)
write.csv(d_treat_tk10_noout$Treatment, "d_treat_tk10_noout.csv")



########check output - old???######
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




library(ggplot2)
library(nlme)
library(lmerTest)
library(visreg)
library(MetBrewer)
library("flextable")
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
title(main="ProB at 23")

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
title(main="ProB at 29")

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
title(main="ProD at 23")

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
title(main="ProD at 27")

#find values
summary(fit1)$adj.r.squared
fit1$coefficients
summary(fit1.rm)$adj.r.squared
fit1.rm$coefficients

#####mixed effects models - b #####
#b - b_evo_nona
#make Treatment, Rep into factors
b_evo_nona$Treatment <- as.factor(b_evo_nona$Treatment)
b_evo_nona$Rep <- as.factor(b_evo_nona$Rep)
#make version where remove weird points
b_evo_rm <- b_evo_nona[-which(b_evo_nona$day%in%c(125,131,138)),]

#fit a mixed effects model
#full model with random effects (slope and intercept) and no temporal correlation fixing
b_rm <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=b_evo_rm, na.action=na.exclude)
#full model with random effects (intercept only) and no temporal correlation fixing
b_rmI <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude)
#no random effects
b_rm0 <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude)


#full model with random effects (slope and intercept) and AR1 correlation
#if add form with day, singular convergence
b_rm_AR1 <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corAR1())

#full model with random effects (intercept) and AR1 correlation
b_rmI_AR1 <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corAR1())
b_rmI_AR1b <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day))
b_rmI_AR1c <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day|Rep))

#no random effects and AR1 correlation
b_rm0_AR1 <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corAR1())
b_rm0_AR1a <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corAR1(form=~1|Rep))
b_rm0_AR1c <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day|Rep))

#full model with random effects (slope and intercept) and ARMA correlation
#if add form with day, singular convergence
b_rm_ARMA <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1))

#full model with random effects (intercept) and ARMA correlation
b_rmI_ARMA <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1))
b_rmI_ARMAb <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(form=~day,p=1,q=1))
b_rmI_ARMAc <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(form=~day|Rep,p=1,q=1))

#no random effects and ARMA correlation
b_rm0_ARMA <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1))
b_rm0_ARMAa <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1,form=~1|Rep))
b_rm0_ARMAc <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1, form=~day|Rep))

#full model with random effects (slope and intercept) and CAR1 correlation
b_rm_CAR1 <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1())
b_rm_CAR1b <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day))
b_rm_CAR1c <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day|Rep))

#full model with random effects (intercept) and CAR1 correlation
b_rmI_CAR1 <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1())
b_rmI_CAR1b <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day))
b_rmI_CAR1c <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day|Rep))

#no random effects and CAR1 correlation
b_rm0_CAR1 <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1())
b_rm0_CAR1a <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~1|Rep))
b_rm0_CAR1c <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day|Rep))

#full model with random effects (slope and intercept) and CompSymm correlation
#won't converge

#full model with random effects (intercept) and CompSymm correlation
b_rmI_CompSymm <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm())
b_rmI_CompSymmb <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day))
b_rmI_CompSymmc <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day|Rep))

#no random effects and CompSymm correlation
b_rm0_CompSymm <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm())
b_rm0_CompSymma <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~1|Rep))
b_rm0_CompSymmb <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day))
b_rm0_CompSymmc <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day|Rep))

#corSymm has convergence problems in all tried
#this one wouldn't finish running
#b_rm0_Symm <- gls(gr ~ day*Treatment, data=b_evo_rm, na.action=na.exclude, correlation=corSymm())


aic_comp_b <- AIC(b_rm,b_rmI,b_rm0,b_rm_AR1,b_rmI_AR1,b_rmI_AR1b,b_rmI_AR1c,b_rm0_AR1,b_rm0_AR1a,b_rm0_AR1c,b_rm_ARMA,b_rmI_ARMA,b_rmI_ARMAb,b_rmI_ARMAc,b_rm0_ARMA,b_rm0_ARMAa,b_rm0_ARMAc,b_rm_CAR1,b_rm_CAR1b,b_rm_CAR1c,b_rmI_CAR1,b_rmI_CAR1b,b_rmI_CAR1c,b_rm0_CAR1,b_rm0_CAR1a,b_rm0_CAR1c,b_rmI_CompSymm,b_rmI_CompSymmb,b_rmI_CompSymmc,b_rm0_CompSymm,b_rm0_CompSymma,b_rm0_CompSymmb,b_rm0_CompSymmc)

which(aic_comp_b$AIC == min(aic_comp_b$AIC))
which(aic_comp_b$AIC < -600)
aic_comp_b[c(8,15,24),]
aic_comp_b[15,]
#         df       AIC
#b_rm0_ARMA  7 -643.10244

anov_comp_b <- anova(b_rm,b_rmI,b_rm0,b_rm_AR1,b_rmI_AR1,b_rmI_AR1b,b_rmI_AR1c,b_rm0_AR1,b_rm0_AR1a,b_rm0_AR1c,b_rm_ARMA,b_rmI_ARMA,b_rmI_ARMAb,b_rmI_ARMAc,b_rm0_ARMA,b_rm0_ARMAa,b_rm0_ARMAc,b_rm_CAR1,b_rm_CAR1b,b_rm_CAR1c,b_rmI_CAR1,b_rmI_CAR1b,b_rmI_CAR1c,b_rm0_CAR1,b_rm0_CAR1a,b_rm0_CAR1c,b_rmI_CompSymm,b_rmI_CompSymmb,b_rmI_CompSymmc,b_rm0_CompSymm,b_rm0_CompSymma,b_rm0_CompSymmb,b_rm0_CompSymmc, test=TRUE)

#save anova table output
save_as_docx(flextable(anov_comp_b),path="tables/tableSup1_b.docx")

which(anov_comp_b$logLik==max(anov_comp_b$logLik))
anov_comp_b[15,]
# Model df       AIC       BIC   logLik     Test  L.Ratio p-value
#b_rm0_ARMA    15  7 -643.1024 -619.2522 328.5512 14 vs 15 136.3066  <.0001
which(anov_comp_b$logLik>300)
anov_comp_b[c(8,15,24),]
anova(b_rm0_AR1, b_rm0_ARMA, b_rm0_CAR1)

anova(b_rm,b_rmI,b_rm0)
#all same
anova(b_rm,b_rm_AR1,b_rm_ARMA,b_rm_CAR1,b_rm_CAR1b,b_rm_CAR1c)
#b_rm_ARMA best
anova(b_rmI,b_rmI_AR1,b_rmI_AR1b,b_rmI_AR1c,b_rmI_ARMA,b_rmI_ARMAb,b_rmI_ARMAc,b_rmI_CAR1,b_rmI_CAR1b,b_rmI_CAR1c,b_rmI_CompSymm,b_rmI_CompSymmb,b_rmI_CompSymmc)
#b_rmI_ARMA best
anova(b_rm0,b_rm0_AR1,b_rm0_AR1a, b_rm0_AR1c,b_rm0_ARMA,b_rm0_ARMAa,b_rm0_ARMAc,b_rm0_CAR1,b_rm0_CAR1a,b_rm0_CAR1c,b_rm0_CompSymm,b_rm0_CompSymma,b_rm0_CompSymmb,b_rm0_CompSymmc)
#b_rm0_ARMA best
anova(b_rm_AR1,b_rmI_AR1,b_rmI_AR1b,b_rmI_AR1c,b_rm0_AR1,b_rm0_AR1a,b_rm0_AR1c)
#b_rm0_AR1 best
anova(b_rm_ARMA,b_rmI_ARMA,b_rmI_ARMAb,b_rmI_ARMAc,b_rm0_ARMA,b_rm0_ARMAa,b_rm0_ARMAc)
#b_rm0_ARMA best
anova(b_rm_CAR1,b_rm_CAR1b,b_rm_CAR1c,b_rmI_CAR1,b_rmI_CAR1b,b_rmI_CAR1c,b_rm0_CAR1,b_rm0_CAR1a,b_rm0_CAR1c)
#b_rm0_CAR1 best
anova(b_rmI_CompSymm,b_rmI_CompSymmb,b_rmI_CompSymmc,b_rm0_CompSymm,b_rm0_CompSymma,b_rm0_CompSymmb,b_rm0_CompSymmc)
#many same/similar, including b_rm0_CompSymmc

#b_rm0_ARMA best. but if want random effects, look at one of the models with that.

#checking plots
plot(b_rm0_ARMA)

plot(b_rm0, resid(., type = "normalized") ~ day | Rep, abline = 0)
plot(b_rm0_ARMA, resid(., type = "normalized") ~ day | Rep, abline = 0)
plot(b_rm0, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
plot(b_rm0_ARMA, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
#these don't work for rm0
plot(ACF(b_rm0, resType = "normalized"), alpha=0.05)
plot(ACF(b_rm0_ARMA, resType = "normalized"), alpha=0.05)

#look at best models
summary(b_rm0_ARMA)
#VarCorr(b_rm0_ARMA)
confint(b_rm0_ARMA, parm=c("day","Treatment29"))   # lmer
intervals(b_rm0_ARMA, which="all", level=0.95)
anova(b_rm0_ARMA)
#fixef(b_rm0_ARMA)
#ranef(b_lme4)

#plot - could also use , type="contrast"
visreg(b_rm0_ARMA, "day", by="Treatment", gg=TRUE, overlay=TRUE)
visreg(b_rm0_ARMA, "day", by="Treatment", overlay=TRUE, ylab="Growth rate", xlab="Day", gg=TRUE)+
  theme_bw() +
  scale_color_manual(values=c("blue","red")) +
  scale_fill_manual(values=c("blue","red")) +
  labs(title="ProB")

#####mixed effects models - d#####
#d - d_evo_nona
#make Treatment, Rep into factors
d_evo_nona$Treatment <- as.factor(d_evo_nona$Treatment)
d_evo_nona$Rep <- as.factor(d_evo_nona$Rep)
#make version where remove weird points - same dates as b
d_evo_rm <- d_evo_nona[-which(d_evo_nona$day%in%c(118,124,131)),]


#fit a mixed effects model
#full model with random effects (slope and intercept) and no temporal correlation fixing
#d_rm <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude)
#won't converge
#full model with random effects (intercept only) and no temporal correlation fixing
d_rmI <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude)
#no random effects
d_rm0 <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude)


#full model with random effects (slope and intercept) and AR1 correlation
d_rm_AR1 <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corAR1())
#d_rm_AR1b <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day))
#d_rm_AR1c <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day|Rep))

#full model with random effects (intercept) and AR1 correlation
d_rmI_AR1 <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corAR1())
d_rmI_AR1b <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day))
d_rmI_AR1c <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day|Rep))

#no random effects and AR1 correlation
d_rm0_AR1 <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corAR1())
d_rm0_AR1a <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~1|Rep))
#d_rm0_AR1b <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day))
#doesn't work
d_rm0_AR1c <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corAR1(form=~day|Rep))

#full model with random effects (slope and intercept) and ARMA correlation
d_rm_ARMA <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1))
#d_rm_ARMAb <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1,form=~day))
#doesn't converge
#d_rm_ARMAc <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1,form=~day|Rep))
#doesn't converge

#full model with random effects (intercept) and ARMA correlation
d_rmI_ARMA <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1))
d_rmI_ARMAb <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(form=~day,p=1,q=1))
d_rmI_ARMAc <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(form=~day|Rep,p=1,q=1))

#no random effects and ARMA correlation
d_rm0_ARMA <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1))
d_rm0_ARMAa <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1,form=~1|Rep))
#d_rm0_ARMAb <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1,form=~day))
#doesn't work
d_rm0_ARMAc <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corARMA(p=1,q=1, form=~day|Rep))

#full model with random effects (slope and intercept) and CAR1 correlation
d_rm_CAR1 <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1())
#d_rm_CAR1b <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day))
#won't converge
#d_rm_CAR1c <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day|Rep))
#won't converge

#full model with random effects (intercept) and CAR1 correlation
d_rmI_CAR1 <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1())
d_rmI_CAR1b <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day))
d_rmI_CAR1c <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day|Rep))

#no random effects and CAR1 correlation
d_rm0_CAR1 <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1())
d_rm0_CAR1a <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~1|Rep))
#d_rm0_CAR1b <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day))
#doesn't work
d_rm0_CAR1c <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCAR1(form=~day|Rep))

#full model with random effects (slope and intercept) and CompSymm correlation
#none worked
#d_rm_CompSymm <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm())
#d_rm_CompSymmb <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day))
#d_rm_CompSymmc <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day|Rep))


#full model with random effects (intercept) and CompSymm correlation
d_rmI_CompSymm <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm())
d_rmI_CompSymmb <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day))
d_rmI_CompSymmc <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day|Rep))

#no random effects and CompSymm correlation
d_rm0_CompSymm <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm())
d_rm0_CompSymma <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~1|Rep))
d_rm0_CompSymmb <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day))
d_rm0_CompSymmc <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corCompSymm(form=~day|Rep))

#corSymm has convergence problems in all tried

#d_rm_Symm <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corSymm())
#d_rm_Symmb <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~day))
#d_rm_Symmc <- lme(gr ~ day*Treatment, random = ~ day|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~day|Rep))

#full model with random effects (intercept) and CompSymm correlation
#d_rmI_Symm <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corSymm())
#d_rmI_Symmb <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~day))
#d_rmI_Symmc <- lme(gr ~ day*Treatment, random = ~ 1|Rep, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~day|Rep))

#no random effects and CompSymm correlation
#d_rm0_Symm <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corSymm())
#d_rm0_Symma <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~1|Rep))
#d_rm0_Symmb <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~day))
#d_rm0_Symmc <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corSymm(form=~day|Rep))

#wouldn't finish running
#d_rm0_Symm <- gls(gr ~ day*Treatment, data=d_evo_rm, na.action=na.exclude, correlation=corSymm())

aic_comp_d <- AIC(d_rmI,d_rm0,d_rm_AR1,d_rmI_AR1,d_rmI_AR1b,d_rmI_AR1c,d_rm0_AR1,d_rm0_AR1a,d_rm0_AR1c,d_rm_ARMA,d_rmI_ARMA,d_rmI_ARMAb,d_rmI_ARMAc,d_rm0_ARMA,d_rm0_ARMAa,d_rm0_ARMAc,d_rm_CAR1,d_rmI_CAR1,d_rmI_CAR1b,d_rmI_CAR1c,d_rm0_CAR1,d_rm0_CAR1a,d_rm0_CAR1c,d_rmI_CompSymm,d_rmI_CompSymmb,d_rmI_CompSymmc,d_rm0_CompSymm,d_rm0_CompSymma,d_rm0_CompSymmb,d_rm0_CompSymmc)

which(aic_comp_d$AIC == min(aic_comp_d$AIC))
aic_comp_d[15,]
#            df       AIC
#d_rm0_ARMAa  7 -286.4654
which(aic_comp_d$AIC < -270)
aic_comp_d[c(8,10,11,15,22),]
#           df       AIC
#d_rm_ARMA   10 -280.4654
#d_rmI_ARMA   8 -284.4654
#d_rm0_ARMAa  7 -286.4654


anov_comp_d <- anova(d_rmI,d_rm0,d_rm_AR1,d_rmI_AR1,d_rmI_AR1b,d_rmI_AR1c,d_rm0_AR1,d_rm0_AR1a,d_rm0_AR1c,d_rm_ARMA,d_rmI_ARMA,d_rmI_ARMAb,d_rmI_ARMAc,d_rm0_ARMA,d_rm0_ARMAa,d_rm0_ARMAc,d_rm_CAR1,d_rmI_CAR1,d_rmI_CAR1b,d_rmI_CAR1c,d_rm0_CAR1,d_rm0_CAR1a,d_rm0_CAR1c,d_rmI_CompSymm,d_rmI_CompSymmb,d_rmI_CompSymmc,d_rm0_CompSymm,d_rm0_CompSymma,d_rm0_CompSymmb,d_rm0_CompSymmc,test=TRUE)

#save anova table output
save_as_docx(flextable(anov_comp_d),path="tables/tableSup1_d.docx")

which(anov_comp_d$logLik==max(anov_comp_d$logLik))
anov_comp_d[15,]
#           Model df       AIC       BIC   logLik Test L.Ratio p-value
#d_rm0_ARMAa    15  7 -286.4654 -262.7739 150.2327
which(anov_comp_d$logLik>145)
anov_comp_d[c(10,11,15),]
#all the same
anova(d_rm0_CAR1a,d_rm_ARMA, d_rmI_ARMA, d_rm0_ARMAa, d_rm0_AR1a)

anova(d_rmI,d_rm0)
#with only fixed effects sig worse
anova(d_rm_AR1,d_rm_ARMA,d_rm_CAR1)
#d_rm_ARMA best
anova(d_rmI,d_rmI_AR1,d_rmI_AR1b,d_rmI_AR1c,d_rmI_ARMA,d_rmI_ARMAb,d_rmI_ARMAc,d_rmI_CAR1,d_rmI_CAR1b,d_rmI_CAR1c,d_rmI_CompSymm,d_rmI_CompSymmb,d_rmI_CompSymmc)
#d_rmI_ARMA best
anova(d_rm0,d_rm0_AR1,d_rm0_AR1a,d_rm0_AR1c,d_rm0_ARMA,d_rm0_ARMAa,d_rm0_ARMAc,d_rm0_CAR1,d_rm0_CAR1a,d_rm0_CAR1c,d_rm0_CompSymm,d_rm0_CompSymma,d_rm0_CompSymmb,d_rm0_CompSymmc)
#d_rm0_ARMAa best
anova(d_rm_AR1,d_rmI_AR1,d_rmI_AR1b,d_rmI_AR1c,d_rm0_AR1,d_rm0_AR1a,d_rm0_AR1c)
#close between d_rm0_AR1a, d_rm_AR1, d_rmI_AR1
anova(d_rm_ARMA,d_rmI_ARMA,d_rmI_ARMAb,d_rmI_ARMAc,d_rm0_ARMA,d_rm0_ARMAa,d_rm0_ARMAc)
#close between d_rm_ARMA, d_rmI_ARMA, d_rm0_ARMAa
anova(d_rm_CAR1,d_rmI_CAR1,d_rmI_CAR1b,d_rmI_CAR1c,d_rm0_CAR1,d_rm0_CAR1a,d_rm0_CAR1c)
#close between most except d_rm0_CAR1
anova(d_rmI_CompSymm,d_rmI_CompSymmb,d_rmI_CompSymmc,d_rm0_CompSymm,d_rm0_CompSymma,d_rm0_CompSymmb,d_rm0_CompSymmc)
#many same/similar, including d_rm0_CompSymmc

#d_rm0_ARMAa has lowest AIC but d_rm_ARMA, d_rmI_ARMA as good.

#checking plots
plot(d_rm0_ARMAa)

plot(d_rm0, resid(., type = "normalized") ~ day | Rep, abline = 0)
plot(d_rmI, resid(., type = "normalized") ~ day | Rep, abline = 0)
plot(d_rm0_ARMAa, resid(., type = "normalized") ~ day | Rep, abline = 0)
#R6 27 still crazy
plot(d_rm0, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
plot(d_rm0_ARMAa, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
plot(d_rm_ARMA, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
plot(d_rmI, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
plot(d_rmI_ARMA, resid(., type = "normalized") ~ fitted(.) | Rep, abline = 0)
#these don't work for rm0
plot(ACF(d_rm_ARMA, resType = "normalized"), alpha=0.05, main="full mixed, ARMA corr")
plot(ACF(d_rmI, resType = "normalized"), alpha=0.05, main="mixed ~1|Rep, no corr")
plot(ACF(d_rmI_ARMA, resType = "normalized"), alpha=0.05, main="mixed ~1|Rep, ARMA corr")
#temporal autocorrelation worse????

#look at best models
summary(d_rm0_ARMAa)
summary(d_rm_ARMA)
summary(d_rmI_ARMA)
#VarCorr(b_rm0_ARMA)
confint(d_rm0_ARMAa, parm=c("day","Treatment27"))
intervals(d_rm0_ARMAa, which="all", level=0.95)
intervals(d_rm_ARMA, which="fixed", level=0.95)
intervals(d_rmI_ARMA, which="fixed", level=0.95)
anova(d_rm0_ARMAa)
anova(d_rm_ARMA)
anova(d_rmI_ARMA)
fixef(d_rm_ARMA)
fixef(d_rmI_ARMA)
ranef(d_rm_ARMA)
ranef(d_rmI_ARMA)

#plot - could also use , type="contrast"
visreg(d_rm0_ARMAa, "day", by="Treatment", gg=TRUE, overlay=TRUE)
visreg(d_rm0_ARMAa, "day", by="Treatment", overlay=TRUE, ylab="Growth rate", xlab="Day", gg=TRUE)+
  theme_bw() +
  scale_color_manual(values=c("blue","red")) +
  scale_fill_manual(values=c("blue","red")) +
  labs(title="ProD")






#######left off here - visualize and do d and omit those others and get estimates for each?###########



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








fit1.rm <- lm(gr~day, data=b_subt23.rm)








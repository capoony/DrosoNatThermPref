library(lme4)
library(afex)
library(emmeans)
library(ggplot2)
library(sjmisc)
require(gridExtra)
library(tidyverse)
library(readxl)
library(ggpubr)

########Basic Data exploration and data shaping ####
DATA<-read_excel("/Users/antonstrunov/Desktop/wild flies_thermal preference/thermal pref_qPCR results_29.03.23.xlsx", sheet="altogether")

### convert to factors
DATA$BioRep=as.factor(DATA$BioRep)
DATA$TechRep=as.factor(DATA$TechRep)
DATA$Run=as.factor(DATA$Run)

## calculate mean Ct across technical replicates and remove samples with SD>=0.5
means=DATA %>%
  group_by(Run,Line,BioRep,Gene,Sex,median_Tp,mean_Tp, WolbStrain) %>%
  dplyr::summarise(Mean=mean(Ct),SD=sd(Ct)) %>%
  filter(SD<0.5)%>%
  select(!SD)
means

## calculate mean Ct across technical replicates and retain samples with SD>=0.3 to get a list of samples to be repeated
means.SD=DATA %>%
  group_by(Run,Line,BioRep,Gene,Sex) %>%
  summarise(Mean=mean(Ct),SD=sd(Ct)) %>%
  filter(SD>=0.5)
means.SD

#write.table(means.SD,file="/Users/Anton/Keepitcool!/qPCR_competition assay fecundity/qPCR_2B_repeated.txt",quote = F,row.names = F)

## calculate delta CT
means.spread=na.omit(spread(means, Gene, Mean))
means.spread$delta=2^(-(means.spread$wd-means.spread$rpl))

## plot Boxplots (direct from runs)
Run.data=subset(means.spread,means.spread$Run!="stock") 

Figure=ggplot(Run.data, aes(x=WolbStrain, y=delta, col=Run))
Figure+geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8),
                   size = 2)+
  theme_bw()+
  stat_summary(
    aes(color = Run),
    fun.data="mean_sdl",  fun.args = list(mult=1),
    geom = "crossbar",  size = 0.4,
    position = position_dodge(0.8))+
  theme(text = element_text(size=15))+
  xlab("")+scale_y_continuous(name="delta Ct", breaks=seq(0,8,2))

## plot Boxplots (stocks)
Stock.data=subset(means.spread,means.spread$Run=="stock") 

Figure=ggplot(Stock.data, aes(x=WolbStrain, y=delta))
Figure+geom_jitter(aes(), position = position_jitter(0.2),
                   size = 4)+
  theme_bw()+
  stat_summary(
    aes(),
    fun.data="mean_sdl",  fun.args = list(mult=1),
    geom = "crossbar",  size = 0.4,
    position = position_dodge(0.8))+
  theme(text = element_text(size=15))+
  xlab("")+scale_y_continuous(name="delta Ct", breaks=seq(0,8,2))

## plot Regression analysis  

Figure2=ggplot(Run.data, aes(x=median_Tp, y=delta, color=WolbStrain))+
  geom_point()+geom_smooth(method=lm)+theme_classic()+
  facet_grid(WolbStrain~.)
Figure2

Figure3=ggplot(Run.data, aes(x=mean_Tp, y=delta, color=WolbStrain))+
  geom_point()+geom_smooth(method=lm)+theme_classic()+
  facet_grid(WolbStrain~.)
Figure3

### Statistics

### Direct from runs

### How Wolbachia strain affects titer in an individual? Directly from runs.

options(contrasts = c("contr.sum", "contr.poly"))

Titer.full=lmer(delta ~ WolbStrain + Sex + (1|Run), data=Run.data)
Titer.null.strain=lmer(delta ~ Sex + (1|Run), data=Run.data)
Titer.null.sex=lmer(delta ~ WolbStrain + (1|Run), data=Run.data)


summary(Titer.full)

anova(Titer.full, Titer.null.strain, type=3, test.statistic = "F")
anova(Titer.full, Titer.null.sex, type=3, test.statistic = "F")


emmeans(Titer.full, pairwise ~ WolbStrain)

### How Wolbachia strain affects titer in an individual? In stocks.

options(contrasts = c("contr.sum", "contr.poly"))

Titer.full2=lmer(delta ~ WolbStrain + Sex + (1|Line), data=Stock.data)
Titer.null2.strain=lmer(delta ~ Sex + (1|Line), data=Stock.data)
Titer.null2.sex=lmer(delta ~ WolbStrain + (1|Line), data=Stock.data)


summary(Titer.full2)

anova(Titer.full2, Titer.null2.strain, type=3, test.statistic = "F")
anova(Titer.full2, Titer.null2.sex, type=3, test.statistic = "F")


emmeans(Titer.full2, pairwise ~ WolbStrain)

### Is there a correlation between Wolbachia titer and Tp?

ggscatter(means.spread, x = "mean_Tp", y = "delta", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Temperature (°C)", ylab = "Wolbachia titer")

Figure=ggscatter(Run.data, x = "mean_Tp", y = "delta", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Temperature (°C)", ylab = "Wolbachia titer")+facet_grid(WolbStrain~.)+theme(text = element_text(size=15))
Figure

## Is data normally distributed?

ggplot(means.spread, aes(x=delta)) + geom_histogram(binwidth=.5)
ggplot(means.spread, aes(x=mean_Tp)) + geom_histogram(binwidth=.5)

shapiro.test(means.spread$mean_Tp)
shapiro.test(means.spread$delta)

ggqqplot(means.spread$mean_Tp, ylab = "mean_Tp")
ggqqplot(means.spread$delta, ylab = "delta")

res <- cor.test(means.spread$mean_Tp, means.spread$delta, 
                method = "pearson")
res

res2 <- cor.test(means.spread$mean_Tp, means.spread$delta,  method="kendall")
res2

res3 <-cor.test(means.spread$mean_Tp, means.spread$delta,  method = "spearman", exact = FALSE)
res3

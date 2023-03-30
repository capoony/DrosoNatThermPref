library(tidyverse)
library(readxl)
library(readr)
library(lme4)
library(afex)
library(emmeans)
library(multcomp)
library(caret)
library(leaps)
library(car)
library(sjPlot)

DATA<-read_excel("/Users/antonstrunov/Desktop/wild flies_thermal preference/Tp_data_transformed.xlsx")
#DATA$replica=as.factor(DATA$replica)
#DATA$exp_date=as.factor(DATA$exp_date)
#DATA$country=as.numeric(DATA$country)
#DATA$Infection=as.numeric(DATA$Infection)
#DATA$Wolbachia_strain=as.numeric(DATA$Wolbachia_strain)
#DATA$antibiotic_treatment=as.numeric(DATA$antibiotic_treatment)
#DATA$country=as.numeric(DATA$country)
#DATA$Run=as.factor(DATA$Run)
#DATA$slab_lane=as.factor(DATA$slab_lane)
#DATA$ambient_temp=as.factor(DATA$ambient_temp)
#DATA$humidity=as.factor(DATA$humidity)
#DATA$RepTime=paste0(DATA$exp_date,DATA$time)


cat("**** Summary Table ****\n")

means=DATA %>%
  group_by(Wolbachia_strain) %>%
  dplyr::summarise(Mean=mean(TempEst),SD=sd(TempEst), SE=SD/sqrt(n()), Median=median(TempEst))
means

ggplot(DATA,aes(x=Wolbachia_strain,y=TempEst,col=Wolbachia_strain)) + geom_jitter() + geom_violin(alpha=0.2)+theme_classic()+
  theme(text = element_text(size=10))+
  xlab("Infection type")+scale_y_continuous(name="Temperature (°C)", breaks=seq(10,36,1))


#scale_colour_manual(values=c("#999999","#0072B2","#D55E00", "#009E73"), labels=labels)


### Finding the best model to fit the data

set.seed(123)
train.control <- trainControl(method = "cv", number = 10)
step.model <- train(TempEst~Wolbachia_strain+country+Infection+antibiotic_treatment+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA,
                    method = "leapSeq", 
                    tuneGrid = data.frame(nvmax = 1:8),
                    trControl = train.control
)
step.model$results
step.model$bestTune
summary(step.model$finalModel)


cat("\n**** Linear mixed model ****\n")

options(contrasts = c("contr.sum","contr.poly"))

### Wolbachia strain (wMel+/wMel- and other pairs)

LMM1=lmer(TempEst~Wolbachia_strain+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM1.null.infection=lmer(TempEst~(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
summary(LMM1)

anova(LMM1, LMM1.null.infection, type=3, test.statistic = "F")
emmeans(LMM1, pairwise ~ Wolbachia_strain)

### Wolbachia infection (+/-)

LMM2=lmer(TempEst~Infection+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM2.null.infection=lmer(TempEst~(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
summary(LMM2)

anova(LMM2, LMM2.null.infection, type=3, test.statistic = "F")

### Wolbachia infection with all fixed effects

LMM3=lmer(TempEst~Infection+country+antibiotic_treatment+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM3.null.infection=lmer(TempEst~country+antibiotic_treatment+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM3.null.antibiotic=lmer(TempEst~country+Infection+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM3.null.country=lmer(TempEst~Infection+antibiotic_treatment+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)

summary(LMM3)

anova(LMM3, LMM3.null.infection, type=3, test.statistic = "F")
anova(LMM3, LMM3.null.antibiotic, type=3, test.statistic = "F")
anova(LMM3, LMM3.null.country, type=3, test.statistic = "F")

### Wolbachia status (wMel, wMelCS, uninfected, rif_treated)

LMM4=lmer(TempEst~Wolbachia_status+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM4.null.status=lmer(TempEst~(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
summary(LMM4)

anova(LMM4, LMM4.null.status, type=3, test.statistic = "F")
emmeans(LMM4, pairwise ~ Wolbachia_status)

### Wolbachia status by country (wMel Portugal vs wMel Finland)

LMM5=lmer(TempEst~Wolbachia_status*country+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM5.null.status=lmer(TempEst~country+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM5.null.country=lmer(TempEst~Wolbachia_status+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM5.null.interaction=lmer(TempEst~country+Wolbachia_status+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)

summary(LMM5)

anova(LMM5, LMM5.null.status, type=3, test.statistic = "F")
anova(LMM5, LMM5.null.country, type=3, test.statistic = "F")
anova(LMM5, LMM5.null.interaction, type=3, test.statistic = "F")

emmeans(LMM5, pairwise ~ Wolbachia_status*country)

### Country (Portugal vs Finland)

No_antibiotic=subset(DATA,DATA$antibiotic_treatment!="yes")
  
LMM6=lmer(TempEst~country+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=No_antibiotic)
LMM6.null.country=lmer(TempEst~(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=No_antibiotic)
summary(LMM6)

anova(LMM6, LMM6.null.country, type=3, test.statistic = "F")

### Antibiotic treatment

LMM7=lmer(TempEst~antibiotic_treatment+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM7.null.antibiotic=lmer(TempEst~(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
summary(LMM7)

anova(LMM7, LMM7.null.antibiotic, type=3, test.statistic = "F")

### Antibiotic treatment by infection

LMM8=lmer(TempEst~antibiotic_treatment*Infection+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM8.null.antibiotic=lmer(TempEst~Infection+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM8.null.infection=lmer(TempEst~antibiotic_treatment+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)
LMM8.null.interaction=lmer(TempEst~antibiotic_treatment+Infection+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=DATA)

summary(LMM8)

anova(LMM8, LMM8.null.antibiotic, type=3, test.statistic = "F")
anova(LMM8, LMM8.null.infection, type=3, test.statistic = "F")
anova(LMM8, LMM8.null.interaction, type=3, test.statistic = "F")

emmeans(LMM8, pairwise ~ antibiotic_treatment*Infection)

### Country by strain

LMM9=lmer(TempEst~country*Wolbachia_strain+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=No_antibiotic)
LMM9.null.country=lmer(TempEst~Wolbachia_strain+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=No_antibiotic)
LMM9.null.strain=lmer(TempEst~country+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=No_antibiotic)
LMM9.null.interaction=lmer(TempEst~country+Wolbachia_strain+(1|Run)+(1|slab_lane)+(1|replica)+(1|ambient_temp)+(1|humidity), data=No_antibiotic)

summary(LMM9)

anova(LMM9, LMM9.null.country, type=3, test.statistic = "F")
anova(LMM9, LMM9.null.strain, type=3, test.statistic = "F")
anova(LMM9, LMM9.null.interaction, type=3, test.statistic = "F")

emmeans(LMM9, pairwise ~ Wolbachia_strain*country)
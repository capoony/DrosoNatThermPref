library(lme4)
library(afex)
library(emmeans)
library(ggplot2)
library(sjmisc)
library(plyr)
require(gridExtra)
library(grid)
library(tidyverse)
library(readxl)
library(car)
library(cowplot)
library(knitr)
library(emmeans)

### set working directory and make new output folder

# setwd("/media/inter/mkapun/projects/")
setwd("/Users/martinkapun/Documents/GitHub/")

### load Data
DATA <- read_excel("DrosoNatThermPref/data/thermal pref_qPCR results.xlsx",
    sheet = "altogether",
    na = "NA"
)

### revalue and convert to factors
DATA$BioRep <- as.factor(DATA$BioRep)
DATA$WolbType[DATA$WolbType == "cs"] <- "wMelCS"
DATA$WolbType[DATA$WolbType == "mel"] <- "wMel"
DATA$Tissue[DATA$Tissue == "h"] <- "head"
DATA$Tissue[DATA$Tissue == "o"] <- "ovaries"
DATA$LineID <- as.factor(paste(DATA$TempComp,
    DATA$Line,
    sep = "_"
))
DATA$TechRep <- as.factor(DATA$TechRep)
DATA$TempFec <- as.factor(DATA$TempFec)
DATA$TempComp <- as.factor(DATA$TempComp)
DATA$WolbType <- as.factor(DATA$WolbType)

DATA$median_Tp <- as.factor(DATA$median_Tp)

means <- DATA %>%
    group_by(WolbType, Line, Tp_id, BioRep, Sex, Gene, median_Tp) %>%
    summarise(Mean = mean(Ct), SD = sd(Ct), N = n()) %>%
    select(!SD)
means

means.spread <- na.omit(spread(means, Gene, Mean))
means.spread$delta <- 2^(-(means.spread$wd - means.spread$rpl))
# means.spread <- subset(means.spread, means.spread$fecundity != 0)

Figure <- ggplot(means.spread, aes(x = median_Tp, y = delta)) +
    facet_grid(Line ~ .) +
    theme_bw() +
    geom_point() +
    stat_summary(fun.data = mean_cl_normal) +
    geom_smooth(method = "lm", formula = y ~ x) +
    theme(axis.title.y = element_text(size = 22, angle = 90)) +
    theme(axis.title.x = element_text(size = 22, angle = 00)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.text = element_text(size = 20)) +
    theme(legend.title = element_text(size = 20)) +
    theme(strip.text = element_text(size = 20))
Figure

ggsave("DrosoNatThermPref/analyses/TpDelta.pdf",
    Figure,
    width = 5,
    height = 8
)

ggsave("DrosoNatThermPref/analyses/TpDelta.png",
    Figure,
    width = 5,
    height = 8
)

sink("DrosoNatThermPref/analyses/TpDelta.stats")
print("##### Delta vs. Tp #####")

options(contrasts = c("contr.sum", "contr.poly"))
res <- lmer(median_Tp ~ Line * delta * Sex + (1 | BioRep),
    data = means.spread
)
anova(res, type = 3)
sink()

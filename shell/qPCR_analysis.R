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
DATA <- read_excel("DrosoNatThermPref/data/thermal pref_qPCR results_13.03.23.xlsx",
    sheet = "altogether"
)

means <- DATA %>%
    group_by(WolbType, Line, Tp_id, BioRep, Sex, Gene, median_Tp) %>%
    summarise(Mean = mean(Ct), SD = sd(Ct), N = n()) %>%
    select(!SD)
means

means.spread <- na.omit(spread(means, Gene, Mean))
means.spread$delta <- 2^(-(means.spread$wd - means.spread$rpl))
# means.spread <- subset(means.spread, !means.spread$Line %in% c("ak7", "ak9"))
# means.spread$Line[means.spread$Line == "re1"] <- "wMel+1"
# means.spread$Line[means.spread$Line == "re10"] <- "wMelCS+2"

Figure <- ggplot(means.spread, aes(x = Line, y = delta)) +
    theme_bw() +
    geom_boxplot() +
    theme(axis.title.y = element_text(size = 22, angle = 90)) +
    theme(axis.title.x = element_text(size = 22, angle = 00)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.text = element_text(size = 20)) +
    theme(legend.title = element_text(size = 20)) +
    theme(strip.text = element_text(size = 20)) +
    ylab(expression(2^("-" ~ Delta ~ italic("Ct")))) +
    xlab("Line ID")
Figure

ggsave("DrosoNatThermPref/analyses/DeltaFull.pdf",
    Figure,
    width = 8,
    height = 3
)

ggsave("DrosoNatThermPref/analyses/DeltaFull.png",
    Figure,
    width = 8,
    height = 3
)

sink("DrosoNatThermPref/analyses/TpDelta.stats")
print("##### Delta vs. Tp #####")

options(contrasts = c("contr.sum", "contr.poly"))
res <- lmer(median_Tp ~ Line * delta * Sex + (1 | BioRep),
    data = means.spread
)
anova(res, type = 3)
sink()

library(lme4)
library(afex)
library(emmeans)
library(ggplot2)
library(car)
require(gridExtra)
library(tidyverse)
library(readxl)
library(ggpubr)

setwd("D:/GitHub/DrosoNatThermPref/analyses")

######## Basic Data exploration and data shaping ####
DATA <- read_excel("../data/thermal pref_qPCR results_29.03.23_MK.xlsx",
    sheet = "altogether"
)


### convert to factors
DATA$BioRep <- as.factor(DATA$BioRep)
DATA$TechRep <- as.factor(DATA$TechRep)
DATA$Run <- as.factor(DATA$Run)
levels(factor(DATA$WolbStrain))

color1 <- c("firebrick3", "blue3", "blue3", "blue3")
color2 <- c("firebrick3", "firebrick3", "blue3", "blue3", "blue3", "blue3")
DATA$WolbStrain <- factor(DATA$WolbStrain, levels = c("wMelCS_Portugal_1", "wMelCS_Portugal_2", "wMel_Portugal_1", "wMel_Portugal_2", "wMel_Finland_1", "wMel_Finland_2"))


## calculate mean Ct across technical replicates and remove samples with SD>=0.5
means <- DATA %>%
    group_by(Run, Line, BioRep, Gene, Sex, median_Tp, mean_Tp, WolbStrain) %>%
    dplyr::summarise(Mean = mean(Ct), SD = sd(Ct)) %>%
    select(!SD)
means

## calculate mean Ct across technical replicates and retain samples with SD>=0.3 to get a list of samples to be repeated
means.SD <- DATA %>%
    group_by(Run, Line, BioRep, Gene, Sex) %>%
    summarise(Mean = mean(Ct), SD = sd(Ct)) %>%
    filter(SD >= 0.5)
means.SD

# write.table(means.SD,file="/Users/Anton/Keepitcool!/qPCR_competition assay fecundity/qPCR_2B_repeated.txt",quote = F,row.names = F)

## calculate delta CT
means.spread <- na.omit(spread(means, Gene, Mean))
means.spread$delta <- 2^(-(means.spread$wd - means.spread$rpl))

## plot RegressionPlots (direct from runs)
Run.data <- subset(means.spread, means.spread$Run != "stock")

qPCR.reg <- ggplot(Run.data, aes(x = mean_Tp, y = delta, color = WolbStrain, fill = WolbStrain)) +
    geom_point(show.legend = FALSE) +
    geom_smooth(
        method = lm,
        show.legend = FALSE
    ) +
    theme_bw() +
    ylab(expression(2^("-" ~ Delta ~ italic("Ct")))) +
    xlab("Mean Thermal Preference (°C)") +
    facet_grid(. ~ WolbStrain) +
    scale_color_manual(values = color1) +
    scale_fill_manual(values = color1) +
    stat_cor(label.y = 7.8, col = "black")


ggsave(
    "MK/qPCR_regressions.pdf",
    qPCR.reg,
    width = 8,
    height = 5
)

ggsave(
    "MK/qPCR_regressions.png",
    qPCR.reg,
    width = 8,
    height = 5
)

Run.data$IDRep <- paste0(Run.data$BioRep, Run.data$Sex)
res <- glm(delta ~ WolbStrain * mean_Tp,
    data = Run.data
)
Anova(res, type = 3)

Stock.data <- subset(means.spread, means.spread$Run == "stock")


qPCR.bp <- ggplot(Stock.data, aes(x = WolbStrain, y = delta, color = WolbStrain, fill = WolbStrain)) +
    geom_jitter(
        show.legend = FALSE,
        alpha = 0.2
    ) +
    geom_boxplot(
        show.legend = FALSE,
        alpha = 0.4
    ) +
    theme_bw() +
    ylab(expression(2^("-" ~ Delta ~ italic("Ct")))) +
    xlab("") +
    scale_color_manual(values = color2) +
    scale_fill_manual(values = color2) +
    stat_cor(label.y = 7.8, col = "black")
qPCR.bp

ggsave(
    "MK/qPCR_boxplot.pdf",
    qPCR.bp,
    width = 8,
    height = 5
)

ggsave(
    "MK/qPCR_boxplot.png",
    qPCR.bp,
    width = 8,
    height = 5
)

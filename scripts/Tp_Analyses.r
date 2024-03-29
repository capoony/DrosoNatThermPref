library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(knitr)
library(car)

setwd("D:/GitHub/DrosoNatThermPref")

dir.create("analyses/Tp/")

DATA <- read_excel("data/Strunov_DrosNatThermPref_2023_Tp.xlsx")
DATA$Run <- factor(DATA$Run)
DATA$slab_lane <- factor(DATA$slab_lane)
DATA$replica <- factor(DATA$replica)
DATA$ID2 <- paste(DATA$Wolbachia_status, DATA$country, DATA$BioRep, sep = "_")
DATA$RepID <- paste0(DATA$Wolbachia_strain, DATA$replica)

## reorder Strains
DATA$Wolbachia_strain <- factor(DATA$Wolbachia_strain, levels = c("Uninf_Po_1_-", "Uninf-AB_Po_1_-", "Uninf_Po_2_-", "Uninf-AB_Po_2_-", "wMel_Po_1_+", "wMel_Po_1_-", "wMel_Po_2_+", "wMel_Po_2_-", "wMel_Fi_1_+", "wMel_Fi_1_-", "wMel_Fi_2_+", "wMel_Fi_2_-", "wMelCS_Po_1_+", "wMelCS_Po_1_-", "wMelCS_Po_2_+", "wMelCS_Po_2_-"))

## Define Strain Colors
color <- c("darkgreen", "lightgreen", "darkgreen", "lightgreen", "blue3", "lightsteelblue", "blue3", "lightsteelblue", "blue3", "lightsteelblue", "blue3", "lightsteelblue", "firebrick3", "lightsalmon", "firebrick3", "lightsalmon")

DATA2 <- DATA
DATA2$Wolbachia_status[DATA2$Wolbachia_status == "Uninf-AB"] <- "Uninf"
DATA2$ID2 <- paste(DATA2$Wolbachia_status, DATA2$country, DATA2$BioRep, sep = "_")


means <- DATA2 %>%
    group_by(ID2) %>%
    dplyr::summarise(Mean = mean(TempEst), SD = sd(TempEst), SE = SD / sqrt(n()), Median = median(TempEst))


write.table(
    file = "analyses/Tp/TP_means.txt",
    means,
    quote = F,
    row.names = F
)

## now plot averages
TP.plot <- ggplot(DATA, aes(x = Wolbachia_strain, y = TempEst, col = Wolbachia_strain, fill = Wolbachia_strain)) +
    facet_grid(. ~ country, scales = "free_x", space = "free") +
    geom_jitter(alpha = 0.1, show.legend = FALSE) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    theme_bw() +
    theme(text = element_text(size = 10)) +
    xlab("SampleID (Wolbachiastrain,Origin,Infectionstatus)") +
    scale_y_continuous(name = "Temperature (°C)", breaks = seq(10, 36, 1)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    "analyses/Tp/TP_plot.pdf",
    TP.plot,
    width = 8,
    height = 6
)

ggsave(
    "analyses/Tp/TP_plot.png",
    TP.plot,
    width = 8,
    height = 6
)

## now plot by replicates
TP.replicas <- ggplot(DATA, aes(x = as.factor(replica), y = TempEst, col = Wolbachia_strain, fill = Wolbachia_strain)) +
    facet_grid(. ~ Wolbachia_strain, scales = "free_x", space = "free") +
    geom_jitter(alpha = 0.1, show.legend = FALSE) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    theme_bw() +
    theme(text = element_text(size = 10)) +
    xlab("SampleID (Wolbachiastrain,Origin,Infectionstatus)") +
    scale_y_continuous(name = "Temperature (°C)", breaks = seq(10, 36, 1)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    "analyses/Tp/TP_replicas.pdf",
    TP.replicas,
    width = 20,
    height = 6
)

ggsave(
    "analyses/Tp/TP_replicas.png",
    TP.replicas,
    width = 20,
    height = 6
)

## now plot by replicates
TP.humidity <- ggplot(DATA, aes(x = humidity, y = TempEst, col = Wolbachia_strain, fill = Wolbachia_strain)) +
    geom_jitter(
        show.legend = FALSE,
        alpha = 0.2
    ) +
    geom_smooth(
        method = lm,
        show.legend = FALSE
    ) +
    theme_bw() +
    theme(text = element_text(size = 10)) +
    xlab("Humidity (%)") +
    scale_y_continuous(name = "Temperature (°C)", breaks = seq(10, 36, 1)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
TP.humidity

ggsave(
    "analyses/Tp/TP_humidity.pdf",
    TP.humidity,
    width = 6,
    height = 4
)

ggsave(
    "analyses/Tp/TP_humidity.png",
    TP.humidity,
    width = 6,
    height = 4
)

sink("analyses/Tp/Stats.txt")
### Finding the best model to fit the data

## adjust contrasts to fit Type - III ANOVAs in R, see here https://www.r-bloggers.com/2011/03/anova-%E2%80%93-type-iiiiii-ss-explained/
## essentially this uses sum contrasts to compare each group against grand mean.
options(contrasts = c("contr.sum", "contr.poly"))

#### First test, is there a difference between lines and AB treatment

RES.line <- lmer(TempEst ~ line * antibiotic_treatment + ambient_temp + humidity + (1 | Replicate) + (1 | slab_lane / Run),
    data = DATA
)

cat("#### Test by Lines ####")
print(kable(Anova(RES.line,
    type = 3
)))

cat("\n#### PostHoc Test ####")
print(kable(pairs(emmeans(RES.line, ~ antibiotic_treatment | line))))

#### Then test, if there is a difference between Wolbachia status and AB treatment
RES.Wolb <- lmer(TempEst ~ Wolbachia_status * antibiotic_treatment + country * antibiotic_treatment + ambient_temp + humidity + (1 | line) + (1 | Replicate) + (1 | slab_lane / Run) + (1 | BioRep / country),
    data = DATA2
)

cat("\n#### Test by Wolbachia Status ####")

print(kable(Anova(RES.Wolb,
    type = 3
)))

RES.Rep <- lmer(TempEst ~ Replicate + ambient_temp + humidity + (1 | line) + (1 | slab_lane / Run) + (1 | BioRep / country),
    data = DATA2
)

cat("\n#### Test by Wolbachia Status ####")

print(kable(Anova(RES.Rep)))

cat("\n#### PostHoc Test ####")
print(kable(pairs(emmeans(RES.Wolb, ~ antibiotic_treatment | Wolbachia_status))))

sink()


means <- DATA2 %>%
    group_by(Wolbachia_status, antibiotic_treatment) %>%
    dplyr::summarise(Mean = mean(TempEst), SD = sd(TempEst), SE = SD / sqrt(n()), Median = median(TempEst))
means

TP <- ggplot(means, aes(x = Wolbachia_status, y = Mean, col = antibiotic_treatment, group = antibiotic_treatment)) +
    geom_line(linewidth = 2, position = position_dodge(width = 0.75)) +
    geom_point(position = position_dodge(width = 0.75)) +
    geom_pointrange(aes(ymin = Mean - SE, ymax = Mean + SE), linewidth = 2, position = position_dodge(width = 0.75)) +
    theme_bw() +
    theme(text = element_text(size = 10)) +
    xlab("SampleID (Wolbachiastrain,Origin,Infectionstatus)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
TP

ggsave(
    "analyses/Tp/TP_interaction.pdf",
    TP,
    width = 8,
    height = 6
)

ggsave(
    "analyses/Tp/TP_interaction.png",
    TP,
    width = 8,
    height = 6
)


## what about replicates?? Does it matter to remove clustered samples?

means2 <- DATA2 %>%
    group_by(RepID) %>%
    dplyr::summarise(
        Mean = mean(TempEst),
        SD = sd(TempEst),
        SE = SD / sqrt(n()),
        Median = median(TempEst),
        N = n(),
        N90 = sum(TempEst < min(TempEst) + 2)
    ) %>%
    ## identify samples with more than 80% within 2°C from minimum where the median is < 20°C
    filter(Median < 20 & N90 / N > 0.80)

## remove those samples
DATA3 <- DATA2 %>%
    filter(!RepID %in% means2$RepID)

RES.line <- lmer(TempEst ~ line * antibiotic_treatment + ambient_temp + humidity + (1 | Replicate) + (1 | slab_lane / Run),
    data = DATA3
)

sink("analyses/Tp/Stats_reduced.txt")
cat("#### Test by Lines ####")
print(kable(Anova(RES.line,
    type = 3
)))

cat("\n#### PostHoc Test ####")
print(kable(pairs(emmeans(RES.line, ~ antibiotic_treatment | line))))

cat("\n#### Test by Wolbachia Status ####")

#### Then test, if there is a difference between Wolbachia status and AB treatment
RES.Wolb <- lmer(TempEst ~ Wolbachia_status * antibiotic_treatment + country * antibiotic_treatment + ambient_temp + humidity + (1 | line) + (1 | Replicate) + (1 | slab_lane / Run) + (1 | BioRep / country),
    data = DATA3
)

cat("\n#### Test by Wolbachia Status ####")

print(kable(Anova(RES.Wolb,
    type = 3
)))
cat("\n#### PostHoc Test ####")
print(kable(pairs(emmeans(RES.Wolb, ~ antibiotic_treatment | Wolbachia_status))))

sink()

TP.plot <- ggplot(DATA3, aes(x = Wolbachia_strain, y = TempEst, col = Wolbachia_strain, fill = Wolbachia_strain)) +
    facet_grid(. ~ country, scales = "free_x", space = "free") +
    geom_jitter(alpha = 0.1, show.legend = FALSE) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    theme_bw() +
    theme(text = element_text(size = 10)) +
    xlab("SampleID (Wolbachiastrain,Origin,Infectionstatus)") +
    scale_y_continuous(name = "Temperature (°C)", breaks = seq(10, 36, 1)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    "analyses/Tp/TP_plot_reduced.pdf",
    TP.plot,
    width = 8,
    height = 6
)

ggsave(
    "analyses/Tp/TP_plot_reduced.png",
    TP.plot,
    width = 8,
    height = 6
)

TP.replicas <- ggplot(DATA3, aes(x = as.factor(replica), y = TempEst, col = Wolbachia_strain, fill = Wolbachia_strain)) +
    facet_grid(. ~ Wolbachia_strain, scales = "free_x", space = "free") +
    geom_jitter(alpha = 0.1, show.legend = FALSE) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    theme_bw() +
    theme(text = element_text(size = 10)) +
    xlab("SampleID (Wolbachiastrain,Origin,Infectionstatus)") +
    scale_y_continuous(name = "Temperature (°C)", breaks = seq(10, 36, 1)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    "analyses/Tp/TP_replicas_reduced.pdf",
    TP.replicas,
    width = 20,
    height = 6
)

ggsave(
    "analyses/Tp/TP_replicas_reduced.png",
    TP.replicas,
    width = 20,
    height = 6
)

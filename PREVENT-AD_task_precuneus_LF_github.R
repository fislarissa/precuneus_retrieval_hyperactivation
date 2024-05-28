##########################################################################################################################

#Analysis of PREVENT-AD data - baseline fMRI and follow-up PET with APOE genotype & cognition
#Larissa Fischer & Eóin N. Molloy - Multimodal Neuroimaging Lab, DZNE Magdeburg

#########################################################################################################################

#Import packages:
library(readxl)
library(writexl)
library(dplyr)
library(lme4)
library(lmtest)
library(sjPlot)
library(car)
library(psych)
library(lmerTest)
library(lubridate)
library(MuMIn)
library(sjstats)
library(ggpubr)
library(ggplot2)
library(plotly)
library(devtools)
library(Lahman)
library(ppcor)
library(emmeans)
library(flextable)

#########################################################################################################################
#set working directory:
setwd("/Users/your/path")

#load data
load("/.../base.RData")
load("/.../prec.RData")
load("/.../rbans.RData")
load("/.../corrhr.RData")

#########################################################################################################################
#overview cohort
describe(base$Age_Baseline_Years) 
describe(base$Education_Baseline_Years)
summary(base$Sex) #1 = male
summary(base$RBANS_delayed_memory_index_score_Baseline) 

#overview time-points with rbans
counts <- table(rbans$Time)
count_df <- as.data.frame(counts)
colnames(count_df) <- c("time-point", "count")
count_df

#overview time-points with > 1 fMRI
counts <- table(prec$Time)
count_df <- as.data.frame(counts)
colnames(count_df) <- c("time-point", "count")
count_df

#########################################################################################################################
#subsets and theme
base_m <- subset(base, Sex == 1)
base_f <- subset(base, Sex == 0)
base_carrier <- subset(base, APOE4_Group == "B")
base_non_carrier <- subset(base, APOE4_Group == "A")

prec_m <- subset(prec, Sex == 1)
prec_f <- subset(prec, Sex == 0)

corrhr_m <- subset(corrhr, Sex == 1)
corrhr_f <- subset(corrhr, Sex == 0)

PREC_theme <- function() {
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
}

#########################################################################################################################
#Research question i)
#########################################################################################################################
#Assessment of baseline and longitudinal precuneus activation and APOE genotype
#effects on baseline activation?
baseprec <- lm(fMRI_Precuneus_Baseline ~ APOE4_Group+Age_Baseline_Months + Education_Baseline_Years + Sex + Prec_GM_Volume, data=base)
summary(baseprec) 
shapiro.test(rstandard(baseprec)) #normally distributed
bptest(baseprec) #homoscedasticity
vif(baseprec)   #no multicollinearity (a vif score over 5 is a problem)
tab_model(baseprec, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "baseprec.doc")
AIC(baseprec)
describeBy(base$fMRI_Precuneus_Baseline, base$Sex)

#effects of change in precuneus activation over time?
#first only Time, then add Time*APOE4_Group and Time*Sex to model
PREC <- lmer(PREC_Bilat ~ Time*Sex + Time*APOE4_Group +Age + Education + PREC_GMV + (1 |Subject), data=prec, REML = F)
summary(PREC)
confint(PREC)
AIC(PREC)
tab_model(PREC, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECwobistdu.doc")
#males
PREC_m <- lmer(PREC_Bilat ~ Time +Age + Education + PREC_GMV + (1 |Subject), data=prec_m, REML = F)
summary(PREC_m)
confint(PREC_m)
AIC(PREC_m)
tab_model(PREC_m, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC.doc")
#females
PREC_f <- lmer(PREC_Bilat ~ Time +Age + Education + PREC_GMV + (1 |Subject), data=prec_f, REML = F)
summary(PREC_f)
confint(PREC_f)
AIC(PREC_f)
tab_model(PREC_f, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC.doc")

#FIGURE APOE
plot_PREC_APOE <- ggplot(prec, aes(as.factor(Time), PREC_Bilat, colour = APOE4_Group)) + 
  geom_jitter(width = 0.25) +
  geom_smooth(aes(as.numeric(Time), PREC_Bilat), method = lm) +
  xlab ("Session") + ylab ("Precuneus BOLD-Signal") +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = expression(italic("APOE4") * " group")) +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme()
plot_PREC_APOE

#SUPP FIGURE SEX
plot_PREC_sex <- ggplot(prec, aes(as.factor(Time), PREC_Bilat, colour = Sex)) + 
  geom_jitter(width = 0.25) +
  geom_smooth(aes(as.numeric(Time), PREC_Bilat), method = lm) +
  xlab ("Session") + ylab ("Precuneus BOLD-Signal") +
  scale_color_manual(labels = c("Female", "Male"), values = c("#FFA500", "#0000FF")) + labs(color = "Sex") +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme()
plot_PREC_sex

#########################################################################################################################
### Research question ii)
#########################################################################################################################
### Assessment of baseline and longitudinal precuneus activation and AD pathology

#transform AD pathology scale
#amyloid distribution
ggdensity(base$Amyloid_WB_BoxCox)
ggqqplot(base$Amyloid_WB_BoxCox)
shapiro.test(base$Amyloid_WB_BoxCox) #amyloid not normally distributed, box-cox best approximation (transformation as for tau)
#tau distribution
ggdensity(base$Entorhinal_Tau_PET_LR_MEAN)
ggqqplot(base$Entorhinal_Tau_PET_LR_MEAN)
shapiro.test(base$Entorhinal_Tau_PET_LR_MEAN) #tau not normally distributed
#Box-Cox correction:
#compute a linear model and pass it to the boxcox function:
mod <- boxcox(lm(base$Entorhinal_Tau_PET_LR_MEAN ~ 1))
#calculate lambda
lambda <- mod$x[which.max(mod$y)]
lambda
#apply lambda to data
base$Entorhinal_Tau_PET_LR_MEAN_BoxCox <- (base$Entorhinal_Tau_PET_LR_MEAN ^ lambda - 1) / lambda
shapiro.test(base$Entorhinal_Tau_PET_LR_MEAN_BoxCox) #tau normally distributed


### effects of baseline activation?
#BL prec on amyloid:
PRECvsAMY <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY)
shapiro.test(rstandard(PRECvsAMY)) #not normally distributed
bptest(PRECvsAMY) #homoscedasticity
vif(PRECvsAMY)   #no multicollinearity
summary(PRECvsAMY)
AIC(PRECvsAMY)
tab_model(PRECvsAMY, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY.doc")
#FIGURE
plot_PRECvsAMY <- ggplot(base, aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox), method = lm) +
  xlab ("Precuneus BOLD-Signal at baseline") + ylab ("Whole Brain Amyloid SURV") +
  PREC_theme()
plot_PRECvsAMY

#BL prec on tau:
PRECvsTAU <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years +Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU)
shapiro.test(rstandard(PRECvsTAU)) #not normally distributed
bptest(PRECvsTAU) #homoscedasticity
vif(PRECvsTAU)   #no multicollinearity
summary(PRECvsTAU)
AIC(PRECvsTAU)
tab_model(PRECvsTAU, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsTAU.doc")

describeBy(base$Entorhinal_Tau_PET_LR_MEAN, base$Sex)
PRECvsTAU_f <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base_f) 
anova(PRECvsTAU_f)


### effects of slope of activation over time?
#model for slope extraction
PREC <- lmer(PREC_Bilat ~ Time + (1 + Time|Subject), data=prec, REML = F)
anova(PREC, type = 3)
# Get the FC slope differences
beta <- coef(PREC) #fixed + random effects ( mean slope + individual deviations)
coef(summary(PREC))[ , "Estimate"]
slope_values <- ranef(PREC) #random effect: individual deviation from mean slope per subject
slope_values <- as.data.frame(slope_values)
slope_prec <- slope_values$condval


#prec slope on amyloid:
PRECvsAMY_SLOPE <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY_SLOPE) #  effect slope
shapiro.test(rstandard(PRECvsAMY_SLOPE)) #not normally distributed
bptest(PRECvsAMY_SLOPE) #homoscedasticity
vif(PRECvsAMY_SLOPE)   #no multicollinearity
summary(PRECvsAMY_SLOPE)
AIC(PRECvsAMY_SLOPE)
tab_model(PRECvsAMY_SLOPE, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_SLOPE.doc")
#FIGURE
plot_PRECvsAMY_SLOPE <- ggplot(base, aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox), method = lm) +
  xlab ("Precuneus BOLD-Signal over time (slope)") + ylab ("Whole Brain Amyloid SURV") +
  PREC_theme()
plot_PRECvsAMY_SLOPE

#prec slope on tau:
PRECvsTAU_SLOPE <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years+Sex+ BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_SLOPE) 
shapiro.test(rstandard(PRECvsTAU_SLOPE)) # almost normally distributed
bptest(PRECvsTAU_SLOPE) #homoscedasticity
vif(PRECvsTAU_SLOPE)   #no multicollinearity
summary(PRECvsTAU_SLOPE)
AIC(PRECvsTAU_SLOPE)
tab_model(PRECvsTAU_SLOPE, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsTAU_SLOPE.doc")


#is baseline activation and slope correlated?
nonan <- subset(base, Subject != "sub-MTL0274" & Subject != "sub-MTL0387" & Subject != "sub-MTL0460" & Subject != "sub-MTL0461" & Subject != "sub-MTL0489" & Subject != "sub-MTL0491" & Subject != "sub-MTL0494" & Subject != "sub-MTL0531" & Subject != "sub-MTL0537" & Subject != "sub-MTL0543" & Subject != "sub-MTL0550" & Subject != "sub-MTL0595" & Subject != "sub-MTL0599" & Subject != "sub-MTL0608")
cor.test(nonan$fMRI_Precuneus_Baseline, nonan$fMRI_Precuneus_Slope)



#########################################################################################################################
### effects of APOE genotype?  (Research question iv))
### effects of baseline activation?
#BL prec on amyloid:
PRECvsAMY_GEN <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Baseline*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY_GEN) 
shapiro.test(rstandard(PRECvsAMY_GEN)) #not normally distributed
bptest(PRECvsAMY_GEN) #homoscedasticity
vif(PRECvsAMY_GEN)   #no multicollinearity
summary(PRECvsAMY_GEN)
AIC(PRECvsAMY_GEN)
tab_model(PRECvsAMY_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_GEN.doc")
#FIGURE
plot_PRECvsAMY_GEN <- ggplot(base, aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox, colour = APOE4_Group)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox), method = lm) +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = expression(italic("APOE4") * " group")) +
  xlab ("Precuneus BOLD-Signal at baseline") + ylab ("Whole Brain Amyloid SURV") +
  PREC_theme()
plot_PRECvsAMY_GEN

#carrier
PRECvsAMY_GEN_carrier <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base_carrier)
anova(PRECvsAMY_GEN_carrier) 
coefficients_summary2 <- summary(PRECvsAMY_GEN_carrier)$coefficients
coefficients_summary2
AIC(PRECvsAMY_GEN_carrier)
tab_model(PRECvsAMY_GEN_carrier, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_GEN_carrier.doc")


#BL prec on tau:
PRECvsTAU_GEN <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Baseline*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years +Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_GEN) 



### effects of slope of activation over time?
#prec slope on amyloid:
PRECvsAMY_SLOPE_GEN <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Slope*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY_SLOPE_GEN) 
shapiro.test(rstandard(PRECvsAMY_SLOPE_GEN)) #not normally distributed
bptest(PRECvsAMY_SLOPE_GEN) #homoscedasticity
vif(PRECvsAMY_SLOPE_GEN)   #no multicollinearity
summary(PRECvsAMY_SLOPE_GEN)
AIC(PRECvsAMY_SLOPE_GEN)
tab_model(PRECvsAMY_SLOPE_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_SLOPE_GEN.doc")
#FIGURE
plot_PRECvsAMY_SLOPE_GEN <- ggplot(base, aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox, colour = APOE4_Group)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox), method = lm) +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = expression(italic("APOE4") * " group")) +
  xlab ("Precuneus BOLD-Signal over time (slope)") + ylab ("Whole Brain Amyloid SURV") +
  PREC_theme()
plot_PRECvsAMY_SLOPE_GEN

#carrier
PRECvsAMY_SLOPE_GEN_carrier <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base_carrier)
anova(PRECvsAMY_SLOPE_GEN_carrier) 
summary(PRECvsAMY_SLOPE_GEN_carrier)
AIC(PRECvsAMY_SLOPE_GEN_carrier)
tab_model(PRECvsAMY_SLOPE_GEN_carrier, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_SLOPE_GEN_carrier.doc")
          
#prec slope on tau:
PRECvsTAU_SLOPE_GEN <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Slope*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years+Sex+ BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_SLOPE_GEN) 



#Now we test for genotype group differences in whole brain amyloid and entorhinal tau:
#Group differences in PET Amyloid?
amy_GEN <- lm(Amyloid_WB_BoxCox ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex, data=base)
summary(amy_GEN)
AIC(amy_GEN)
tab_model(amy_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "amy_GEN.doc")
#FIGURE
amy_GEN_boxplot <- ggplot(base, aes(x = APOE4_Group, y = Amyloid_WB_BoxCox, fill = APOE4_Group)) + 
  geom_boxplot(color = c("#1A85FF", "#D41159"), fill = NA) + 
  geom_jitter(width = 0.25, aes(color = APOE4_Group)) + 
  scale_fill_manual(values = c("Non-carrier" = "#1A85FF", "Carrier" = "#D41159"), labels = c("A", "B")) +
  scale_color_manual(values = c("Non-carrier" = "#1A85FF", "Carrier" = "#D41159")) +
  labs(
    x = expression(italic("APOE4") * " Group"),
    y = "Whole Brain Amyloid SUVR"
  ) +
  scale_x_discrete(labels = c("A" = "Non-carrier", "B"= "Carrier")) +
  PREC_theme() 
amy_GEN_boxplot

#Group differences in PET Entorhinal Tau?
tau_GEN <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex, data=base)
summary(tau_GEN)
AIC(tau_GEN)
tab_model(tau_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "tau_GEN.doc")
#FIGURE
tau_GEN_boxplot <- ggplot(base, aes(x = APOE4_Group, y = Entorhinal_Tau_PET_LR_MEAN_BoxCox, fill = APOE4_Group)) + 
  geom_boxplot(color = c("#1A85FF", "#D41159"), fill = NA) +  
  geom_jitter(width = 0.25, aes(color = APOE4_Group)) +  
  scale_fill_manual(values = c("Non-carrier" = "#1A85FF", "Carrier" = "#D41159"), labels = c("A", "B")) +
  scale_color_manual(values = c("Non-carrier" = "#1A85FF", "Carrier" = "#D41159")) +
  labs(
    x = expression(italic("APOE4") * " Group"),
    y = "Entorhinal Tau SUVR"
  ) +
  scale_x_discrete(labels = c("A" = "Non-carrier", "B"= "Carrier")) +
  PREC_theme() 
tau_GEN_boxplot

#########################################################################################################################
# Research question iii) & iv)
#########################################################################################################################

### COGNITION
################################################################################
### Baseline memory performance
#Note --> NaN have been subset out for subjects with missing data points
base$Sex <- as.numeric(base$Sex)
#BL prec activation - BL corrected hit rate fmri task 
prec_v_CorrHR <- base[,c("CorrHR_Baseline", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_CorrHR, method = "spearman") 

#BL prec activation - BL RBANS
nonan <- subset(base, Subject != "sub-MTL0384" & Subject != "sub-MTL0390")
nonan$Sex <- as.numeric(nonan$Sex)
prec_v_RBANS <- nonan[,c("RBANS_delayed_memory_index_score_Baseline", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_RBANS, method = "spearman") 

#BL prec activation - BL Moca
nonan <- subset(base, Subject != "sub-MTL0550")
nonan$Sex <- as.numeric(nonan$Sex)
prec_v_moca <- nonan[,c("MOCA_Total_Score", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_moca, method = "spearman") 



### longitudinal memory performance
#corrected hit rate over time

CorrHR <- lmer(HR_corr ~ Time + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
anova(CorrHR, type = 3)
summary(CorrHR)
confint(CorrHR)
AIC(CorrHR)
tab_model(CorrHR, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "CorrHR.doc")
describeBy(corrhr_f$HR_corr, corrhr_f$Time)

#RBANS delayed memory index over time
RBANS <- lmer(delayed_memory_index_score ~ Time + Age_Baseline_Months + Sex_String + Education_Baseline_Years + (1 |Subject), data=rbans, REML = F)
anova(RBANS, type = 3)
Anova(RBANS)
summary(RBANS)
confint(RBANS)
AIC(RBANS)
tab_model(RBANS, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RBANS.doc")



### COGNITION, PRECUNEUS ACTIVATION, APOE
################################################################################
rbans$Time <- as.factor(rbans$Time)
corrhr$Time <- as.factor(corrhr$Time)

# longitudinal corrHR - BL activation, time, apoe, age, sex, education

PREC_GEN_TIME_CorrHR <- lmer(HR_corr ~ Time*APOE4_Group + Time*PREC_Bilat_Baseline + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
anova(PREC_GEN_TIME_CorrHR) 
summary(PREC_GEN_TIME_CorrHR)
confint(PREC_GEN_TIME_CorrHR)
AIC(PREC_GEN_TIME_CorrHR)
tab_model(PREC_GEN_TIME_CorrHR, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC_GEN_TIME_CorrHR.doc")

PREC_GEN_TIME_CorrHR2 <- lmer(HR_corr ~ Time*APOE4_Group*PREC_Bilat_Baseline + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
anova(PREC_GEN_TIME_CorrHR2) 
summary(PREC_GEN_TIME_CorrHR2)
confint(PREC_GEN_TIME_CorrHR2)
AIC(PREC_GEN_TIME_CorrHR2)
tab_model(PREC_GEN_TIME_CorrHR2, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC_GEN_TIME_CorrHR.doc")

#post hoc corrHR
emm_act_corrhr <- emmeans(PREC_GEN_TIME_CorrHR,  ~ Time * APOE4_Group)
contrast_results_corrhr <- contrast(emm_act_corrhr, method = "pairwise", by = "Time")
post_hoc_corrhr <- summary(contrast_results_corrhr, infer = c(TRUE, TRUE))
post_hoc_corrhr <- as.data.frame(post_hoc_corrhr) %>% mutate_if(is.numeric, ~ round(.,3))
#table #report
post_hoc_corrhr_to_save <- post_hoc_corrhr %>% regulartable() %>% autofit()
word_corrhr <- read_docx() %>%
  body_add_flextable(post_hoc_corrhr_to_save) %>%
  print(target = "post_hoc_corrhr.docx")
#Fig
emmip(emm_act_corrhr, APOE4_Group ~ Time, CIs = TRUE) +
  labs(y = "Corrected Hit Rate",x = "Session") +
  scale_color_manual(values = c("A" = "#1A85FF", "B" = "#D41159"),labels = c("A" = "Non-carrier", "B" = "Carrier"))+ labs(color = expression(italic("APOE4") * " group")) +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme() 


# longitudinal RBANS - BL activation, time, apoe, age, sex, education

PREC_GEN_TIME_RBANS3 <- lmer(delayed_memory_index_score ~ Time*APOE4_Group + Time*fMRI_Precuneus_Baseline + Age_Baseline_Months + Sex + Education_Baseline_Years + (1 |Subject), data=rbans, REML = F)
anova(PREC_GEN_TIME_RBANS3)
Anova(PREC_GEN_TIME_RBANS3)
summary(PREC_GEN_TIME_RBANS3)
tab_model(PREC_GEN_TIME_RBANS3, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC_GEN_TIME_RBANS3.doc")

PREC_GEN_TIME_RBANS4 <- lmer(delayed_memory_index_score ~ Time*APOE4_Group*fMRI_Precuneus_Baseline + Age_Baseline_Months + Sex + Education_Baseline_Years + (1 |Subject), data=rbans, REML = F)
anova(PREC_GEN_TIME_RBANS4) 
Anova(PREC_GEN_TIME_RBANS4)
summary(PREC_GEN_TIME_RBANS4)
tab_model(PREC_GEN_TIME_RBANS4, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC_GEN_TIME_RBANS4.doc")

# post hoc RBANS
#emm <- emmeans(PREC_GEN_TIME_RBANS4, pairwise ~ Time * APOE4_Group * fMRI_Precuneus_Baseline)
#emmip(emm, APOE4_Group ~ Time | fMRI_Precuneus_Baseline, CIs = TRUE)
lower_quartile <- quantile(rbans$fMRI_Precuneus_Baseline, 0.25)
upper_quartile <- quantile(rbans$fMRI_Precuneus_Baseline, 0.75)
emm_act <- emmeans(PREC_GEN_TIME_RBANS4, ~ Time *APOE4_Group | fMRI_Precuneus_Baseline, at = list(fMRI_Precuneus_Baseline = c(lower_quartile, upper_quartile)))
summary(emm_act)
#contrast_results <- contrast(emm_act, method = "pairwise")
contrast_results <- contrast(emm_act, method = "pairwise", by = c("Time", "fMRI_Precuneus_Baseline"))
post_hoc_rbans <- summary(contrast_results, infer = c(TRUE, TRUE))
post_hoc_rbans <- as.data.frame(post_hoc_rbans) %>% mutate_if(is.numeric, ~ round(.,3))
#table #report
post_hoc_rbans_to_save <- post_hoc_rbans %>% regulartable() %>% autofit()
word <- read_docx() %>%
  body_add_flextable(post_hoc_rbans_to_save) %>%
  print(target = "post_hoc_rbans.docx")
#SUPP FIG
emmip(emm_act, APOE4_Group ~ Time | fMRI_Precuneus_Baseline, CIs = TRUE, at = list(fMRI_Precuneus_Baseline = c(lower_quartile, upper_quartile))) +
  labs(y = "RBANS delayed Memory Index Score",x = "Session",title = "Activity at lower and upper quartile") +
  scale_color_manual(values = c("A" = "#1A85FF", "B" = "#D41159"),labels = c("A" = "Non-carrier", "B" = "Carrier"))+ labs(color = expression(italic("APOE4") * " group")) +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme() 

#FIGURE
mean_prec_bilat <- mean(rbans$fMRI_Precuneus_Baseline, na.rm = TRUE)
rbans <- rbans %>%
  mutate(APOE_ACT_group = case_when(
    APOE4_Group == "A" & fMRI_Precuneus_Baseline > mean_prec_bilat ~ "A_high",
    APOE4_Group == "A" & fMRI_Precuneus_Baseline <= mean_prec_bilat ~ "A_low",
    APOE4_Group == "B" & fMRI_Precuneus_Baseline > mean_prec_bilat ~ "B_high",
    APOE4_Group == "B" & fMRI_Precuneus_Baseline <= mean_prec_bilat ~ "B_low",
    TRUE ~ NA_character_
  ))

plot_COG_APOE_PREC_TIME <- ggplot(rbans, aes(x = as.factor(Time), y = delayed_memory_index_score, colour = APOE_ACT_group)) +
  geom_point(shape = 1, color = "white") +
  geom_smooth(aes(x = as.numeric(Time), y = delayed_memory_index_score, fill = APOE_ACT_group), method = lm, alpha = 0.1, size = 1.5) + 
  xlab("Session") + 
  ylab("RBANS delayed Memory Index Score") +
  scale_color_manual(labels = c("Non-carrier & high activation", "Non-carrier & low activation", "Carrier & high activation", "Carrier & low activation"), values = c("#1A85FF", "#80C1FF", "#D41159", "#FF6F61")) + 
  scale_fill_manual(values = c("#1A85FF", "#80C1FF", "#D41159", "#FF6F61"), guide = "none") + 
  labs(color = "Group") + 
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  coord_cartesian(ylim = c(100, 110)) +
  PREC_theme()
plot_COG_APOE_PREC_TIME



##########################################################################################################################

#Analysis of PREVENT-AD data - baseline fMRI and follow-up PET with APOE genotype & cognition
#Larissa Fischer & Eóin N. Molloy - Multimodal Neuroimaging Lab, DZNE Magdeburg

#update after peer-review process

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
#overview cohort
describe(base$Age_Baseline_Years) 
describe(base$Education_Baseline_Years)
summary(base$Sex)
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

#overview amyloid positive
base$Amyloid_Group <- factor(ifelse(base$Amyloid_WB > 1.39, 1, 0))
summary(base$Amyloid_Group)
table(base$Amyloid_Group, base$APOE4_Group)
chi_test_Abeta_APOE <- chisq.test(base$Amyloid_Group, base$APOE4_Group)
print(chi_test_Abeta_APOE)
t.test(base$BaseMRI_to_Amy_PET_Months ~ base$APOE4_Group)

#RBANS memory and fMRI
#BL
cor.test(base$RBANS_delayed_memory_index_score_Baseline, base$CorrHR_Baseline, method = "pearson")
ggplot(base, aes(x = CorrHR_Baseline, y = RBANS_delayed_memory_index_score_Baseline)) +
  geom_point() +  
  geom_smooth(method = "lm", col = "blue") +
  labs(x = "fMRI Task Baseline", y = "RBANS Delayed Memory Index Score Baseline") +
  theme_minimal()
#slope
rbans_1 <- rbans_1 %>%
  rename(Slope_rbans = Slope)
merged_cog_slopes <- corrhr_1 %>%
  left_join(rbans_1, by = "Subject")
cor.test(merged_cog_slopes$Slope, merged_cog_slopes$Slope_rbans, method = "pearson")
ggplot(merged_cog_slopes, aes(x = Slope, y = Slope_rbans)) +
  geom_point() +  
  geom_smooth(method = "lm", col = "blue") + 
  labs(x = "fMRI Task Slope", y = "RBANS Delayed Memory Index Score Slope") +
  theme_minimal()

#RBANS and fNRI task performance per timepoint
describeBy(rbans$delayed_memory_index_score, rbans$Session)
fmri_performance_sample <- subset(fmri_performance, Subject %in% base$Subject)
fmri_performance_sample$Subject<- droplevels(fmri_performance_sample$Subject)
describeBy(fmri_performance_sample$Retr_corrHR, fmri_performance_sample$Session)
describeBy(fmri_performance_sample$Retr_HR, fmri_performance_sample$Session)
fmri_performance_sample$Old_hit <- fmri_performance_sample$ObjRec_Fam + fmri_performance_sample$SourceRecoll + fmri_performance_sample$SourceMisattr
describeBy(fmri_performance_sample$Old_hit, fmri_performance_sample$Session)
fmri_performance_sample$Retr_MR <- fmri_performance_sample$Old_miss/48 #miss rate
describeBy(fmri_performance_sample$Old_miss, fmri_performance_sample$Session)
fmri_performance_sample$Retr_CRR <- fmri_performance_sample$New_CR/48 #correct rejection rate
describeBy(fmri_performance_sample$New_CR, fmri_performance_sample$Session)
describe(fmri_performance_sample$Old_miss)

#overview sample size development with criteria
#start: 443 subjects
fmri_performance <- subset(fmri_performance, Exclude_anat == 0)
fmri_performance$Subject<- droplevels(fmri_performance$Subject) #442 Subjects
fmri_performance <- subset(fmri_performance, Exclude_func == 0)
fmri_performance$Subject<- droplevels(fmri_performance$Subject) #383 Subjects
fmri_performance <- subset(fmri_performance, Exclude_nolog == 0)
fmri_performance$Subject<- droplevels(fmri_performance$Subject) #374 Subjects, 1104 observations
fmri_performance <- subset(fmri_performance, Retr_PercMissing != 100)
fmri_performance$Subject<- droplevels(fmri_performance$Subject) #374 Subjects, 1046 observations
fmri_performance <- subset(fmri_performance, Retr_corrHR > 0.2)
fmri_performance$Subject<- droplevels(fmri_performance$Subject) #358 Subjects, 991 observations
fmri_performance <- subset(fmri_performance, (ObjRec_Fam + SourceRecoll + SourceMisattr) > 10 & New_CR > 10)
fmri_performance$Subject<- droplevels(fmri_performance$Subject) #still 358 Subjects (but 2 sessions less)

subjects_in_FMRI_not_in_base <- setdiff(fmri_performance$Subject, base$Subject)
subjects_in_base_not_in_FMRI <- setdiff(base$Subject, fmri_performance$Subject)
length(subjects_in_FMRI_not_in_base) #193
length(subjects_in_base_not_in_FMRI) #0 (sanity check)



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
describeBy(base$fMRI_Precuneus_Baseline, base$Sex)

#effects of change in precuneus activation over time?
#first only Time, then add Time*APOE4_Group and Time*Sex to model
PREC <- lmer(PREC_Bilat ~ Time*Sex + Time*APOE4_Group +Age + Education + PREC_GMV + (1 |Subject), data=prec, REML = F)
summary(PREC)
confint(PREC)
tab_model(PREC, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECwobistdu.doc")
#male participants
PREC_m <- lmer(PREC_Bilat ~ Time +Age + Education + PREC_GMV + (1 |Subject), data=prec_m, REML = F)
summary(PREC_m)
confint(PREC_m)
tab_model(PREC_m, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC.doc")
#female participants
PREC_f <- lmer(PREC_Bilat ~ Time +Age + Education + PREC_GMV + (1 |Subject), data=prec_f, REML = F)
summary(PREC_f)
confint(PREC_f)
tab_model(PREC_f, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PREC.doc")

#FIGURE 2A precuneus activity & APOE
plot_PREC_APOE <- ggplot(prec, aes(as.factor(Time), PREC_Bilat, colour = APOE4_Group)) + 
  geom_jitter(width = 0.25) +
  geom_smooth(aes(as.numeric(Time), PREC_Bilat), method = lm) +
  xlab ("Session") + ylab ("Precuneus BOLD-Signal") +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = expression(italic("APOE4") * " group")) +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme()
plot_PREC_APOE

#extended data table 3-1  precuneus activity & sex
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

#First we test for genotype group differences in whole brain amyloid and entorhinal tau:
#Group differences in PET Amyloid?
amy_GEN <- lm(Amyloid_WB_BoxCox ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex, data=base)
summary(amy_GEN)
tab_model(amy_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "amy_GEN.doc")
#FIGURE 2B
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
tab_model(tau_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "tau_GEN.doc")
#FIGURE 2C
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


#PRECUNEUS ACTIVITY & PATHOLOGY
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
tab_model(PRECvsAMY, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY.doc")
#FIGURE 3A
plot_PRECvsAMY <- ggplot(base, aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox), method = lm) +
  xlab ("Precuneus BOLD-Signal at baseline") + ylab ("Whole Brain Amyloid SUVR") +
  PREC_theme()
plot_PRECvsAMY

#BL prec on tau:
PRECvsTAU <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years +Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU)
shapiro.test(rstandard(PRECvsTAU)) #not normally distributed
bptest(PRECvsTAU) #homoscedasticity
vif(PRECvsTAU)   #no multicollinearity
summary(PRECvsTAU)
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
tab_model(PRECvsAMY_SLOPE, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_SLOPE.doc")
#FIGURE 3B
plot_PRECvsAMY_SLOPE <- ggplot(base, aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox), method = lm) +
  xlab ("Precuneus BOLD-Signal over time (slope)") + ylab ("Whole Brain Amyloid SUVR") +
  PREC_theme()
plot_PRECvsAMY_SLOPE

#prec slope on tau:
PRECvsTAU_SLOPE <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years+Sex+ BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_SLOPE) 
shapiro.test(rstandard(PRECvsTAU_SLOPE)) # almost normally distributed
bptest(PRECvsTAU_SLOPE) #homoscedasticity
vif(PRECvsTAU_SLOPE)   #no multicollinearity
summary(PRECvsTAU_SLOPE)
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
tab_model(PRECvsAMY_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_GEN.doc")
#FIGURE 3C
plot_PRECvsAMY_GEN <- ggplot(base, aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox, colour = APOE4_Group)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Baseline, Amyloid_WB_BoxCox), method = lm) +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = expression(italic("APOE4") * " group")) +
  xlab ("Precuneus BOLD-Signal at baseline") + ylab ("Whole Brain Amyloid SUVR") +
  PREC_theme()
plot_PRECvsAMY_GEN

#carrier
PRECvsAMY_GEN_carrier <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base_carrier)
anova(PRECvsAMY_GEN_carrier) 
coefficients_summary2 <- summary(PRECvsAMY_GEN_carrier)$coefficients
coefficients_summary2
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
tab_model(PRECvsAMY_SLOPE_GEN, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_SLOPE_GEN.doc")
#FIGURE 3D
plot_PRECvsAMY_SLOPE_GEN <- ggplot(base, aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox, colour = APOE4_Group)) + 
  geom_point() +
  geom_smooth(aes(fMRI_Precuneus_Slope, Amyloid_WB_BoxCox), method = lm) +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + labs(color = expression(italic("APOE4") * " group")) +
  xlab ("Precuneus BOLD-Signal over time (slope)") + ylab ("Whole Brain Amyloid SUVR") +
  PREC_theme()
plot_PRECvsAMY_SLOPE_GEN

#carrier
PRECvsAMY_SLOPE_GEN_carrier <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base_carrier)
anova(PRECvsAMY_SLOPE_GEN_carrier) 
summary(PRECvsAMY_SLOPE_GEN_carrier)
tab_model(PRECvsAMY_SLOPE_GEN_carrier, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "PRECvsAMY_SLOPE_GEN_carrier.doc")
          
#prec slope on tau:
PRECvsTAU_SLOPE_GEN <- lm(Entorhinal_Tau_PET_LR_MEAN_BoxCox ~ fMRI_Precuneus_Slope*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years+Sex+ BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_SLOPE_GEN) 



#########################################################################################################################
# Research question iii) & iv)
#########################################################################################################################

### COGNITION
################################################################################
### Baseline memory performance
#Note --> NaN have been subset out for subjects with missing data points
#BL prec activation - BL corrected hit rate fmri task 
prec_v_CorrHR <- base[,c("CorrHR_Baseline", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_CorrHR, method = "spearman") 

#BL prec activation - BL RBANS
nonan <- subset(base, Subject != "sub-MTL0384" & Subject != "sub-MTL0390")
prec_v_RBANS <- nonan[,c("RBANS_delayed_memory_index_score_Baseline", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_RBANS, method = "spearman") 


### longitudinal memory performance
#corrected hit rate over time

CorrHR <- lmer(HR_corr ~ Time + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
summary(CorrHR)
confint(CorrHR)
tab_model(CorrHR, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "CorrHR.doc")

#RBANS delayed memory index over time
RBANS <- lmer(delayed_memory_index_score ~ Time + Age_Baseline_Months + Sex_String + Education_Baseline_Years + (1 |Subject), data=rbans, REML = F)
summary(RBANS)
confint(RBANS)
tab_model(RBANS, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "RBANS.doc")



### COGNITION, PRECUNEUS ACTIVATION, APOE
################################################################################
rbans$Session <- as.factor(rbans$Session)
corrhr$Session <- as.factor(corrhr$Session)


# corrected hit rate
################################################################################
# longitudinal corrHR - BL activation, time (session), apoe, age, sex, education

PREC_GEN_TIME_CorrHR <- lmer(HR_corr ~ Session*APOE4_Group + Session*PREC_Bilat_Baseline + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
anova(PREC_GEN_TIME_CorrHR) 
summary(PREC_GEN_TIME_CorrHR)
confint(PREC_GEN_TIME_CorrHR)

PREC_GEN_TIME_CorrHR2 <- lmer(HR_corr ~ Session*APOE4_Group*PREC_Bilat_Baseline + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
anova(PREC_GEN_TIME_CorrHR2) 
summary(PREC_GEN_TIME_CorrHR2)
confint(PREC_GEN_TIME_CorrHR2)

#post hoc corrHR
emm_act_corrhr <- emmeans(PREC_GEN_TIME_CorrHR,  ~ Session * APOE4_Group)
contrast_results_corrhr <- contrast(emm_act_corrhr, method = "pairwise", by = "Session")
post_hoc_corrhr <- summary(contrast_results_corrhr, infer = c(TRUE, TRUE))
post_hoc_corrhr <- as.data.frame(post_hoc_corrhr) %>% mutate_if(is.numeric, ~ round(.,3))
#table #report
post_hoc_corrhr_to_save <- post_hoc_corrhr %>% regulartable() %>% autofit()
word_corrhr <- read_docx() %>%
  body_add_flextable(post_hoc_corrhr_to_save) %>%
  print(target = "post_hoc_corrhr.docx")

#FIGURE 4 A 
plot_COG_APOE_TIME <- ggplot(corrhr, aes(x = Time, y = HR_corr, colour = APO4_Group)) +
  geom_point(shape = 1, color = "white") +
  geom_smooth(aes(x = Time, y = HR_corr, fill = APO4_Group), method = lm, alpha = 0.1, size = 1.5) + 
  xlab("Time") + 
  ylab("Corrected hit rate") +
  scale_color_manual(labels = c("Non-carrier", "Carrier"), values = c("#1A85FF", "#D41159")) + 
  scale_fill_manual(values = c("#1A85FF", "#D41159"), guide = "none") + 
  labs(color = "Group") +
  scale_x_continuous(breaks = c(-10, 120, 365, 750, 1490), 
                     labels = c("BL", "3 Months", "12 Months", "24 Months", "48 Months")) +
  coord_cartesian(ylim = c(0.5, 0.8)) +
  PREC_theme()
plot_COG_APOE_TIME

#extended data Fig.4-1 A 
emmip(emm_act_corrhr, APOE4_Group ~ Time, CIs = TRUE) +
  labs(y = "Corrected Hit Rate",x = "Session") +
  scale_color_manual(values = c("A" = "#1A85FF", "B" = "#D41159"),labels = c("A" = "Non-carrier", "B" = "Carrier"))+ labs(color = expression(italic("APOE4") * " group")) +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme() 



# RBANS
################################################################################
# longitudinal RBANS - BL activation, time, apoe, age, sex, education

PREC_GEN_TIME_RBANS3 <- lmer(delayed_memory_index_score ~ Session*APOE4_Group + Session*fMRI_Precuneus_Baseline + Age_Baseline_Months + Sex + Education_Baseline_Years + (1 |Subject), data=rbans, REML = F)
anova(PREC_GEN_TIME_RBANS3)
summary(PREC_GEN_TIME_RBANS3)

PREC_GEN_TIME_RBANS4 <- lmer(delayed_memory_index_score ~ Session*APOE4_Group*fMRI_Precuneus_Baseline + Age_Baseline_Months + Sex + Education_Baseline_Years + (1 |Subject), data=rbans, REML = F)
anova(PREC_GEN_TIME_RBANS4) 
summary(PREC_GEN_TIME_RBANS4)

# post hoc RBANS
#emm <- emmeans(PREC_GEN_TIME_RBANS4, pairwise ~ Session * APOE4_Group * fMRI_Precuneus_Baseline)
#emmip(emm, APOE4_Group ~ Session | fMRI_Precuneus_Baseline, CIs = TRUE)
lower_quartile <- quantile(rbans$fMRI_Precuneus_Baseline, 0.25)
upper_quartile <- quantile(rbans$fMRI_Precuneus_Baseline, 0.75)
emm_act <- emmeans(PREC_GEN_TIME_RBANS4, ~ Session *APOE4_Group | fMRI_Precuneus_Baseline, at = list(fMRI_Precuneus_Baseline = c(lower_quartile, upper_quartile)))
summary(emm_act)
#contrast_results <- contrast(emm_act, method = "pairwise")
contrast_results <- contrast(emm_act, method = "pairwise", by = c("Session", "fMRI_Precuneus_Baseline"))
post_hoc_rbans <- summary(contrast_results, infer = c(TRUE, TRUE))
post_hoc_rbans <- as.data.frame(post_hoc_rbans) %>% mutate_if(is.numeric, ~ round(.,3))
#table #report
post_hoc_rbans_to_save <- post_hoc_rbans %>% regulartable() %>% autofit()
word <- read_docx() %>%
  body_add_flextable(post_hoc_rbans_to_save) %>%
  print(target = "post_hoc_rbans.docx")

#FIGURE 4 B
mean_prec_bilat <- mean(rbans$fMRI_Precuneus_Baseline, na.rm = TRUE)
rbans <- rbans %>%
  mutate(APOE_ACT_group = case_when(
    APOE4_Group == "A" & fMRI_Precuneus_Baseline > mean_prec_bilat ~ "A_high",
    APOE4_Group == "A" & fMRI_Precuneus_Baseline <= mean_prec_bilat ~ "A_low",
    APOE4_Group == "B" & fMRI_Precuneus_Baseline > mean_prec_bilat ~ "B_high",
    APOE4_Group == "B" & fMRI_Precuneus_Baseline <= mean_prec_bilat ~ "B_low",
    TRUE ~ NA_character_
  ))

plot_COG_APOE_PREC_TIME <- ggplot(rbans, aes(x = Time, y = delayed_memory_index_score, colour = APOE_ACT_group)) +
  geom_point(shape = 1, color = "white") +
  geom_smooth(aes(x = Time, y = delayed_memory_index_score, fill = APOE_ACT_group), method = lm, alpha = 0.1, size = 1.5) + 
  xlab("Time") + 
  ylab("RBANS delayed memory index score") +
  scale_color_manual(labels = c("Non-carrier & high activation", "Non-carrier & low activation", "Carrier & high activation", "Carrier & low activation"), values = c("#1A85FF", "#80C1FF", "#D41159", "#FF6F61")) + 
  scale_fill_manual(values = c("#1A85FF", "#80C1FF", "#D41159", "#FF6F61"), guide = "none") + 
  labs(color = "Group")+
  scale_x_continuous(breaks = c(-10, 120, 365, 750, 1490), 
                     labels = c("BL", "3 Months", "12 Months", "24 Months", "48 Months")) +
  coord_cartesian(ylim = c(100, 110)) +
  PREC_theme()
plot_COG_APOE_PREC_TIME

#extended data Fig.4-1 B 
emmip(emm_act, APOE4_Group ~ Time | fMRI_Precuneus_Baseline, CIs = TRUE, at = list(fMRI_Precuneus_Baseline = c(lower_quartile, upper_quartile))) +
  labs(y = "RBANS delayed Memory Index Score",x = "Session",title = "Activity at lower and upper quartile") +
  scale_color_manual(values = c("A" = "#1A85FF", "B" = "#D41159"),labels = c("A" = "Non-carrier", "B" = "Carrier"))+ labs(color = expression(italic("APOE4") * " group")) +
  scale_x_discrete(labels = c("1" = "Baseline", "2" = "3 Months", "3" = "12 Months", "4" = "24 Months", "5" = "48 Months")) +
  PREC_theme() 


#comparison between slopes
################################################################################
rbans$Time_scaled <- scale(rbans$Time)
corrhr$Time_scaled <- scale(corrhr$Time)


####corrhr 
#interaction --> 2 slopes different?
groups <- c("A", "B")
slopes_df <- data.frame()
# Loop through each group
for (group in groups) {
  corrhr_group <- subset(corrhr, APOE4_Group == group)
  COG_APOE_TIME <- lmer(HR_corr ~ Time_scaled + (1 + Time_scaled| Subject), data = corrhr_group, REML = FALSE) 
  # Extract slopes
  slopes <- coef(COG_APOE_TIME)$Subject
  slopes <- data.frame(slopes)
  slopes$Subject <- rownames(slopes)
  slopes$Group <- group
  slopes_df <- rbind(slopes_df, slopes)
}
colnames(slopes_df)[which(colnames(slopes_df) == "Time_scaled")] <- "Slope"
corrhr <- merge(corrhr, slopes_df[, c("Subject", "Slope")], by = "Subject", all.x = TRUE)

corrhr_1 <- corrhr  %>% filter(Session == "1")
describeBy(corrhr_1$Slope, corrhr_1$APOE4_Group)
anova_results <- aov(Slope ~ APOE4_Group + Sex + Education + Age, data = corrhr_1)
anova_results

SLOPE_CORRHR_APOE <- lm(Slope ~ APOE4_Group + Age + Sex + Education, data = corrhr_1)
anova(SLOPE_CORRHR_APOE)
summary(SLOPE_CORRHR_APOE)
tab_model(SLOPE_CORRHR_APOE, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "SLOPE_CORRHR_APOE.doc")


####RBANS
#interaction, non-carriers low activation slope different from other 3 group slopes?
#get slopes (carrier - high act, ...) #this time mean slope + individual difference
rbans$APOE_ACT_group <- as.factor(rbans$APOE_ACT_group)
rbans_A_low <- subset(rbans, APOE_ACT_group == "A_low")
rbans_A_high <- subset(rbans, APOE_ACT_group == "A_high")
rbans_B_low <- subset(rbans, APOE_ACT_group == "B_low") 
rbans_B_high <-subset(rbans, APOE_ACT_group == "B_high")

groups <- c("A_low", "A_high", "B_low", "B_high")
slopes_df <- data.frame()
# Loop through each group
for (group in groups) {
  rbans_group <- subset(rbans, APOE_ACT_group == group)
  COG_APOE_PREC_TIME <- lmer(delayed_memory_index_score ~ Time_scaled + (1 + Time_scaled | Subject), data = rbans_group, REML = FALSE) 
  # Extract slopes
  slopes <- coef(COG_APOE_PREC_TIME)$Subject
  slopes <- data.frame(slopes)
  slopes$Subject <- rownames(slopes)
  slopes$Group <- group
  slopes_df <- rbind(slopes_df, slopes)
}
colnames(slopes_df)[which(colnames(slopes_df) == "Time_scaled")] <- "Slope"
rbans <- merge(rbans, slopes_df[, c("Subject", "Slope")], by = "Subject", all.x = TRUE)

rbans_1 <- rbans  %>% filter(Visit_Label_RBANS == "RBANS_BL00")
describeBy(rbans_1$Slope, rbans_1$APOE_ACT_group)
SLOPE_RBANS_APOE_ACT <- lm(Slope ~ as.numeric(APOE_ACT_group) + Age_Baseline_Years + Sex + Education_Baseline_Years, data = rbans_1)
summary(SLOPE_RBANS_APOE_ACT)
tab_model(SLOPE_RBANS_APOE_ACT, df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "SLOPE_RBANS_APOE_ACT.doc")
anova_results <- aov(Slope ~ APOE_ACT_group + Age_Baseline_Years + Sex + Education_Baseline_Years, data = rbans_1)
summary(anova_results)
tukey_results <- TukeyHSD(anova_results, "APOE_ACT_group")
print(tukey_results)



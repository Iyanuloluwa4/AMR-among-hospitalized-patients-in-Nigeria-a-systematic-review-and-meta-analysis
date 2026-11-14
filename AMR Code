#==============================
# INSTALLATION OF PACKAGES
#==============================
install.packages(c("meta", "metafor", "readr", "dplyr", "googlesheets4"))

rm(list = ls())
gc()

library(meta)
library(grid)
library(metafor)
library(readr)
library(dplyr)
library(googlesheets4)

#================================
# Data Handling
#================================

amr_url <-'https://docs.google.com/spreadsheets/d/1vRPDuOw_9ikflj8yvlNHxxr_XCCKxY3kAhboR0m-TJE/export?format=csv'

amr_data <- read_csv(amr_url)

#amr_data <- amr_data %>%
  #filter(Original == 'O')

#View(amr_data)


#==============================
# Recoding and Data Cleansing
#==============================

amr_data <- amr_data %>%
  mutate(
    study = `First author, year of publication`,
    last_research_year = `Last year of research`,
    phase = `NAP phase (1.0 vs 2.0)`,
    design = `Study design`,
    region = `Geopolitical zone`,
    age_group = `Age group(Years)`,
    sample_size = `Sample size`,
    isolate_sample_size = `Sample size for isolate(n)`,
    in_patient_percent = `Percentage of in-patient (%)`,
    average_age = `Average age`,
    ethnicity = Ethnicity,
    facility_type = `Facility type`,
    isolates = `Bacteria isolates`,
    antibiotic = `Antibiotic tested`,
    anti_class = `Antibiotic_class`,
    specimen_type = `Specimen type`,
    ward = Ward,
    sub_specialty = `Sub-Specialty`,
    diagnostic_method = `Diagnostic method`,
    hiv = `HIV (n)`,
    ckd = `CKD at presentation (n)`,
    self_med = `Self-medication practice (%)`,
    prevalence = `Prevalence of AMR (%)`,
    resistance_rate = `Resistance rate (per antibiotic)`,
    susceptibility_rate = `Susceptibility rate`,
    risk_of_bias_score = `Risk of bias score (JBI)`,
    risk_of_bias = `Risk of bias classification`,
    days_in_hospital = `Days spent in hospital (Days)`,
    sepsis = `Sepsis (n)`,
    rob1 = `Appropriate Sample Size`,
    rob2 = `Study Participants Sampling`,
    rob3 = `Adequate Sample Size`,
    rob4 = `Detailed Study Subjects`,
    rob5 = `Sufficient Coverage of Sample`,
    rob6 = `Valid Methods used`,
    rob7 = `Standard Measurement`,
    rob8 = `Appropriate Statistical Analysis`,
    rob9 = `Adequate Response Rate`,
    study_id = study_id,
    effect_id = effect_id,
    ast = `Antimicrobial susceptibility measure`
  )

amr_data <- amr_data %>%
  mutate(
    susceptibility_rate = as.numeric(susceptibility_rate),
    resistance_rate = as.numeric(resistance_rate),
    prevalence = as.numeric(prevalence),
    sample_size = as.numeric(sample_size)
  )

amr_data$no_of_cases<-((amr_data$resistance_rate/100) * (amr_data$sample_size))

#==================================
# Effect size calculation
#==================================

escalc_data <- escalc(
  measure = 'PLO',
  xi = no_of_cases,
  ni = sample_size,
  data = amr_data
)

gc()

#==========================================
# Meta-analysis result [FOR AMR PREVALENCE]
#==========================================

meta_result <- rma(
  yi = yi,
  vi = vi,
  data = escalc_data,
  method = 'REML'
)

result <- predict(meta_result, transf = transf.ilogit)
print(result)

#============================================
# Meta-regression result [FOR AMR PREVALENCE]
#============================================
escalc_data$phase <- factor(escalc_data$phase)

metareg_result <- rma(
  yi = yi,
  vi = vi,
  mods = ~ region,
  data = escalc_data,
  method = 'REML'
)

print(metareg_result)


#===========================================================================
# Meta-analysis result [MULTIVARIATE ANALYSIS FOR ORGANISM-LEVEL RESISTANCE]
#===========================================================================

meta_result <- rma.mv(
  yi = yi,
  V = vi,
  random = ~1 | study_id/effect_id,
  data = escalc_data,
  method = 'REML'
)


#===============================================
# Meta-regression result [MULTIVARIATE ANALYSIS]
#===============================================
escalc_data$phase <- factor(escalc_data$phase)

metareg_result <- rma.mv(
  yi = yi,
  V = vi,
  mods = ~ phase,
  random = ~1 | study_id/effect_id,
  data = escalc_data,
  method = 'REML'
)

result <- summary(metareg_result)
predict(metareg_result, transf = transf.ilogit)
print(result)

#==================================
# Forest Plot FOR AMR PREVALENCE
#==================================
pdf ('amrforestnew.pdf', width = 15, height = 28)

forest(
  meta_result,
  slab = study,
  xlab = "Prevalence of AMR in Hospitalized Patients in Nigeria (%)",
  transf = transf.ilogit,
  refline = NA,
  alim = c(0, 1),
  at = seq(0, 1, 0.1),
  cex = 0.9,
  cex.lab = 1.1,
  cex.axis = 1,
  col = "steelblue",
  pch = 15,
  efac = 1.2,
  addfit = TRUE,
  digits = 2
)

dev.off()

#================================
#Funnel Plot
#===============================
png("funnel_amr.png", width= 6000, height=5000, res= 900)
funnel(meta_result)
dev.off()

#==============================
#Sensitivity Analyses
#=============================

egger_test <- regtest(meta_result, model = 'rma')
print(egger_test)

trimfill(meta_result)

leave1out(meta_result)












#====================================
#SuBgroup Analyses for AMR Prevalence
#====================================

#Region

NC <- escalc_data %>%
  filter(region == 'North-central')
NC_result <- rma(yi = yi, vi = vi, data = NC, method = 'REML')
print(NC_result)
predict(NC_result, transf = transf.ilogit)

SS <- escalc_data %>%
  filter(region == 'South-south')
SS_result <- rma(yi = yi, vi = vi, data = SS, method = 'REML')
print(SS_result)
predict(SS_result, transf = transf.ilogit)

NW <- escalc_data %>%
  filter(region == 'North-west')
NW_result <- rma(yi = yi, vi = vi, data = NW, method = 'REML')
print(NW_result)
predict(NW_result, transf = transf.ilogit)

SE <- escalc_data %>%
  filter(region == 'South-east')
SE_result <- rma(yi = yi, vi = vi, data = SE, method = 'REML')
print(SE_result)
predict(SE_result, transf = transf.ilogit)

NE <- escalc_data %>%
  filter(region == 'North-east')
NE_result <- rma(yi = yi, vi = vi, data = NE, method = 'REML')
print(NE_result)
predict(NE_result, transf = transf.ilogit)


#Facility Type & Ward

Private <- escalc_data%>%
  filter(facility_type == 'Private' | facility_type == 'private')
Private_result <- rma(yi = yi, vi = vi, data = Private, method = 'REML')
print(Private_result)
predict(Private_result, transf = transf.ilogit) 

Obs <- escalc_data%>%
  filter(Ward == 'Obs/gyn')
Obs_result <- rma(yi = yi, vi = vi, data = Obs, method = 'REML')
print(Obs_result)
predict(Obs_result, transf = transf.ilogit)#

ICU <- escalc_data%>%
  filter(Ward == 'ICU')
ICU_result <- rma(yi = yi, vi = vi, data = ICU, method = 'REML')
print(ICU_result)
predict(ICU_result, transf = transf.ilogit)

Pediatric <- escalc_data%>%
  filter(Ward == 'Neonatal')
Pediatric_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Pediatric, method = 'REML')
print(Pediatric_result)
predict(Pediatric_result, transf = transf.ilogit)

Neonatal <- escalc_data%>%
  filter(Ward == 'Neonatal')
Neonatal_result <- rma(yi = yi, vi = vi, data = Neonatal, method = 'REML')
print(Neonatal_result)
predict(Neonatal_result, transf = transf.ilogit)



#===============================================
#SuBgroup Analyses for ISOLATE-LEVEL RESISTANCE
#===============================================
#NAP Phase

nophase <- escalc_data %>%
  filter (!phase %in% c('1','2', 'Pre NAP'))
phase1_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = nophase, method = 'REML')
print(phase1_result)
predict(phase1_result, transf = transf.ilogit)

phase1 <- escalc_data %>%
  filter(phase == '1')
phase1_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = phase1, method = 'REML')
print(phase1_result)
predict(phase1_result, transf = transf.ilogit)

phase2 <- escalc_data %>%
  filter(phase == '2')
phase2_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = phase2, method = 'REML')
print(phase2_result)
predict(phase2_result, transf = transf.ilogit)

phase_prenap <- escalc_data %>%
  filter(phase == 'Pre NAP')
phaseprenap_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = phase_prenap, method = 'REML')
print(phaseprenap_result)
predict(phaseprenap_result, transf = transf.ilogit)


#AST_METHODS

automated <- escalc_data %>%
  filter (ast =='Automated')
automated_result <- rma.mv(yi = yi, V = vi, data = automated, random = ~1 | study_id/effect_id, method = 'REML')
View(automated)
print(automated_result)
predict(automated_result, transf = transf.ilogit)

traditional <- escalc_data %>%
  filter (ast == 'Traditional')
traditional_result <- rma.mv(yi = yi, V = vi, data = traditional, random = ~1 | study_id/effect_id, method = 'REML')
View(traditional)
print(traditional_result)
predict(traditional_result, transf = transf.ilogit)

#Region

SW <- escalc_data %>%
  filter(region == 'South-south' | region == 'South-South' | region == 'South South')
SW_result <- rma.mv(yi = yi, V = vi,  random = ~1 | study_id/effect_id, data = SW, method = 'REML')
print(SW_result)
predict(SW_result, transf = transf.ilogit)


#Facility Type & Ward

Secondary <- escalc_data%>%
  filter(facility_type == 'Tertiary')
Secondary_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Secondary, method = 'REML')
print(Secondary_result)
predict(Secondary_result, transf = transf.ilogit)

Surgical <- escalc_data%>%
  filter(Ward == 'ICU')
Surgical_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Surgical, method = 'REML')
print(Surgical_result)
predict(Surgical_result, transf = transf.ilogit)


#By Specimen and Isolates

Wound_swab <- escalc_data%>%
  filter(specimen_type %in% c('Wound swab', 'Swab','Wound swab, Aspirate'))
Wound_swab_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Wound_swab, method = 'REML')
print(Wound_swab_result)
predict(Wound_swab_result, transf = transf.ilogit)

Urine <- escalc_data%>%
  filter(specimen_type == 'Ear swab' | specimen_type == 'Nasal swab' | specimen_type == 'Throat swab' | specimen_type == 'ENT')
Urine_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Urine, method = 'REML')
print(Urine_result)
predict(Urine_result, transf = transf.ilogit)

Stool <- escalc_data%>%
  filter(specimen_type == 'Stool')
Stool_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Stool, method = 'REML')
print(Stool_result)
predict(Stool_result, transf = transf.ilogit)

Blood <- escalc_data%>%
  filter(specimen_type == 'Blood' | specimen_type == 'Ascitic fluid, Blood')
Blood_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Blood, method = 'REML')
print(Blood_result)
predict(Blood_result, transf = transf.ilogit)

unique(escalc_data$Antibiotic_class)

bac <- escalc_data%>%
  filter(isolates == 'Proteus species' | isolates == 'Proteus mirabilis')
bac_result <- rma(yi = yi, vi = vi, data = bac, method = 'REML')
print(bac_result)
predict(bac_result, transf = transf.ilogit)


#Risk of Bias

rob <- escalc_data%>%
  filter(
    risk_of_bias == 'High'
  )
Kleb_result <- rma.mv(yi = yi, V = vi, random = ~1 | study_id/effect_id, data = Kleb, method = 'REML')
print(Kleb_result)
predict(Kleb_result, transf = transf.ilogit)


#Antibiotic_class tested
unique(amr_data$anti_class)
lactam <- escalc_data%>%
  #filter(anti_class == 'Cephalosporin' | anti_class == 'Penicillin' | anti_class == 'Carbapenem' | anti_class == 'Beta lactam' | anti_class == 'Beta-lactam' | anti_class == 'Monobactam' | anti_class == 'Beta lacta' |anti_class == 'Cephalosporin, Beta lactam' | anti_class == 'Carbapenem, Beta lactam' )
  filter(anti_class == 'Carbapenem' | anti_class == 'Cephalosporin' | anti_class == 'Monobactam' )
lactam_result <- rma.mv(yi = yi, V = vi, data = lactam, random = ~1 | study_id/effect_id, method = 'REML')
print(lactam_result)
predict(lactam_result, transf = transf.ilogit)

unique(amr_data$anti_class)

Sulfonamide <- escalc_data%>%
  filter(anti_class == 'Aminoglycoside' | anti_class == 'Aminoglycosides')
Sulfonamide_result <- rma(yi = yi, vi = vi, data = Sulfonamide, method = 'REML')
print(Sulfonamide_result)
predict(Sulfonamide_result, transf = transf.ilogit)

#Antibiotic

Ampicillin  <- escalc_data%>%
  filter(Antibiotic_class == 'Ofloxacin')
Ampicillin_result <- rma(yi = yi, vi = vi, data = Ampicillin, method = 'REML')
print(Ampicillin_result)
predict(Ampicillin_result, transf = transf.ilogit)

print(amr_data$Antibiotic_class)

Meropenem <- escalc_data%>%
  filter(antibiotic == 'Meropenem')
Meropenem_result <- rma(yi = yi, vi = vi, data = Meropenem, method = 'REML')
print(Meropenem_result)
predict(Meropenem_result, transf = transf.ilogit)

in_patient <- escalc_data%>%
  filter(in_patient_percent == '100')
in_patient_result <- rma(yi = yi, vi = vi, data = in_patient, method = 'REML')
print(in_patient_result)
predict(in_patient_result, transf = transf.ilogit)

table(amr_data$age_group)









#===================
#HIV PRev on AMR Prev
#====================

escalc_datasamplesize<-as.numeric(escalcdatasample_size)
escalc_datahiv <- as.numeric(escalcdatahiv)
hiv_pct <- (escalc_datahiv/escalcdatasample_size) * 100

hiv_result <- rma(
  yi = yi,
  vi = vi,
  mods = ~ hiv_pct,
  data = escalc_data,
  method = "REML"
)

print(hiv_result)

#==============================================
#Exploring Heterogeneity [using Meta-regression]
#================================================

meta_regression_result <- rma.mv(
  yi = yi,
  V = vi,
  mods = ~ last_research_year + sample_size + risk_of_bias_score,
  random = ~1 | study_id/effect_id,
  data = escalc_data,
  method = "REML"
)

print(meta_regression_result)

meta_regression_result <- rma(
  yi = yi,
  vi = vi,
  mods = ~ last_research_year + sample_size + risk_of_bias_score,
  data = escalc_data,
  method = "REML"
)

print(meta_regression_result)


#================================
#SECONDARY/PRIVATE VS TERTIARY
#===============================
sub_amr_data <-
  amr_data %>%
  filter (specimen_type %in% c('Blood', 'Urine', 'Stool', 'ENT', 'Wound swab'), facility_type %in% c('Tertiary', 'Secondary', 'Private'))

table(sub_amr_data$facility_type)
table(sub_amr_data$specimen_type)

chi_specimen_type <- table(
  sub_amr_data$facility_type,
  sub_amr_data$specimen_type
)

chisq.test(chi_specimen_type, simulate.p.value = TRUE,
           B = 10000)
nrow(chi_specimen_type)

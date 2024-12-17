library(tidyverse)
library(here)
library(DirichletReg)
set.seed(111)
N = 1000
p_female = rnorm(N, 0.8, (0.818-0.741)/(2*qnorm(0.975)))
BMI = rnorm(N, 38, 5/qnorm(0.975))
p_ht_bl = rnorm(N, 0.35, 0.1*0.35)
p_smoke = rnorm(N, 0.125, 0.1*0.125)
age = rnorm(N, 45, 0.6/qnorm(0.975))
SBP_BL = rnorm(N, 125, 1/qnorm(0.975))
r = 0.03
dis_BMI = rnorm(N, -0.00333, 0.00222/qnorm(0.975)) 
dis_age = rnorm(N, -0.0007, 0.0002/qnorm(0.975))
HbA1c_BL = rnorm(N, 5.7, (5.85-5.1)/(2*qnorm(0.975)))

CVD_probs = rdirichlet(N, 1000*c(0.55, 0.23, 0.22))
mort_MI = rnorm(N, 0.08, (0.0985-0.052742)/(2*qnorm(0.975)))
MI_HF_young = rnorm(N, 0.0994, 0.1*0.0994)
MI_HF_mid = rnorm(N, 0.1648, 0.1*0.1648)
MI_HF_old = rnorm(N, 0.268, 0.1*0.268)
MI_rec_young = 1 - MI_HF_young - mort_MI
MI_rec_mid = 1 - MI_HF_mid - mort_MI
MI_rec_old = 1 - MI_HF_old - mort_MI

mort_Stroke = rnorm(N, 0.08, (0.103-0.0561)/(2*qnorm(0.975)))
Stroke_rec = 1-mort_Stroke

rr_MI = rgamma(N, 1.58*24, 24)
rr_Stroke = rgamma(N, 3.13*7.25, 7.25)
rr_Other = rgamma(N, 1.9*11, 11)
rr_HF = rgamma(N, 1.82*10, 10)
rr_DM = rgamma(N, 1.15*45000, 45000)

p_recur_stroke = rnorm(N, 0.12, 0.12*0.1) 
p_recur_MI_f = rnorm(N, 0.0723, 0.0723*0.1)
p_recur_MI_m = rnorm(N, 0.0813, 0.0813*0.1)
  
p_toHF_young = rnorm(N, 0.012, 0.012*0.1)
p_toHF_mid = rnorm(N, 0.031, 0.031*0.1)
p_toHF_old = rnorm(N, 0.08, 0.08*0.1)

w_red_SEM = rnorm(N, -0.1375, 0.0125/qnorm(0.975))
w_red_LIR = rnorm(N, -0.05, 0.011/qnorm(0.975))
w_red_PT = rnorm(N, -0.0910, (0.111-0.0715)/(2*qnorm(0.975)))
w_red_BN = rnorm(N, -0.0460, (0.0604-0.0311)/(2*qnorm(0.975)))

h_red_SEM = rnorm(N, -0.3, 0.03/qnorm(0.975))
h_red_LIR = rnorm(N, -0.2, 0.03/qnorm(0.975))
h_red_PT = rnorm(N, 0, 0.0302/qnorm(0.975))
h_red_BN = rnorm(N, 0, 0.0302/qnorm(0.975))

c_Other = rnorm(N, 14279, (18498 - 10056)/(2*qnorm(0.975)))
c_Stroke1 = rgamma(N, 17316*0.001, 0.001)
c_Stroke2 = rgamma(N, 6500*0.018, 0.018)
c_MI1 = rgamma(N, 38788.29*0.00075, 0.00075)
c_MI2 = rgamma(N, 3117*0.035, 0.035)
c_HF1 = rgamma(N, 27030*0.001, 0.001)
c_HF2 = rgamma(N, 15605*0.0015, 0.0015)
c_DM = rgamma(N, 11425*0.0017, 0.0017)

c_SEM = rgamma(N, 13618*0.1, 0.1)
c_LIR_1 = rgamma(N, 11309*0.1, 0.1)
c_LIR_2 = rgamma(N, 11760*0.1, 0.1)
c_PT_1 = rgamma(N, 1355*0.1, 0.1)
c_PT_2 = rgamma(N, 1465*0.1, 0.1)
c_BN_1 = rgamma(N, 2034*0.1, 0.1)
c_BN_2 = rgamma(N, 2095*0.1, 0.1)

disc_LSM = rbeta(N, 0.025*1300, (1-0.025)*1300)
disc_SEM = rbeta(N, 0.042*1310, (1-0.042)*1310)
disc_LIR = rbeta(N, 0.058*2450, (1-0.058)*2450)
disc_PT = rbeta(N, 0.058*510, (1-0.058)*510)
disc_BN = rbeta(N, 0.053*600, (1-0.053)*600)

mort_data = read.csv(here("data", "all_cause_death_probabilities.csv"))
death_prob = function(age, f_prev) {
  age = round(age)
  probs = map2_dbl(age, f_prev, ~ {
    data = mort_data %>% filter(Age == .x)
    return(data$d_prob_f * .y + data$d_prob_m * (1 - .y))
  })
  
  return(probs)
}

trans_prob_matrix = function(age, BMI, SBP, HbA1C, prev_female, prev_smoke, prev_ht, CVD_types, surv_Stroke, surv_MI, MI_to_HF, mort_Stroke, mort_MI, rr_Other, rr_MI, rr_Stroke, rr_HF, rr_DM, rec_Stroke, rec_MI_m, rec_MI_f, toHF) {
  death_p = death_prob(age, prev_female)
  DM_p = 0.00000146 * exp(1.87 * HbA1C) * 0.0197 * exp(0.101*BMI)
  rec_MI = prev_female*rec_MI_f + (1-prev_female)*rec_MI_m
  
  female_risk_CVD = 1 - 0.94833 ^ (exp(2.72107*(log(age)-3.8686) + 0.51125*(log(BMI)-log(28)) + 2.81291*(log(SBP)*c(1,0,1,0) - 4.24) + 2.88267*(log(SBP)*c(0,1,0,1) - 0.5826) + 0.61868*(c(1,1,0,0)-0.3423) + 0.77763*(0-0.0376)))
  CVD_prob_female = sum(c(prev_smoke*(1-prev_ht), prev_smoke*prev_ht, (1-prev_smoke)*(1-prev_ht), (1-prev_smoke)*prev_ht)*(1-exp(-female_risk_CVD/10)))
  male_risk_CVD = 1 - 0.8843 ^ (exp(3.113*(log(age)-3.856) + 0.7928*(log(BMI)-log(28)) + 1.8551*(log(SBP)*c(1,0,1,0) - 4.3544) + 1.9267*(log(SBP)**c(0,1,0,1)  - 0.5019) + 0.7095*(c(1,1,0,0)-0.3522) + 0.5316*(0-0.065)))
  CVD_prob_male = sum(c(prev_smoke*(1-prev_ht), prev_smoke*prev_ht, (1-prev_smoke)*(1-prev_ht), (1-prev_smoke)*prev_ht)*(1-exp(-male_risk_CVD/10)))
  CVD_prob = CVD_prob_female * prev_female + CVD_prob_male * (1-prev_female)
  female_risk_CVD_DM = 1 - 0.94833 ^ (exp(2.72107*(log(age)-3.8686) + 0.51125*(log(BMI)-log(28)) + 2.81291*(log(SBP)*c(1,0,1,0) - 4.24) + 2.88267*(log(SBP)*c(0,1,0,1) - 0.5826) + 0.61868*(c(1,1,0,0)-0.3423) + 0.77763*(1-0.0376)))
  CVD_prob_female_DM = sum(c(prev_smoke*(1-prev_ht), prev_smoke*prev_ht, (1-prev_smoke)*(1-prev_ht), (1-prev_smoke)*prev_ht)*(1-exp(-female_risk_CVD_DM/10)))
  male_risk_CVD_DM = 1 - 0.8843 ^ (exp(3.113*(log(age)-3.856) + 0.7928*(log(BMI)-log(28)) + 1.8551*(log(SBP)*c(1,0,1,0) - 4.3544) + 1.9267*(log(SBP)**c(0,1,0,1)  - 0.5019) + 0.7095*(c(1,1,0,0)-0.3522) + 0.5316*(1-0.065)))
  CVD_prob_male_DM = sum(c(prev_smoke*(1-prev_ht), prev_smoke*prev_ht, (1-prev_smoke)*(1-prev_ht), (1-prev_smoke)*prev_ht)*(1-exp(-male_risk_CVD_DM/10)))
  CVD_prob_DM = CVD_prob_female_DM * prev_female + CVD_prob_male_DM * (1-prev_female)
  
  probs = matrix(c(
    (1-CVD_prob-death_p)*(1-DM_p), CVD_prob*CVD_types[1]*(1-DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(1-DM_p), 0, CVD_prob*CVD_types[3]*surv_MI*(1-DM_p), 0, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(1-DM_p), 0, 0, 0, 0, (1-CVD_prob-death_p)*(DM_p), CVD_prob*CVD_types[1]*(DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(DM_p), 0, CVD_prob*CVD_types[3]*surv_MI*(DM_p), 0, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(DM_p), 0, 0, 0, 0, death_p + CVD_prob*CVD_types[2]*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From Healthy
    0, (1-CVD_prob*(1-CVD_types[1])-death_p*rr_Other)*(1-DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(1-DM_p), 0, CVD_prob*CVD_types[3]*surv_MI*(1-DM_p), 0, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(1-DM_p), 0, 0, 0, 0, 0, (1-CVD_prob*(1-CVD_types[1])-death_p*rr_Other)*(DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(DM_p), 0, CVD_prob*CVD_types[3]*surv_MI*(DM_p), 0, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(DM_p), 0, 0, 0, 0,  death_p*rr_Other + CVD_prob*CVD_types[2]*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From CVD
    0, 0, rec_Stroke*surv_Stroke*(1-DM_p), (1-rec_Stroke-CVD_prob*CVD_types[3]-death_p*rr_Stroke)*(1-DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*surv_MI*(1-DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(1-DM_p), 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), (1-rec_Stroke-CVD_prob*CVD_types[3]-death_p*rr_Stroke)*(DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*surv_MI*(DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(DM_p), 0, 0, death_p*rr_Stroke + rec_Stroke*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From Stroke 1
    0, 0, rec_Stroke*surv_Stroke*(1-DM_p), (1-rec_Stroke-CVD_prob*CVD_types[3]-death_p*rr_Stroke)*(1-DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*surv_MI*(1-DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(1-DM_p), 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), (1-rec_Stroke-CVD_prob*CVD_types[3]-death_p*rr_Stroke)*(DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*surv_MI*(DM_p), 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF*(DM_p), 0, 0, death_p*rr_Stroke + rec_Stroke*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From Stroke 2
    0, 0, 0, 0, rec_MI*surv_MI*(1-DM_p), (1-rec_MI-CVD_prob*CVD_types[2]-toHF-death_p*rr_MI)*(1-DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(1-DM_p), 0, 0, (toHF + CVD_prob*CVD_types[3]*MI_to_HF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, rec_MI*surv_MI*(DM_p), (1-rec_MI-CVD_prob*CVD_types[2]-toHF-death_p*rr_MI)*(DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(DM_p), 0, 0, (toHF + CVD_prob*CVD_types[3]*MI_to_HF)*(DM_p), 0, 0, 0, 0, death_p*rr_MI + rec_MI*mort_MI + CVD_prob*CVD_types[2]*mort_Stroke, # From MI1
    0, 0, 0, 0, rec_MI*surv_MI*(1-DM_p), (1-rec_MI-CVD_prob*CVD_types[2]-toHF-death_p*rr_MI)*(1-DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(1-DM_p), 0, 0, (toHF + CVD_prob*CVD_types[3]*MI_to_HF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, rec_MI*surv_MI*(DM_p), (1-rec_MI-CVD_prob*CVD_types[2]-toHF-death_p*rr_MI)*(DM_p), CVD_prob*CVD_types[2]*surv_Stroke*(DM_p), 0, 0, (toHF + CVD_prob*CVD_types[3]*MI_to_HF)*(DM_p), 0, 0, 0, 0, death_p*rr_MI + rec_MI*mort_MI + CVD_prob*CVD_types[2]*mort_Stroke, # From MI2
    0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(1-DM_p), rec_MI*surv_MI*(1-DM_p), (1-rec_Stroke-rec_MI-toHF-death_p*rr_MI*rr_Stroke)*(1-DM_p), 0, 0, (rec_MI*MI_to_HF + toHF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), rec_MI*surv_MI*(DM_p), (1-rec_Stroke-rec_MI-toHF-death_p*rr_MI*rr_Stroke)*(DM_p), 0, 0, (rec_MI*MI_to_HF + toHF)*(DM_p), 0, 0, death_p*rr_MI*rr_Stroke + rec_Stroke*mort_Stroke + rec_MI*mort_MI, # From Stroke 1 & MI 2
    0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(1-DM_p), rec_MI*surv_MI*(1-DM_p), (1-rec_Stroke-rec_MI-toHF-death_p*rr_MI*rr_Stroke)*(1-DM_p), 0, 0, (rec_MI*MI_to_HF + toHF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), rec_MI*surv_MI*(DM_p), (1-rec_Stroke-rec_MI-toHF-death_p*rr_MI*rr_Stroke)*(DM_p), 0, 0, (rec_MI*MI_to_HF + toHF)*(DM_p), 0, 0, death_p*rr_MI*rr_Stroke + rec_Stroke*mort_Stroke + rec_MI*mort_MI, # From Stroke 2 & MI 1
    0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(1-DM_p), rec_MI*surv_MI*(1-DM_p), (1-rec_Stroke-rec_MI-toHF-death_p*rr_MI*rr_Stroke)*(1-DM_p), 0, 0, (rec_MI*MI_to_HF + toHF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), rec_MI*surv_MI*(DM_p), (1-rec_Stroke-rec_MI-toHF-death_p*rr_MI*rr_Stroke)*(DM_p), 0, 0, (rec_MI*MI_to_HF + toHF)*(DM_p), 0, 0, death_p*rr_MI*rr_Stroke + rec_Stroke*mort_Stroke + rec_MI*mort_MI, # From Stroke 2 & MI 2
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*CVD_types[2]-death_p*rr_HF)*(1-DM_p), 0, CVD_prob*CVD_types[2]*surv_Stroke*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*CVD_types[2]-death_p*rr_HF)*(DM_p), 0, CVD_prob*CVD_types[2]*surv_Stroke*(DM_p), 0, death_p*rr_HF + CVD_prob*CVD_types[2]*mort_Stroke, # From HF 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*CVD_types[2]-death_p*rr_HF)*(1-DM_p), 0, CVD_prob*CVD_types[2]*surv_Stroke*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*CVD_types[2]-death_p*rr_HF)*(DM_p), 0, CVD_prob*CVD_types[2]*surv_Stroke*(DM_p), 0, death_p*rr_HF + CVD_prob*CVD_types[2]*mort_Stroke, # From HF 2
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(1-DM_p), (1-rec_Stroke-death_p*rr_Stroke*rr_HF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), (1-rec_Stroke-death_p*rr_Stroke*rr_HF)*(DM_p), death_p*rr_Stroke*rr_HF + rec_Stroke*mort_Stroke, # From Stroke 2 & HF 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(1-DM_p), (1-rec_Stroke-death_p*rr_Stroke*rr_HF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), (1-rec_Stroke-death_p*rr_Stroke*rr_HF)*(DM_p), death_p*rr_Stroke*rr_HF + rec_Stroke*mort_Stroke, # From Stroke 1 & HF 2
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(1-DM_p), (1-rec_Stroke-death_p*rr_Stroke*rr_HF)*(1-DM_p), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke*(DM_p), (1-rec_Stroke-death_p*rr_Stroke*rr_HF)*(DM_p), death_p*rr_Stroke*rr_HF + rec_Stroke*mort_Stroke, # From Stroke 2 & HF 2
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob-death_p*rr_DM), CVD_prob*CVD_types[1], CVD_prob*CVD_types[2]*surv_Stroke, 0, CVD_prob*CVD_types[3]*surv_MI, 0, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF, 0, 0, 0, 0, death_p*rr_DM + CVD_prob*CVD_types[2]*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*(1-CVD_types[1])-death_p*rr_DM*rr_Other), CVD_prob*CVD_types[2]*surv_Stroke, 0, CVD_prob*CVD_types[3]*surv_MI, 0, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF, 0, 0, 0, 0,  death_p*rr_DM*rr_Other + CVD_prob*CVD_types[2]*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From CVD with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, (1-rec_Stroke-CVD_prob*CVD_types[3]-death_p*rr_DM*rr_Stroke), 0, 0, 0, CVD_prob*CVD_types[3]*surv_MI, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF, 0, 0, death_p*rr_DM*rr_Stroke + rec_Stroke*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From Stroke 1 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, (1-rec_Stroke-CVD_prob*CVD_types[3]-death_p*rr_DM*rr_Stroke), 0, 0, 0, CVD_prob*CVD_types[3]*surv_MI, 0, 0, 0, CVD_prob*CVD_types[3]*MI_to_HF, 0, 0, death_p*rr_DM*rr_Stroke + rec_Stroke*mort_Stroke + CVD_prob*CVD_types[3]*mort_MI, # From Stroke 2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_MI*surv_MI, (1-rec_MI-CVD_prob*CVD_types[2]-toHF-death_p*rr_DM*rr_MI), CVD_prob*CVD_types[2]*surv_Stroke, 0, 0, (toHF + CVD_prob*CVD_types[3]*MI_to_HF), 0, 0, 0, 0, death_p*rr_DM*rr_MI + rec_MI*mort_MI + CVD_prob*CVD_types[2]*mort_Stroke, # From MI1 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_MI*surv_MI, (1-rec_MI-CVD_prob*CVD_types[2]-toHF-death_p*rr_DM*rr_MI), CVD_prob*CVD_types[2]*surv_Stroke, 0, 0, (toHF + CVD_prob*CVD_types[3]*MI_to_HF), 0, 0, 0, 0, death_p*rr_DM*rr_MI + rec_MI*mort_MI + CVD_prob*CVD_types[2]*mort_Stroke, # From MI2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, rec_MI*surv_MI, (1-rec_Stroke-rec_MI-toHF-death_p*rr_DM*rr_MI*rr_Stroke), 0, 0, (rec_MI*MI_to_HF + toHF), 0, 0, death_p*rr_DM*rr_MI*rr_Stroke + rec_Stroke*mort_Stroke + rec_MI*mort_MI, # From Stroke 1 & MI 2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, rec_MI*surv_MI, (1-rec_Stroke-rec_MI-toHF-death_p*rr_DM*rr_MI*rr_Stroke), 0, 0, (rec_MI*MI_to_HF + toHF), 0, 0, death_p*rr_DM*rr_MI*rr_Stroke + rec_Stroke*mort_Stroke + rec_MI*mort_MI, # From Stroke 2 & MI 1 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, rec_MI*surv_MI, (1-rec_Stroke-rec_MI-toHF-death_p*rr_DM*rr_MI*rr_Stroke), 0, 0, (rec_MI*MI_to_HF + toHF), 0, 0, death_p*rr_DM*rr_MI*rr_Stroke + rec_Stroke*mort_Stroke + rec_MI*mort_MI, # From Stroke 2 & MI 2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*CVD_types[2]-death_p*rr_DM*rr_HF), 0, CVD_prob*CVD_types[2]*surv_Stroke, 0, death_p*rr_DM*rr_HF + CVD_prob*CVD_types[2]*mort_Stroke, # From HF 1 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1-CVD_prob*CVD_types[2]-death_p*rr_DM*rr_HF), 0, CVD_prob*CVD_types[2]*surv_Stroke, 0, death_p*rr_DM*rr_HF + CVD_prob*CVD_types[2]*mort_Stroke, # From HF 2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, (1-rec_Stroke-death_p*rr_DM*rr_Stroke*rr_HF), death_p*rr_DM*rr_Stroke*rr_HF + rec_Stroke*mort_Stroke, # From Stroke 2 & HF 1 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, (1-rec_Stroke-death_p*rr_DM*rr_Stroke*rr_HF), death_p*rr_DM*rr_Stroke*rr_HF + rec_Stroke*mort_Stroke, # From Stroke 1 & HF 2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rec_Stroke*surv_Stroke, (1-rec_Stroke-death_p*rr_DM*rr_Stroke*rr_HF), death_p*rr_DM*rr_Stroke*rr_HF + rec_Stroke*mort_Stroke, # From Stroke 2 & HF 2 with DM
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  # From death
  ), byrow = F, nrow = 29)
  
  if(!all(rowSums(probs) == 1)) {
    error("Probabilities don't sum to 1")
  }
  else {
    return(probs)
  }
}

qol_array = function(age, BMI, dis_age, dis_BMI) {
  return(c(
    0.9442 - dis_age*age - dis_BMI*(BMI-21.75), # Healthy
    0.959*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # CVD
    0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.19, # Stroke 1
    0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2
    0.955*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.15, # MI 1
    0.955*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # MI 2
    0.955*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.19, # Stroke 1 & MI 2
    0.955*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.15, # Stroke 2 & MI 1
    0.955*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 & MI 2
    0.93*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # HF 1
    0.93*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # HF 2
    0.93*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 & HF 1
    0.93*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.19, # Stroke 1 & HF 2
    0.93*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 & HF 2
    0.962*0.9442 - dis_age*age - dis_BMI*(BMI-21.75), # DM
    0.962*0.959*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # CVD with DM
    0.962*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.19, # Stroke 1 with DM
    0.962*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 with DM
    0.962*0.955*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.15, # MI 1 with DM
    0.962*0.955*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # MI 2 with DM
    0.962*0.955*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.19, # Stroke 1 & MI 2 with DM
    0.962*0.955*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.15, # Stroke 2 & MI 1 with DM
    0.962*0.955*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 & MI 2 with DM
    0.962*0.93*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # HF 1 with DM
    0.962*0.93*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # HF 2 with DM
    0.962*0.93*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 & HF 1 with DM
    0.962*0.93*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)) - 0.5*0.19, # Stroke 1 & HF 2 with DM
    0.962*0.93*0.943*(0.9442 - dis_age*age - dis_BMI*(BMI-21.75)), # Stroke 2 & HF 2 with DM
    0 # Death
  ))
}

cost_comorbidity_array = function(c_other, c_Stroke1, c_Stroke2, c_MI1, c_MI2, c_HF1, c_HF2, c_DM) {
  return(c(
    0, # Healthy
    c_other, # CVD
    c_Stroke1+c_Stroke2, # Stroke 1
    c_Stroke2, # Stroke 2
    c_MI1+c_MI2, # MI 1
    c_MI2, # MI 2
    c_Stroke1+c_Stroke2+c_MI2, # Stroke 1 & MI 2
    c_Stroke2+c_MI1+c_MI2, #Stroke 2 & MI 1
    c_Stroke2+c_MI2, #Stroke 2 & MI 2
    c_HF1, # HF 1
    c_HF2, # HF 2
    c_Stroke2+c_HF1, # Stroke 2 & HF 1
    c_Stroke1+c_Stroke2+c_HF2, # Stroke 1 & HF 2
    c_Stroke2+c_HF2, # Stroke 2 & HF 2
    c_DM, # DM
    c_DM+c_other, # CVD with DM
    c_DM+c_Stroke1+c_Stroke2, # Stroke 1 with DM
    c_DM+c_Stroke2, # Stroke 2 with DM
    c_DM+c_MI1+c_MI2, # MI 1 with DM
    c_DM+c_MI2, # MI 2 with DM
    c_DM+c_Stroke1+c_Stroke2+c_MI2, # Stroke 1 & MI 2 with DM
    c_DM+c_Stroke2+c_MI1+c_MI2, #Stroke 2 & MI 1 with DM
    c_DM+c_Stroke2+c_MI2, #Stroke 2 & MI 2 with DM
    c_DM+c_HF1, # HF 1 with DM
    c_DM+c_HF2, # HF 2 with DM
    c_DM+c_Stroke2+c_HF1, # Stroke 2 & HF 1 with DM
    c_DM+c_Stroke1+c_Stroke2+c_HF2, # Stroke 1 & HF 2 with DM
    c_DM+c_Stroke2+c_HF2, # Stroke 2 & HF 2 with DM
    0 # Death
  ))
}

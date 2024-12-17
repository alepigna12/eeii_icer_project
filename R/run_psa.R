library(tidyverse)
library(here)
set.seed(111)
N = 1000
p_female = rnorm(N, 0.8, (0.818-0.741)/(2*qnorm(0.975)))
BMI = rnorm(N, 38, 5/qnorm(0.975))
p_ht_bl = 0.35
p_smoke = 0.125
age = rnorm(N, 45, 0.6/qnorm(0.975))
SBP_BL = rnorm(N, 125, 1/qnorm(0.975))
r = 0.03
dis_BMI = rnorm(N, -0.00333, 0.00222/qnorm(0.975)) 
dis_age = rnorm(N, -0.0007, 0.0002/qnorm(0.975))
HbA1c_BL = rnorm(N, 5.7, (5.85-5.1)/(2*qnorm(0.975)))

CVD_probs = c(0.22, 0.23, 0.55)
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

mort_data = read.csv(here("data", "all_cause_death_probabilities.csv"))
death_prob = function(age, f_prev) {
  age = round(age)
  probs = map2_dbl(age, f_prev, ~ {
    data = mort_data %>% filter(Age == .x)
    return(data$d_prob_f * .y + data$d_prob_m * (1 - .y))
  })
  
  return(probs)
}

trans_prob_matrix = function(age, BMI, SBP, HbA1C, prev_female, prev_smoke, prev_ht) {
  death_p = death_prob(age, prev_female)
  DM_p = 0.00000146 * exp(1.87 * HbA1C) * 0.0197 * exp(0.101*BMI)
  prev_ht = 0.35
  age=45
  BMI=38
  SBP=125
  female_risk_CVD = 1 - 0.94833 ^ (exp(2.72107*(log(age)-3.8686) + 0.51125*(log(BMI)-log(28)) + 2.81291*(log(SBP)*c(1,0,1,0) - 4.24) + 2.88267*(log(SBP)*c(0,1,0,1) - 0.5826) + 0.61868*(c(1,1,0,0)-0.3423) + 0.77763*(0-0.0376)))
  CVD_prob_female = sum(c(prev_smoke*(1-prev_ht), prev_smoke*prev_ht, (1-prev_smoke)*(1-prev_ht), (1-prev_smoke)*prev_ht)*(1-exp(-female_risk_CVD/10)))
  male_risk_CVD = 1 - 0.8843 ^ (exp(3.113*(log(age)-3.856) + 0.7928*(log(BMI)-log(28)) + 1.8551*(log(SBP)*c(1,0,1,0) - 4.3544) + 1.9267*(log(SBP)**c(0,1,0,1)  - 0.5019) + 0.7095*(c(1,1,0,0)-0.3522) + 0.5316*(0-0.065)))
  CVD_prob_male = sum(c(prev_smoke*(1-prev_ht), prev_smoke*prev_ht, (1-prev_smoke)*(1-prev_ht), (1-prev_smoke)*prev_ht)*(1-exp(-male_risk_CVD/10)))
  CVD_prob = CVD_prob_female * prev_female + CVD_prob_male * (1-prev_female)
}

library(tidyverse)
library(here)
library(kableExtra)
library(webshot2)
library(gridExtra)
results = readRDS(here("results", "psa.rds"))

ce_50 = 100*c(
  mean((results$SEM_QALY_disc - results$LSM_QALY_disc)*50000 - (results$SEM_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$LIR_QALY_disc - results$LSM_QALY_disc)*50000 - (results$LIR_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$PT_QALY_disc - results$LSM_QALY_disc)*50000 - (results$PT_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$BN_QALY_disc - results$LSM_QALY_disc)*50000 - (results$BN_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$LIR_QALY_disc)*50000 - (results$SEM_cost_disc-results$LIR_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$PT_QALY_disc)*50000 - (results$SEM_cost_disc-results$PT_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$BN_QALY_disc)*50000 - (results$SEM_cost_disc-results$BN_cost_disc) > 0)
)

ce_100 = 100*c(
  mean((results$SEM_QALY_disc - results$LSM_QALY_disc)*100000 - (results$SEM_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$LIR_QALY_disc - results$LSM_QALY_disc)*100000 - (results$LIR_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$PT_QALY_disc - results$LSM_QALY_disc)*100000 - (results$PT_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$BN_QALY_disc - results$LSM_QALY_disc)*100000 - (results$BN_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$LIR_QALY_disc)*100000 - (results$SEM_cost_disc-results$LIR_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$PT_QALY_disc)*100000 - (results$SEM_cost_disc-results$PT_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$BN_QALY_disc)*100000 - (results$SEM_cost_disc-results$BN_cost_disc) > 0)
)

ce_150 = 100*c(
  mean((results$SEM_QALY_disc - results$LSM_QALY_disc)*150000 - (results$SEM_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$LIR_QALY_disc - results$LSM_QALY_disc)*150000 - (results$LIR_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$PT_QALY_disc - results$LSM_QALY_disc)*150000 - (results$PT_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$BN_QALY_disc - results$LSM_QALY_disc)*150000 - (results$BN_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$LIR_QALY_disc)*150000 - (results$SEM_cost_disc-results$LIR_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$PT_QALY_disc)*150000 - (results$SEM_cost_disc-results$PT_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$BN_QALY_disc)*150000 - (results$SEM_cost_disc-results$BN_cost_disc) > 0)
)

ce_200 = 100*c(
  mean((results$SEM_QALY_disc - results$LSM_QALY_disc)*200000 - (results$SEM_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$LIR_QALY_disc - results$LSM_QALY_disc)*200000 - (results$LIR_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$PT_QALY_disc - results$LSM_QALY_disc)*200000 - (results$PT_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$BN_QALY_disc - results$LSM_QALY_disc)*200000 - (results$BN_cost_disc-results$LSM_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$LIR_QALY_disc)*200000 - (results$SEM_cost_disc-results$LIR_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$PT_QALY_disc)*200000 - (results$SEM_cost_disc-results$PT_cost_disc) > 0),
  mean((results$SEM_QALY_disc - results$BN_QALY_disc)*200000 - (results$SEM_cost_disc-results$BN_cost_disc) > 0)
)

trt = c("Semaglutide", "Liraglutide", "Phentermine/Topiramate", "Bupropion/Naltrexone", rep("Semaglutide",3))
comp = c(rep("Lifestyle modification", 4), "Liraglutide", "Phentermine/Topiramate", "Bupropion/Naltrexone")
table4.7 = as.data.frame(cbind(trt, comp, ce_50, ce_100, ce_150, ce_200)) %>%
  kable(col.names = c("Treatment", "Comparator", "Cost Effective at $50,000 per QALY Gained (%)", "Cost Effective at $100,000 per QALY Gained (%)", "Cost Effective at $150,000 per QALY Gained (%)", "Cost Effective at $200,000 per QALY Gained (%)")) %>%
  kable_styling()

save_kable(table4.7, here("results", "table4_7.html"))
webshot(here("results/table4_7.html"), file = here("results", "table4_7.png"))

wtp = seq(0, 300000, 100)
results$inc_cost_SEM = results$SEM_cost_disc-results$LSM_cost_disc
results$inc_QALY_SEM = results$SEM_QALY_disc-results$LSM_QALY_disc
acc_SEM = sapply(wtp, function(w) mean(w * results$inc_QALY_SEM - results$inc_cost_SEM > 0))
results$inc_cost_LIR = results$LIR_cost_disc-results$LSM_cost_disc
results$inc_QALY_LIR = results$LIR_QALY_disc-results$LSM_QALY_disc
acc_LIR = sapply(wtp, function(w) mean(w * results$inc_QALY_LIR - results$inc_cost_LIR > 0))
results$inc_cost_PT = results$PT_cost_disc-results$LSM_cost_disc
results$inc_QALY_PT = results$PT_QALY_disc-results$LSM_QALY_disc
acc_PT = sapply(wtp, function(w) mean(w * results$inc_QALY_PT - results$inc_cost_PT > 0))
results$inc_cost_BN = results$BN_cost_disc-results$LSM_cost_disc
results$inc_QALY_BN = results$BN_QALY_disc-results$LSM_QALY_disc
acc_BN = sapply(wtp, function(w) mean(w * results$inc_QALY_BN - results$inc_cost_BN > 0))

plot1 = ggplot(data = results, aes(x = inc_QALY_SEM, y = inc_cost_SEM)) +
  geom_point(color="blue") +
  labs(
    title = "Incremental Cost-Effectiveness Scatterplot",         
    x = "Incremental Effectiveness (QALY gained)",  
    y = "Incremental Cost (USD)"   
  ) +
  xlim(-.5, 2) +  
  ylim(0, 400000) 

plot2 = ggplot(data = as.data.frame(cbind(wtp, acc_SEM)), aes(x = wtp, y = acc_SEM)) +
  geom_line(color="blue") +
  labs(
    title = "Acceptability Curve",         
    x = "Cost-Effectiveness Threshold (USD/QALY)",  
    y = "Probability Cost Effective"   
  ) +
  xlim(0, 300000) +  
  ylim(0, 1) 

fig = grid.arrange(plot1, plot2, ncol = 2, top = "Figure E4. Results of Probabilistic Sensitivity Analysis: Semaglutide v. Lifestyle Modification")

ggsave(here("results", "figureE4_SEM.png"), fig, width = 10, height = 6)

plot1 = ggplot(data = results, aes(x = inc_QALY_LIR, y = inc_cost_LIR)) +
  geom_point(color="blue") +
  labs(
    title = "Incremental Cost-Effectiveness Scatterplot",         
    x = "Incremental Effectiveness (QALY gained)",  
    y = "Incremental Cost (USD)"   
  ) +
  xlim(-.2, 1) +  
  ylim(0, 350000) 

plot2 = ggplot(data = as.data.frame(cbind(wtp, acc_LIR)), aes(x = wtp, y = acc_LIR)) +
  geom_line(color="blue") +
  labs(
    title = "Acceptability Curve",         
    x = "Cost-Effectiveness Threshold (USD/QALY)",  
    y = "Probability Cost Effective"   
  ) +
  xlim(0, 300000) +  
  ylim(0, 1) 

fig = grid.arrange(plot1, plot2, ncol = 2, top = "Figure E4. Results of Probabilistic Sensitivity Analysis: Liraglutide v. Lifestyle Modification")

ggsave(here("results", "figureE4_LIR.png"), fig, width = 10, height = 6)

plot1 = ggplot(data = results, aes(x = inc_QALY_PT, y = inc_cost_PT)) +
  geom_point(color="blue") +
  labs(
    title = "Incremental Cost-Effectiveness Scatterplot",         
    x = "Incremental Effectiveness (QALY gained)",  
    y = "Incremental Cost (USD)"   
  ) +
  xlim(-.5, 1.25) +  
  ylim(-30000, 40000) 

plot2 = ggplot(data = as.data.frame(cbind(wtp, acc_PT)), aes(x = wtp, y = acc_PT)) +
  geom_line(color="blue") +
  labs(
    title = "Acceptability Curve",         
    x = "Cost-Effectiveness Threshold (USD/QALY)",  
    y = "Probability Cost Effective"   
  ) +
  xlim(0, 300000) +  
  ylim(0, 1) 

fig = grid.arrange(plot1, plot2, ncol = 2, top = "Figure E4. Results of Probabilistic Sensitivity Analysis: Phentermine/Topiramate v. Lifestyle Modification")

ggsave(here("results", "figureE4_PT.png"), fig, width = 10, height = 6)

plot1 = ggplot(data = results, aes(x = inc_QALY_BN, y = inc_cost_BN)) +
  geom_point(color="blue") +
  labs(
    title = "Incremental Cost-Effectiveness Scatterplot",         
    x = "Incremental Effectiveness (QALY gained)",  
    y = "Incremental Cost (USD)"   
  ) +
  xlim(-0.75, 0.8) +  
  ylim(0, 60000) 

plot2 = ggplot(data = as.data.frame(cbind(wtp, acc_BN)), aes(x = wtp, y = acc_BN)) +
  geom_line(color="blue") +
  labs(
    title = "Acceptability Curve",         
    x = "Cost-Effectiveness Threshold (USD/QALY)",  
    y = "Probability Cost Effective"   
  ) +
  xlim(0, 300000) +  
  ylim(0, 1) 

fig = grid.arrange(plot1, plot2, ncol = 2, top = "Figure E4. Results of Probabilistic Sensitivity Analysis: Naltrexone/Bupropion v. Lifestyle Modification")

ggsave(here("results", "figureE4_BN.png"), fig, width = 10, height = 6)

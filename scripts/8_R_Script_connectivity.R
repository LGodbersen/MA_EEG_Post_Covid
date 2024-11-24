## R Script for connectivity analyses
# overview
# 1. load packages
# 2. load data from MATLAB .csv table
# 3. modify dataset
  # 3.1 age match
  # 3.2 add a normalized gcc and cpl
# 4. have a look at the number of epochs
# 5. have a look at the values (SW, GCC, CPL, FC) per threshold 
# 6. Is the threshold itself different?
# 7. remove outliers
# 8. check requirements for icoh
# 9. permutation tests
# 10. compare normal, overall, random
# 11. boxplots
# 12. Correlations with behavioral data
  # 12.1 TMTA, B-A
  # 12.2 FACIT
# 13. corr of sw and coherence
# 14. exploratory: gender differences (you need the R Script power for this as well)

#--------- 1. load packages------------
library(tidyverse)
library(car)
library(readr)
library(skimr)
library(hrbrthemes)
library(viridis)
library(ggExtra)
library(ggdist)
library(ggsignif)
library(ggpubr)
#-------- 2. load data from MATLAB .csv table-----------
# load table with graph measures: these are 4s epochs and were 0.1 high pass filtered
graph_measures_roi <- read_csv("data/analysis_connectivity_icoh/01/table_graph_measures_ROI.csv")
graph_measures_roi_random <- read_csv("data/analysis_connectivity_icoh/01/table_graph_measures_ROI_random.csv")
graph_measures_roi_overall <- read_csv("data/analysis_connectivity_icoh/01/table_graph_measures_ROI_overall.csv")

number_of_epochs_con <- read_csv("data/analysis_connectivity_icoh/01/number_of_epochs_01.csv")

graph_measures_roi <- merge(graph_measures_roi,number_of_epochs_con)
graph_measures_roi_random <- merge(graph_measures_roi_random,number_of_epochs_con)
graph_measures_roi_overall <- merge(graph_measures_roi_overall,number_of_epochs_con)

# ------ 3. modify dataset ---------
# 3.1  age match the participants
age_match <- read.delim("C:/Users/Lara Godbersen/Documents/GitHub/Masters-thesis/data/PuG/matched_participants_conn.tsv",sep="\t")

graph_measures_roi_agematch <- graph_measures_roi%>%
  filter(graph_measures_roi$participant_id %in% age_match$participant_id)

graph_measures_roi_overall_agematch <- graph_measures_roi_overall%>%
  filter(graph_measures_roi_overall$participant_id %in% age_match$participant_id)

graph_measures_roi_random_agematch <- graph_measures_roi_random%>%
  filter(graph_measures_roi_random$participant_id %in% age_match$participant_id)

# modify table (f.ex. add tmt b-a)
graph_measures_roi_agematch <- graph_measures_roi_agematch%>%
  mutate(facit_f_FS = as.numeric(facit_f_FS),
         tmt_b_minus_a = tmt_b_time-tmt_a_time)

graph_measures_roi_random_agematch <- graph_measures_roi_random_agematch%>%
  mutate(facit_f_FS = as.numeric(facit_f_FS),
         tmt_b_minus_a = tmt_b_time-tmt_a_time)

graph_measures_roi_overall_agematch <- graph_measures_roi_overall_agematch%>%
  mutate(facit_f_FS = as.numeric(facit_f_FS),
         tmt_b_minus_a = tmt_b_time-tmt_a_time)


# 3.2 add a normalized gcc and cpl
graph_measures_roi_agematch<- graph_measures_roi_agematch%>%
  mutate(gcc_delta_norm_01 = gcc_delta_01/gcc_rand_delta_01,
         gcc_delta_norm_02 = gcc_delta_02/gcc_rand_delta_02,
         gcc_delta_norm_03 = gcc_delta_03/gcc_rand_delta_03,
         gcc_delta_norm_04 = gcc_delta_04/gcc_rand_delta_04,
         gcc_delta_norm_05 = gcc_delta_05/gcc_rand_delta_05,
         gcc_delta_norm_06 = gcc_delta_06/gcc_rand_delta_06,
         gcc_delta_norm_07 = gcc_delta_07/gcc_rand_delta_07,
         gcc_delta_norm_08 = gcc_delta_08/gcc_rand_delta_08,
         gcc_delta_norm_09 = gcc_delta_09/gcc_rand_delta_09,
         cpl_delta_norm_01 = cpl_delta_01/cpl_rand_delta_01,
         cpl_delta_norm_02 = cpl_delta_02/cpl_rand_delta_02,
         cpl_delta_norm_03 = cpl_delta_03/cpl_rand_delta_03,
         cpl_delta_norm_04 = cpl_delta_04/cpl_rand_delta_04,
         cpl_delta_norm_05 = cpl_delta_05/cpl_rand_delta_05,
         cpl_delta_norm_06 = cpl_delta_06/cpl_rand_delta_06,
         cpl_delta_norm_07 = cpl_delta_07/cpl_rand_delta_07,
         cpl_delta_norm_08 = cpl_delta_08/cpl_rand_delta_08,
         cpl_delta_norm_09 = cpl_delta_09/cpl_rand_delta_09,
         gcc_beta_norm_01 = gcc_beta_01/gcc_rand_beta_01,
         gcc_beta_norm_02 = gcc_beta_02/gcc_rand_beta_02,
         gcc_beta_norm_03 = gcc_beta_03/gcc_rand_beta_03,
         gcc_beta_norm_04 = gcc_beta_04/gcc_rand_beta_04,
         gcc_beta_norm_05 = gcc_beta_05/gcc_rand_beta_05,
         gcc_beta_norm_06 = gcc_beta_06/gcc_rand_beta_06,
         gcc_beta_norm_07 = gcc_beta_07/gcc_rand_beta_07,
         gcc_beta_norm_08 = gcc_beta_08/gcc_rand_beta_08,
         gcc_beta_norm_09 = gcc_beta_09/gcc_rand_beta_09,
         cpl_beta_norm_01 = cpl_beta_01/cpl_rand_beta_01,
         cpl_beta_norm_02 = cpl_beta_02/cpl_rand_beta_02,
         cpl_beta_norm_03 = cpl_beta_03/cpl_rand_beta_03,
         cpl_beta_norm_04 = cpl_beta_04/cpl_rand_beta_04,
         cpl_beta_norm_05 = cpl_beta_05/cpl_rand_beta_05,
         cpl_beta_norm_06 = cpl_beta_06/cpl_rand_beta_06,
         cpl_beta_norm_07 = cpl_beta_07/cpl_rand_beta_07,
         cpl_beta_norm_08 = cpl_beta_08/cpl_rand_beta_08,
         cpl_beta_norm_09 = cpl_beta_09/cpl_rand_beta_09)

graph_measures_roi_random_agematch<- graph_measures_roi_random_agematch%>%
  mutate(gcc_delta_norm_01 = gcc_delta_01/gcc_rand_delta_01,
         gcc_delta_norm_02 = gcc_delta_02/gcc_rand_delta_02,
         gcc_delta_norm_03 = gcc_delta_03/gcc_rand_delta_03,
         gcc_delta_norm_04 = gcc_delta_04/gcc_rand_delta_04,
         gcc_delta_norm_05 = gcc_delta_05/gcc_rand_delta_05,
         gcc_delta_norm_06 = gcc_delta_06/gcc_rand_delta_06,
         gcc_delta_norm_07 = gcc_delta_07/gcc_rand_delta_07,
         gcc_delta_norm_08 = gcc_delta_08/gcc_rand_delta_08,
         gcc_delta_norm_09 = gcc_delta_09/gcc_rand_delta_09,
         cpl_delta_norm_01 = cpl_delta_01/cpl_rand_delta_01,
         cpl_delta_norm_02 = cpl_delta_02/cpl_rand_delta_02,
         cpl_delta_norm_03 = cpl_delta_03/cpl_rand_delta_03,
         cpl_delta_norm_04 = cpl_delta_04/cpl_rand_delta_04,
         cpl_delta_norm_05 = cpl_delta_05/cpl_rand_delta_05,
         cpl_delta_norm_06 = cpl_delta_06/cpl_rand_delta_06,
         cpl_delta_norm_07 = cpl_delta_07/cpl_rand_delta_07,
         cpl_delta_norm_08 = cpl_delta_08/cpl_rand_delta_08,
         cpl_delta_norm_09 = cpl_delta_09/cpl_rand_delta_09,
         gcc_beta_norm_01 = gcc_beta_01/gcc_rand_beta_01,
         gcc_beta_norm_02 = gcc_beta_02/gcc_rand_beta_02,
         gcc_beta_norm_03 = gcc_beta_03/gcc_rand_beta_03,
         gcc_beta_norm_04 = gcc_beta_04/gcc_rand_beta_04,
         gcc_beta_norm_05 = gcc_beta_05/gcc_rand_beta_05,
         gcc_beta_norm_06 = gcc_beta_06/gcc_rand_beta_06,
         gcc_beta_norm_07 = gcc_beta_07/gcc_rand_beta_07,
         gcc_beta_norm_08 = gcc_beta_08/gcc_rand_beta_08,
         gcc_beta_norm_09 = gcc_beta_09/gcc_rand_beta_09,
         cpl_beta_norm_01 = cpl_beta_01/cpl_rand_beta_01,
         cpl_beta_norm_02 = cpl_beta_02/cpl_rand_beta_02,
         cpl_beta_norm_03 = cpl_beta_03/cpl_rand_beta_03,
         cpl_beta_norm_04 = cpl_beta_04/cpl_rand_beta_04,
         cpl_beta_norm_05 = cpl_beta_05/cpl_rand_beta_05,
         cpl_beta_norm_06 = cpl_beta_06/cpl_rand_beta_06,
         cpl_beta_norm_07 = cpl_beta_07/cpl_rand_beta_07,
         cpl_beta_norm_08 = cpl_beta_08/cpl_rand_beta_08,
         cpl_beta_norm_09 = cpl_beta_09/cpl_rand_beta_09)

graph_measures_roi_overall_agematch<- graph_measures_roi_overall_agematch%>%
  mutate(gcc_delta_norm_01 = gcc_delta_01/gcc_rand_delta_01,
         gcc_delta_norm_02 = gcc_delta_02/gcc_rand_delta_02,
         gcc_delta_norm_03 = gcc_delta_03/gcc_rand_delta_03,
         gcc_delta_norm_04 = gcc_delta_04/gcc_rand_delta_04,
         gcc_delta_norm_05 = gcc_delta_05/gcc_rand_delta_05,
         gcc_delta_norm_06 = gcc_delta_06/gcc_rand_delta_06,
         gcc_delta_norm_07 = gcc_delta_07/gcc_rand_delta_07,
         gcc_delta_norm_08 = gcc_delta_08/gcc_rand_delta_08,
         gcc_delta_norm_09 = gcc_delta_09/gcc_rand_delta_09,
         cpl_delta_norm_01 = cpl_delta_01/cpl_rand_delta_01,
         cpl_delta_norm_02 = cpl_delta_02/cpl_rand_delta_02,
         cpl_delta_norm_03 = cpl_delta_03/cpl_rand_delta_03,
         cpl_delta_norm_04 = cpl_delta_04/cpl_rand_delta_04,
         cpl_delta_norm_05 = cpl_delta_05/cpl_rand_delta_05,
         cpl_delta_norm_06 = cpl_delta_06/cpl_rand_delta_06,
         cpl_delta_norm_07 = cpl_delta_07/cpl_rand_delta_07,
         cpl_delta_norm_08 = cpl_delta_08/cpl_rand_delta_08,
         cpl_delta_norm_09 = cpl_delta_09/cpl_rand_delta_09,
         gcc_beta_norm_01 = gcc_beta_01/gcc_rand_beta_01,
         gcc_beta_norm_02 = gcc_beta_02/gcc_rand_beta_02,
         gcc_beta_norm_03 = gcc_beta_03/gcc_rand_beta_03,
         gcc_beta_norm_04 = gcc_beta_04/gcc_rand_beta_04,
         gcc_beta_norm_05 = gcc_beta_05/gcc_rand_beta_05,
         gcc_beta_norm_06 = gcc_beta_06/gcc_rand_beta_06,
         gcc_beta_norm_07 = gcc_beta_07/gcc_rand_beta_07,
         gcc_beta_norm_08 = gcc_beta_08/gcc_rand_beta_08,
         gcc_beta_norm_09 = gcc_beta_09/gcc_rand_beta_09,
         cpl_beta_norm_01 = cpl_beta_01/cpl_rand_beta_01,
         cpl_beta_norm_02 = cpl_beta_02/cpl_rand_beta_02,
         cpl_beta_norm_03 = cpl_beta_03/cpl_rand_beta_03,
         cpl_beta_norm_04 = cpl_beta_04/cpl_rand_beta_04,
         cpl_beta_norm_05 = cpl_beta_05/cpl_rand_beta_05,
         cpl_beta_norm_06 = cpl_beta_06/cpl_rand_beta_06,
         cpl_beta_norm_07 = cpl_beta_07/cpl_rand_beta_07,
         cpl_beta_norm_08 = cpl_beta_08/cpl_rand_beta_08,
         cpl_beta_norm_09 = cpl_beta_09/cpl_rand_beta_09)

# ------- 4. have a look at number of epochs -----------------
# number of epochs
graph_measures_roi%>%
  group_by(group)%>%
  summarise(mean_epoch = mean(number_epochs),
            sd_epoch = sd(number_epochs))

t.test(number_epochs~group, data = graph_measures_roi, alternative = "two.sided", paired = FALSE)# 0.4701
#effsize
graph_measures_roi%>%
  ungroup()%>%
  cohens_d(number_epochs~group)# small

graph_measures_roi%>%
  group_by(group)%>%
  summarise(max_epoch = max(number_epochs),
            min_epoch = min(number_epochs))


# do they correlate with the power?
graph_measures_roi%>%
  ggplot(aes(x = number_epochs, y = coh_delta, color = group))+
  geom_point()
cor.test(graph_measures_roi$coh_delta,graph_measures_roi$number_epochs)

graph_measures_roi%>%
  ggplot(aes(x = number_epochs, y = coh_beta, color = group))+
  geom_point()
cor.test(graph_measures_roi$coh_beta,graph_measures_roi$number_epochs)

# ------- 5. have a look at the distribution per threshold ------------
# modify the table -> longer format

# small worldness delta
test_threshold <- graph_measures_roi_agematch[, 1:258]
test_threshold <- test_threshold%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("sw_delta_"),
               names_to = "threshold_sw", 
               values_to = "sw_delta")

test_threshold%>%
  ggplot(aes(x = threshold_sw, y= sw_delta, fill = threshold_sw))+
  #geom_violin(width=1.4, show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey", show.legend = FALSE)

# gcc delta
test_threshold <- graph_measures_roi_agematch[, 236:276]
test_threshold <- test_threshold%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("gcc_delta_norm_"),
               names_to = "threshold_gcc", 
               values_to = "gcc_delta")

test_threshold%>%
  ggplot(aes(x = threshold_gcc, y= gcc_delta, fill = threshold_gcc))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey",show.legend = FALSE)

# cpl delta
test_threshold <- graph_measures_roi_agematch[, 236:296]
test_threshold <- test_threshold%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("cpl_delta_norm_"),
               names_to = "threshold_cpl", 
               values_to = "cpl_delta")

test_threshold%>%
  ggplot(aes(x = threshold_cpl, y= cpl_delta, fill = threshold_cpl))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey", show.legend = FALSE)

# small worldness beta
test_threshold <- graph_measures_roi_agematch[, 1:258]
test_threshold <- test_threshold%>%
  #select(-sw_beta_01)%>%
  pivot_longer(cols = starts_with("sw_beta_"),
               names_to = "threshold_sw", 
               values_to = "sw_beta")

test_threshold%>%
  ggplot(aes(x = threshold_sw, y= sw_beta, fill = group))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey", show.legend = FALSE)

# gcc beta
test_threshold <- graph_measures_roi_agematch[, 1:294]
test_threshold <- test_threshold%>%
  #select(-gcc_beta_norm_01)%>%
  pivot_longer(cols = starts_with("gcc_beta_norm_"),
               names_to = "threshold_gcc", 
               values_to = "gcc_beta")

test_threshold%>%
  ggplot(aes(x = threshold_gcc, y= gcc_beta, fill = threshold_gcc))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey",show.legend = FALSE)

# cpl beta
test_threshold <- graph_measures_roi_agematch[, 236:296]
test_threshold <- test_threshold%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("cpl_beta_norm_"),
               names_to = "threshold_cpl", 
               values_to = "cpl_beta")

test_threshold%>%
  ggplot(aes(x = threshold_cpl, y= cpl_beta, fill = threshold_cpl))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey",show.legend = FALSE)

# now for the coherence
test_coherence <- graph_measures_roi_agematch%>%
  mutate(mean_coh_d = mean(coh_delta),
         sd_coh_d = sd(coh_delta),
         lower_bound = mean_coh_d - 3 * sd_coh_d,
         upper_bound = mean_coh_d + 3 * sd_coh_d) %>%
  filter(coh_delta >= lower_bound & coh_delta <= upper_bound) %>%
  ungroup() 

test_coherence <- test_coherence%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("coh_delta"),
               names_to = "threshold_sw", 
               values_to = "coh_delta_compair")

test_coherence%>%
  ggplot(aes(x = threshold_sw, y= coh_delta_compair, fill = group))+
  #geom_violin(width=1.4, show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey", show.legend = FALSE)

# beta
test_coherence <- graph_measures_roi_agematch%>%
  mutate(mean_coh_b = mean(coh_beta),
         sd_coh_b = sd(coh_beta),
         lower_bound = mean_coh_b - 3 * sd_coh_b,
         upper_bound = mean_coh_b + 3 * sd_coh_b) %>%
  filter(coh_beta >= lower_bound & coh_beta <= upper_bound) %>%
  ungroup() 

test_coherence <- test_coherence%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("coh_beta"),
               names_to = "threshold_sw", 
               values_to = "coh_beta_compair")

test_coherence%>%
  ggplot(aes(x = threshold_sw, y= coh_beta_compair, fill = group))+
  #geom_violin(width=1.4, show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey",show.legend = FALSE)

# ------ 6. Is the threshold itself different? -------
test_threshold <- graph_measures_roi_agematch%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("thr_delta"),
               names_to = "threshold", 
               values_to = "thr_delta")

test_threshold%>%
  ggplot(aes(x = threshold, y= thr_delta, fill = group))+
  #geom_violin(width=1.4, show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey", show.legend = FALSE)

# just an example t-test: you would probably need permutation tests to test it for all of them
t.test(thr_delta_03~group,data = graph_measures_roi_agematch)


# for beta
test_threshold <- graph_measures_roi_agematch%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("thr_beta"),
               names_to = "threshold", 
               values_to = "thr_beta")

test_threshold%>%
  ggplot(aes(x = threshold, y= thr_beta, fill = group))+
  #geom_violin(width=1.4, show.legend = FALSE) +
  geom_boxplot(width=0.1, color="grey", show.legend = FALSE)

# just an example t-test: you would probably need permutation tests to test it for all of them
t.test(thr_beta_03~group,data = graph_measures_roi_agematch)

# ------------------ 7. remove outliers ------------------------------
# coh delta
graph_measures_final <-graph_measures_roi_agematch%>%
  group_by(group) %>%
  mutate(mean_coh_delta = mean(coh_delta),
         sd_coh_delta = sd(coh_delta),
         lower_bound = mean_coh_delta - 3 * sd_coh_delta,
         upper_bound = mean_coh_delta + 3 * sd_coh_delta) %>%
  filter(coh_delta >= lower_bound & coh_delta <= upper_bound) %>%
  ungroup()

# coh beta
graph_measures_final <-graph_measures_final%>%
  group_by(group) %>%
  mutate(mean_coh_beta = mean(coh_beta),
         sd_coh_beta = sd(coh_beta),
         lower_bound = mean_coh_beta - 3 * sd_coh_beta,
         upper_bound = mean_coh_beta + 3 * sd_coh_beta) %>%
  filter(coh_beta >= lower_bound & coh_beta <= upper_bound) %>%
  ungroup()

# coh delta overall
graph_measures_final_overall <-graph_measures_roi_overall_agematch%>%
  group_by(group) %>%
  mutate(mean_coh_delta = mean(coh_delta),
         sd_coh_delta = sd(coh_delta),
         lower_bound = mean_coh_delta - 3 * sd_coh_delta,
         upper_bound = mean_coh_delta + 3 * sd_coh_delta) %>%
  filter(coh_delta >= lower_bound & coh_delta <= upper_bound) %>%
  ungroup()

# coh beta overall
graph_measures_final_overall <-graph_measures_final_overall%>%
  group_by(group) %>%
  mutate(mean_coh_beta = mean(coh_beta),
         sd_coh_beta = sd(coh_beta),
         lower_bound = mean_coh_beta - 3 * sd_coh_beta,
         upper_bound = mean_coh_beta + 3 * sd_coh_beta) %>%
  filter(coh_beta >= lower_bound & coh_beta <= upper_bound) %>%
  ungroup()

# coh delta random
graph_measures_final_random <-graph_measures_roi_random_agematch%>%
  group_by(group) %>%
  mutate(mean_coh_delta = mean(coh_delta),
         sd_coh_delta = sd(coh_delta),
         lower_bound = mean_coh_delta - 3 * sd_coh_delta,
         upper_bound = mean_coh_delta + 3 * sd_coh_delta) %>%
  filter(coh_delta >= lower_bound & coh_delta <= upper_bound) %>%
  ungroup()

# coh beta random
graph_measures_final_random <-graph_measures_final_random%>%
  group_by(group) %>%
  mutate(mean_coh_beta = mean(coh_beta),
         sd_coh_beta = sd(coh_beta),
         lower_bound = mean_coh_beta - 3 * sd_coh_beta,
         upper_bound = mean_coh_beta + 3 * sd_coh_beta) %>%
  filter(coh_beta >= lower_bound & coh_beta <= upper_bound) %>%
  ungroup()

#----------------- 8. check requirements for the coherence --------------------------------
shapiro_df_withPCS <- graph_measures_final%>%
  filter(group == 'withPCS')

shapiro_df_withoutPCS <- graph_measures_final%>%
  filter(group == 'withoutPCS')

# coh delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = coh_delta))+
  geom_histogram(color = "black", fill = "white",bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()

shapiro.test(shapiro_df_withPCS$coh_delta)
shapiro.test(shapiro_df_withoutPCS$coh_delta)
leveneTest(coh_delta~group,data = graph_measures_final)

graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = coh_beta))+
  geom_histogram(color = "black", fill = "white",bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# the distributions do look different

shapiro.test(shapiro_df_withPCS$coh_beta)
shapiro.test(shapiro_df_withoutPCS$coh_beta)
leveneTest(coh_beta~group,data = graph_measures_final)# t test bei beta also okay

# ---------- 9. permutation tests -------------------------
# comparing the small worldness over all thresholds
# to find out if one of the thresholds is worth exploring any further

# Create permutation test function (mit ttest für beta)
permutation_function_con <- function(nsim,df,freq,x){
  res <- numeric(nsim) ## set aside space for results
  sw <- paste('sw_', freq, '_0', x,sep = '')
  df <- df%>%select(group,!!sym(sw))
  for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df))
    bdat <- transform(df,group = group[perm])
    ## compute & store difference in means; store the value
    tt <- t.test(bdat[[sw]]~group,data=bdat,var.equal=FALSE)
    res[i] <- tt$statistic
  }
  obs <- t.test(df[[sw]]~group,data=df,var.equal=FALSE)
  ## append the observed value to the list of results
  res <- c(res,obs$statistic)
  res <<- res
  obs <<- obs$statistic
}

visualize_results <- function(res,obs){
  result_sim <- sum(res > obs)/(length(res)-1)
  Histogramm_pre <- hist(res,
                         ylab = 'Anzahl',
                         xlab = substitute(paste(italic('t'), ' values')),
                         main = substitute(paste(italic('t'), ' values with observed value')))
  Histogramm <- abline(v=obs,col="red")
  return(result_sim)
  return(Histogramm)
}

# Create permutation test function for wilcox test
permutation_function_con_wilcox <- function(nsim,df,freq,x){
  res <- numeric(nsim) ## set aside space for results
  sw <- paste('sw_', freq, '_0', x,sep = '')
  df <- df%>%select(group,!!sym(sw))
  for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df))
    bdat <- transform(df,group = group[perm])
    ## compute & store difference in means; store the value
    w <- wilcox.test(bdat[[sw]]~group,data=bdat,var.equal=FALSE)
    res[i] <- w$statistic
  }
  obs <- wilcox.test(df[[sw]]~group,data=df,var.equal=FALSE)
  ## append the observed value to the list of results
  res <- c(res,obs$statistic)
  res <<- res
  obs <<- obs$statistic
}


# first delta
set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',1)
delta_1 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',2)
delta_2 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',3)
delta_3 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',4)
delta_4 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',5)
delta_5 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',6)
delta_6 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',7)
delta_7 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',8)
delta_8 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'delta',9)
delta_9 <- visualize_results(res,obs)

# beta
set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',1)
beta_1 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',2)
beta_2 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',3)
beta_3 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',4)
beta_4 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',5)
beta_5 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',6)
beta_6 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',7)
beta_7 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',8)
beta_8 <- visualize_results(res,obs)

set.seed(1)
permutation_function_con(1000,graph_measures_final,'beta',9)
beta_9 <- visualize_results(res,obs)

# visualise it
color_palette <- c("with PCS" = '#F59541',
                   "without PCS" = "#02CAF5")

# small worldness beta
test_threshold_final <- graph_measures_final%>%
  #select(-sw_beta_01)%>%
  pivot_longer(cols = starts_with("sw_beta_"),
               names_to = "threshold_sw", 
               values_to = "sw_beta")

test_threshold_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  mutate(threshold_sw = fct_recode(threshold_sw,
                                   "10" = "sw_beta_01",
                                   "20" = "sw_beta_02",
                                   "30" = "sw_beta_03",
                                   "40" = "sw_beta_04",
                                   "50" = "sw_beta_05",
                                   "60" = "sw_beta_06",
                                   "70" = "sw_beta_07",
                                   "80" = "sw_beta_08",
                                   "90" = "sw_beta_09"))%>%
  ggplot(aes(x = threshold_sw, y= sw_beta, colour = group))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.35)+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  labs(x = 'Relative Threshold [in %]',
       y = 'Small World Index beta')


test_threshold_final_delta <- graph_measures_final%>%
  #select(-sw_beta_01)%>%
  pivot_longer(cols = starts_with("sw_delta_"),
               names_to = "threshold_sw", 
               values_to = "sw_delta")

test_threshold_final_delta%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  mutate(threshold_sw = fct_recode(threshold_sw,
                                   "10" = "sw_delta_01",
                                   "20" = "sw_delta_02",
                                   "30" = "sw_delta_03",
                                   "40" = "sw_delta_04",
                                   "50" = "sw_delta_05",
                                   "60" = "sw_delta_06",
                                   "70" = "sw_delta_07",
                                   "80" = "sw_delta_08",
                                   "90" = "sw_delta_09"))%>%
  ggplot(aes(x = threshold_sw, y= sw_delta, colour = group))+
  #geom_violin(width=1.4,show.legend = FALSE) +
  geom_boxplot(width=0.35)+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  labs(x = 'Relative Threshold [in %]',
       y = 'Small World Index delta')


graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_03))+
  geom_histogram(color = "black", fill = "white",bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# the distributions do look different

shapiro.test(shapiro_df_withPCS$sw_beta_03)
shapiro.test(shapiro_df_withoutPCS$sw_beta_03)
leveneTest(sw_beta_03~group,data = graph_measures_final)# t test bei beta also okay


#---------- 10. compare the methods -----------------------------------
graph_measures_final <- graph_measures_final%>%
  mutate(method = 'normal')

graph_measures_final_overall <- graph_measures_final_overall%>%
  mutate(method = 'overall')

graph_measures_final_random <- graph_measures_final_random%>%
  mutate(method = 'random')

combined_table <- rbind(graph_measures_final,graph_measures_final_overall, graph_measures_final_random)

combined_table%>%
  group_by(group, method)%>%
  ggplot(aes(x = method, y = coh_delta, color = group))+
  geom_boxplot(outlier.color = 'black')+
  geom_jitter(alpha = 0.5,position = position_jitterdodge())


combined_table%>%
  group_by(group, method)%>%
  ggplot(aes(x = method, y = coh_beta, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.1, alpha = 0.5)

# put a significant threshold from before here
# delta
combined_table%>%
  group_by(group, method)%>%
  ggplot(aes(x = method, y = sw_delta_05, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.1, alpha = 0.5)

# beta
combined_table%>%
  group_by(group, method)%>%
  ggplot(aes(x = method, y = sw_beta_05, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.1, alpha = 0.5)

# check requirements for SW
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_03))+
  geom_histogram(color = "black", fill = "white",bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()

shapiro.test(shapiro_df_withPCS$sw_beta_03)
shapiro.test(shapiro_df_withoutPCS$sw_beta_03)
leveneTest(sw_beta_03~group,data = graph_measures_final)


graph_measures_final%>%
  group_by(group)%>%
  summarise(mean_age = mean(age),
            sd_age = sd(age),
            mean_epochs = mean(number_epochs))# age of course still good, but withoutPCS has more epochs which is a problem!!


# ------------------11. boxplots ---------------------
sum_data_conn<- graph_measures_final%>%
  group_by(group)%>%
  summarise(mean_coh_d = mean(coh_delta),
            sd_coh_d = sd(coh_delta),
            mean_coh_b = mean(coh_beta),
            sd_coh_b = sd(coh_beta),
            mean_sw_d10 = mean(sw_delta_01),
            sd_sw_d10 = sd(sw_delta_01),
            mean_sw_d20 = mean(sw_delta_02),
            sd_sw_d20 = sd(sw_delta_02),
            mean_sw_d30 = mean(sw_delta_03),
            sd_sw_d30 = sd(sw_delta_03),
            mean_sw_d40 = mean(sw_delta_04),
            sd_sw_d40 = sd(sw_delta_04),
            mean_sw_d50 = mean(sw_delta_05),
            sd_sw_d50 = sd(sw_delta_05),
            mean_sw_d60 = mean(sw_delta_06),
            sd_sw_d60 = sd(sw_delta_06),
            mean_sw_d70 = mean(sw_delta_07),
            sd_sw_d70 = sd(sw_delta_07),
            mean_sw_d80 = mean(sw_delta_08),
            sd_sw_d80 = sd(sw_delta_08),
            mean_sw_d90 = mean(sw_delta_09),
            sd_sw_d90 = sd(sw_delta_09),
            mean_sw_b10 = mean(sw_beta_01),
            sd_sw_b10 = sd(sw_beta_01),
            mean_sw_b20 = mean(sw_beta_02),
            sd_sw_b20 = sd(sw_beta_02),
            mean_sw_b30 = mean(sw_beta_03),
            sd_sw_b30 = sd(sw_beta_03),
            mean_sw_b40 = mean(sw_beta_04),
            sd_sw_b40 = sd(sw_beta_04),
            mean_sw_b50 = mean(sw_beta_05),
            sd_sw_b50 = sd(sw_beta_05),
            mean_sw_b60 = mean(sw_beta_06),
            sd_sw_b60 = sd(sw_beta_06),
            mean_sw_b70 = mean(sw_beta_07),
            sd_sw_b70 = sd(sw_beta_07),
            mean_sw_b80 = mean(sw_beta_08),
            sd_sw_b80 = sd(sw_beta_08),
            mean_sw_b90 = mean(sw_beta_09),
            sd_sw_b90 = sd(sw_beta_09)
            )

# coherence delta
graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = coh_delta, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = function(p) sprintf("p = %.2g", p), test = 'wilcox.test', color = 'black')+
  labs(y = 'mean delta imag(coh)')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )

g4 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"),
         group = fct_relevel(group, "with PCS", "without PCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = coh_delta, color = sex))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2)+
  labs(y = 'mean delta imag(coh)')+
  theme_classic()+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )


wilcox.test(coh_delta~group, data = graph_measures_final, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)
t.test(coh_delta~group, data = graph_measures_final, alternative = "two.sided", paired = FALSE)

graph_measures_final <- graph_measures_final%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = graph_measures_final,coh_delta~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')

# get the effsize
graph_measures_final%>%
  ungroup()%>% # apparently you have to ungroup here, otherwise, wilcox_effsize does not work
  wilcox_effsize(coh_delta~group)


# coherence beta
graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = coh_beta, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = function(p) sprintf("p = %.2g", p), test = 't.test', color = 'black')+
  labs(y = 'mean beta imag(coh)')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )# plot does not have the one sided p value yet

t.test(coh_beta~group, data = graph_measures_final, alternative = "less", paired = FALSE)
wilcox.test(coh_beta~group, data = graph_measures_final, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)
cohens_d(graph_measures_final,coh_beta ~ group)


g3 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"),
         group = fct_relevel(group, "with PCS", "without PCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = coh_beta, color = sex))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2)+
  labs(y = 'mean beta imag(coh)')+
  theme_classic()+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )# plot does not have the one sided p value yet


#small worldness delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = sw_delta_05, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5)


t.test(sw_delta_05~group, data = graph_measures_final, alternative = "two.sided", paired = FALSE)


# small worldness beta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = sw_beta_05, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5)

t.test(sw_beta_05~group, data = graph_measures_final, alternative = "two.sided", paired = FALSE)

g5 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"),
         group = fct_relevel(group, "with PCS", "without PCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = sw_beta_03, color = sex))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2)+
  labs(y = 'mean beta SWI (30%)')+
  theme_classic()+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )# plot does not have the one sided p value yet


g6 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"),
         group = fct_relevel(group, "with PCS", "without PCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = sw_delta_03, color = sex))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2)+
  labs(y = 'mean delta SWI (30%)')+
  theme_classic()+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )# plot does not have the one sided p value yet

# gcc_delta z.b. 05
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = gcc_delta_norm_05, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5)


# gcc_beta z.b. 04
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = gcc_beta_norm_04, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5)

# cpl delta z.b. 05
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = cpl_delta_norm_05, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5)


# cpl norm beta z.b. 03
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = cpl_beta_norm_03, color = group))+
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5)

graph_measures_final%>%
  group_by(group)%>%
  summarise(mean_cpl = mean(cpl_beta_03),
            sd_cpl = sd(cpl_beta_03),
            mean_gcc = mean(gcc_beta_03),
            sd_gcc = sd(gcc_beta_03))

# -------- 12. Correlations with behavioral data ---------
# beta and delta
cor.test(graph_measures_final$coh_beta,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$sw_beta_03,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$sw_beta_03,graph_measures_final$coh_beta, method = 'spearman', exact = FALSE)

## -------------- 12.1 TMT-A-----------------------
# coherence delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = coh_delta,y = tmt_a_time, color = group))+
  stat_cor(aes(color = "Correlation: "),method = "pearson", label.x = 0.055, label.y = 40,hjust=0)+
  geom_point()

cor.test(graph_measures_final$tmt_a_time,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)


# coherence beta
p1 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = coh_beta,y = tmt_a_time, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic()+
  stat_cor(aes(color = "Correlation: "),method = "spearman", p.accuracy = .001, r.accuracy = 0.01, label.x = 0.025, label.y = 45,hjust=0)+
  theme(legend.position = c(0.17, 0.85))+
  labs(x = 'rel. central beta imag(coh)',
       y = 'TMT A [in seconds, 0-200]')+
  theme(text = element_text(size = 16)) 

ggMarginal(p1, type = "densigram")

cor.test(graph_measures_final$tmt_a_time,graph_measures_final$coh_beta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$coh_beta, method = 'spearman', exact = FALSE)

# gcc_ delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = gcc_delta_norm_05,y = tmt_a_time, color = group))+
  geom_point()


# gcc beta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = gcc_beta_norm_04,y = tmt_a_time, color = group))+
  geom_point()


# cpl delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = cpl_delta_norm_05,y = tmt_a_time, color = group))+
  geom_point()

cor.test(graph_measures_roi_agematch$cpl_delta_norm_05,graph_measures_roi_agematch$tmt_a_time)


# cpl beta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = cpl_beta_norm_04,y = tmt_a_time, color = group))+
  geom_point()


# small worldness delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_delta_05,y = tmt_a_time, color = group))+
  geom_point()

# small worldness beta
p2 <- graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_03,y = tmt_a_time, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme(legend.position = c(0.15, 0.85))

ggMarginal(p2, type = "densigram")


# tmta
# beta
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_beta_09, method = 'spearman', exact = FALSE)
# delta
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_a_time,graph_measures_final$sw_delta_09, method = 'spearman', exact = FALSE)

# tmt b-a
# beta
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_beta_09, method = 'spearman', exact = FALSE)
# delta
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$tmt_b_minus_a,graph_measures_final$sw_delta_09, method = 'spearman', exact = FALSE)


##--------------- 12.2 FACIT-F --------------------
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = coh_delta,y = facit_f_FS, color = group))+
  geom_point()

cor.test(graph_measures_final$facit_f_FS,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)


# coherence beta
p4 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = coh_beta,y = facit_f_FS, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic()+
  stat_cor(aes(color = "Correlation: "),method = "spearman", p.accuracy = .001, r.accuracy = 0.01, label.x = 0.023, label.y = 43,hjust=0)+
  theme(legend.position = c(0.165, 0.16))+
  labs(x = 'rel. central beta imag(coh)',
       y = 'FACIT Fatigue Score [0-50]')+
  theme(text = element_text(size = 16))  # Adjust the size here

ggMarginal(p4, type = "densigram")

cor.test(graph_measures_final$facit_f_FS,graph_measures_final$coh_beta, method = 'spearman', exact = FALSE)

# gcc_ delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = gcc_delta_norm_05,y = facit_f_FS, color = group))+
  geom_point()

# gcc beta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = gcc_beta_norm_04,y = facit_f_FS, color = group))+
  geom_point()


# cpl delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = cpl_delta_norm_05,y = facit_f_FS, color = group))+
  geom_point()

# cpl beta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = cpl_beta_norm_04,y = facit_f_FS, color = group))+
  geom_point()

# small worldness delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_delta_05,y = facit_f_FS, color = group))+
  geom_point()

cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_delta_09, method = 'spearman', exact = FALSE)


# small worldness beta
p5 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_03,y = facit_f_FS, color = group))+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic()+
  stat_cor(aes(color = "Correlation: "),method = "pearson", label.x = 0.98, label.y = 10,hjust=0)+
  geom_point()+
  theme(legend.position = c(0.10, 0.27))+
  labs(x = 'Small World Index bei 40 % threshold',
       y = 'FACIT Fatigue Score [0-50]')

ggMarginal(p5, type = "densigram")

cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)

p10 <- graph_measures_final%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_03,y = facit_f_FS, color = group))+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic()+
  stat_cor(aes(color = "Correlation: "),method = "pearson", label.x = 0.99, label.y = 10,hjust=0)+
  geom_point()+
  theme(legend.position = c(0.12, 0.27))+
  labs(x = 'Small World Index bei 30 % threshold',
       y = 'FACIT Fatigue Score [0-50]')

ggMarginal(p10, type = "densigram")

cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$facit_f_FS,graph_measures_final$sw_beta_09, method = 'spearman', exact = FALSE)



# 12. 2 HADS-D -----------------
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$coh_beta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)

# beta sw
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_beta_09, method = 'spearman', exact = FALSE)
# delta sw
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$hads_d_total_score,graph_measures_final$sw_delta_09, method = 'spearman', exact = FALSE)


# 12.3 MoCA --------------------
cor.test(graph_measures_final$moca,graph_measures_final$coh_delta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$coh_beta, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)

# beta 
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_beta_09, method = 'spearman', exact = FALSE)
# delta 
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_01, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_02, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_03, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_04, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_05, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_06, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_07, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_08, method = 'spearman', exact = FALSE)
cor.test(graph_measures_final$moca,graph_measures_final$sw_delta_09, method = 'spearman', exact = FALSE)


# -------- 13. corr of sw and coherence ---------------
# delta
p3 <- graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_delta_05, y = coh_delta, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme(legend.position = c(0.15, 0.85))

ggMarginal(p3, type = "densigram")

cor.test(graph_measures_final$sw_delta_05,graph_measures_final$coh_delta)

# beta
p4 <- graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_04, y = coh_beta, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme(legend.position = c(0.15, 0.85))

ggMarginal(p4, type = "densigram")

cor.test(graph_measures_final$sw_beta_04,graph_measures_final$coh_beta)


# corr of sw and gcc/cpl
# delta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_delta_05, y = gcc_delta_norm_05, color = group))+
  geom_point()

graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_delta_05, y = cpl_delta_norm_05, color = group))+
  geom_point()

cor.test(graph_measures_final$sw_delta_05,graph_measures_final$cpl_delta_norm_05)

# beta
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_04, y = gcc_beta_norm_04, color = group))+
  geom_point()

graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = sw_beta_04, y = cpl_beta_norm_04, color = group))+
  geom_point()

# gcc and cpl
graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = cpl_delta_norm_05, y = gcc_delta_norm_05, color = group))+
  geom_point()

cor.test(graph_measures_final$gcc_delta_norm_05,graph_measures_final$cpl_delta_norm_05)

graph_measures_final%>%
  group_by(group)%>%
  ggplot(aes(x = cpl_beta_norm_09, y = gcc_beta_norm_09, color = group))+
  geom_point()

cor.test(graph_measures_final$gcc_delta_norm_02,graph_measures_final$cpl_delta_norm_02)


# correlation of small worldness and functional connectivity per threshold
graph_measures_final <- graph_measures_final%>%
  mutate(cor_sw_fc_01 = cor(coh_delta_01,sw_delta_01),
         cor_sw_fc_02 = cor(coh_delta_02,sw_delta_02),
         cor_sw_fc_03 = cor(coh_delta_03,sw_delta_03),
         cor_sw_fc_04 = cor(coh_delta_04,sw_delta_04),
         cor_sw_fc_05 = cor(coh_delta_05,sw_delta_05),
         cor_sw_fc_06 = cor(coh_delta_06,sw_delta_06),
         cor_sw_fc_07 = cor(coh_delta_07,sw_delta_07),
         cor_sw_fc_08 = cor(coh_delta_08,sw_delta_08),
         cor_sw_fc_09 = cor(coh_delta_09,sw_delta_09))

test_cor_sw_fc <- graph_measures_final%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("cor_sw_fc_"),
               names_to = "threshold", 
               values_to = "cor_sw_fc")

test_cor_sw_fc%>%
  ggplot(aes(x = threshold, y= cor_sw_fc, fill = threshold))+
  geom_point()

graph_measures_final%>%
  group_by(group)%>%
  mutate(cor_sw_fc_01 = cor(coh_delta_01,sw_delta_01),
         cor_sw_fc_02 = cor(coh_delta_02,sw_delta_02),
         cor_sw_fc_03 = cor(coh_delta_03,sw_delta_03),
         cor_sw_fc_04 = cor(coh_delta_04,sw_delta_04),
         cor_sw_fc_05 = cor(coh_delta_05,sw_delta_05),
         cor_sw_fc_06 = cor(coh_delta_06,sw_delta_06),
         cor_sw_fc_07 = cor(coh_delta_07,sw_delta_07),
         cor_sw_fc_08 = cor(coh_delta_08,sw_delta_08),
         cor_sw_fc_09 = cor(coh_delta_09,sw_delta_09))%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("cor_sw_fc_"),
               names_to = "threshold", 
               values_to = "cor_sw_fc")%>%
  ggplot(aes(x = threshold, y= cor_sw_fc, color = group))+
  geom_point()

graph_measures_final%>%
  group_by(group)%>%
  mutate(cor_sw_fc_01 = cor(coh_beta_01,sw_beta_01),
         cor_sw_fc_02 = cor(coh_beta_02,sw_beta_02),
         cor_sw_fc_03 = cor(coh_beta_03,sw_beta_03),
         cor_sw_fc_04 = cor(coh_beta_04,sw_beta_04),
         cor_sw_fc_05 = cor(coh_beta_05,sw_beta_05),
         cor_sw_fc_06 = cor(coh_beta_06,sw_beta_06),
         cor_sw_fc_07 = cor(coh_beta_07,sw_beta_07),
         cor_sw_fc_08 = cor(coh_beta_08,sw_beta_08),
         cor_sw_fc_09 = cor(coh_beta_09,sw_beta_09))%>%
  #select(-sw_delta_01)%>%
  pivot_longer(cols = starts_with("cor_sw_fc_"),
               names_to = "threshold", 
               values_to = "cor_sw_fc")%>%
  ggplot(aes(x = threshold, y= cor_sw_fc, color = group))+
  geom_point()

# 14. exploratory: gender differences
grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 2)
grid.arrange(g1, g2, g3, nrow = 3)
grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 1)
grid.arrange(g2, arrangeGrob(g3, g4, nrow=2), nrow = 1)

grid.arrange(arrangeGrob(g1, g2, ncol=2), arrangeGrob(g4, g3, ncol=2),arrangeGrob(g6, g5, ncol=2), nrow = 3)

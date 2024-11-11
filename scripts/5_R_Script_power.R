## R Script for power analysis
# ---------- content -----------------------
# 1. load packages
# 2. load data
# 3. summarise mean
# 4. demographics
# 5. outlier removal
  # 5.1 delta (relative and absolute)
  # 5.2 beta (relative)
  # 5.3 aperiodic components
# 6. export tables for topoplots
# 7. check requirements (normality, variances, etc.)
# 8. boxplots and stats
  # 8.1 aperiodic exponent (whole brain)
  # 8.2 aperiodic offset (whole brain)
  # 8.3 rel and abs delta frontal
  # 8.4 rel beta central
  # 8.5 tables of all EEG values
# 9. plot behavioral data and corr tests
  # 9.1 just behavioral data
  # 9.2 corr tests with behav - EEG data
    # 9.2.1 rel delta w TMTA & B-A
    # 9.2.2 rel delta w moca
    # 9.2.3 rel/abs delta w FACIT
    # 9.2.4 rel delta w hads
    # 9.2.5 rel beta w TMTA & B-A
    # 9.2.6 rel beta w FACIT
    # 9.2.7 aperiodic exponent with everything
    # 9.2.8 aperiodic offset with everything
# 10. r squared
# 11. permutation tests

#------------ 1. load packages------------------
library(tidyverse)
library(car)
library(readr)
library(ggdist)
library(ggExtra)# displaying distributions next to plots
library(gridExtra)
library(ggsignif)# displaying stats in plots
library(ggpubr)
library(coin)# need this for z value of wilcox test
library(effsize)# for cohens d
library(rstatix)# for wilcox test
library(dplyr)

#--------------- 2. load data--------------------
# load csv file that I created in MATLAB (has ID, channel, aperiodic offset, aperiodic exponent, abs and rel delta and beta power)
table_power_5 <- read_csv("data/analysis_power/table_power_final_01.csv") # this is the 5s data set with a 0.1 high pass filtering
number_of_epochs_5 <- read_csv("data/analysis_power/number_of_epochs_01.csv")# and load the number of 'good' epochs
number_of_bad_channels <- read_csv("data/analysis_power/number_of_bad_channels.csv")
table_power_5 <- merge(table_power_5, number_of_epochs_5)# put them together
table_power_5 <- merge(table_power_5, number_of_bad_channels)# put them together


table_power_final_01_janka <- read_csv("data/analysis_power/table_power_final_01_janka.csv")
test_table_janka <- table_power_final_01_janka%>%
  group_by(participant_id, group)%>%
  summarise(mean_d = mean(rel_delta),
            mean_b = mean(rel_beta))

# modify table (f.ex. add tmt b-a)
table_power_5 <- table_power_5%>%
  mutate(facit_f_FS = as.numeric(facit_f_FS),
         tmt_b_minus_a = tmt_b_time-tmt_a_time)

# age match the participants
age_match <- read.delim("C:/Users/Lara Godbersen/Documents/GitHub/Masters-thesis/data/PuG/matched_participants_conn.tsv",sep="\t")

table_power_5 <- table_power_5%>%
  filter(table_power_5$participant_id %in% age_match$participant_id)

test_table <- table_power_5%>%
  group_by(participant_id, group)%>%
  summarise(mean_d = mean(rel_delta),
            mean_b = mean(rel_beta))

# Define the channel names you want to select (for delta)
frontal_channels <- c('22','105','11','40','75','39','49','82','48','19','112','25','94','93','83','92','95','96','21','50','10','59','26')

# Filter rows with the specified channel names
table_power_frontal <- table_power_5%>%
  filter(table_power_5$channel %in% frontal_channels)

# Define the channel names you want to select (for beta)
central_channels <- c('85','65','90','66','1','68','3','67','2','70','74','76','81','34','37','42','86','43','87','44','88','45','89','46','77','5','78','6','7','79','8','80','71','35','72','36','73')

# Filter rows with the specified channel names
table_power_central <- table_power_5%>%
  filter(table_power_5$channel %in% central_channels)

#-------3. summarise mean -----------------
df_corr_frontal <- table_power_frontal%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age,years_of_education,sex,hads_d_total_score, number_epochs, moca)%>%
  summarise(mean_delta_power = mean(rel_delta),
            mean_beta_power = mean(rel_beta),
            mean_theta_power = mean(rel_theta),
            mean_alpha_power = mean(rel_alpha),
            mean_aperiodic_exponent = mean(aperiodic_exponent))

df_corr_central <- table_power_central%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,hads_d_total_score, number_epochs,sex)%>%
  summarise(mean_delta_power = mean(rel_delta),
            mean_beta_power = mean(rel_beta),
            mean_theta_power = mean(rel_theta),
            mean_alpha_power = mean(rel_alpha),
            mean_aperiodic_exponent = mean(aperiodic_exponent))

# is the variance different between the groups?
leveneTest(mean_delta_power~group,data = df_corr_frontal)# not significant
leveneTest(mean_beta_power~group,data = df_corr_central)# not significant
#------- 4. demographics-----------------
shapiro_df_withPCS <- df_corr_frontal%>%
  filter(group == 'withPCS')

shapiro_df_withoutPCS <- df_corr_frontal%>%
  filter(group == 'withoutPCS')
# sex
df_corr_frontal%>%
  group_by(group,sex)%>%
  count()# with PCS 16f 7m, without PCS 13f 10m

# age
df_corr_frontal%>%
  group_by(participant_id)%>%
  ggplot(aes(age))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()

df_corr_frontal%>%
  group_by(group)%>%
  summarise(mean_age = mean(age),
            sd_age = sd(age),
            min_age = min(age),
            max_age = max(age))

t.test(age~group, data = df_corr_frontal, alternative = "two.sided")# 0.64
wilcox.test(age~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.72
#effsize
df_corr_frontal%>%
  ungroup()%>%
  wilcox_effsize(age~group)# small

df_corr_frontal <- df_corr_frontal%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_frontal,
                            age~group, 
                            comparisons = list(c('withPCS','withoutPCS')),
                            alternative = 'two.sided')



# years of education
df_corr_frontal%>%
  group_by(participant_id)%>%
  ggplot(aes(years_of_education))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()

t.test(years_of_education~group, data = df_corr_frontal, alternative = "two.sided")# 0.82
wilcox.test(years_of_education~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.85
#effsize
df_corr_frontal%>%
  ungroup()%>%
  wilcox_effsize(years_of_education~group)# small

# FACIT
df_corr_frontal%>%
  group_by(group)%>%
  ggplot(aes(facit_f_FS))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# withoutPCS is very skew

shapiro.test(shapiro_df_withPCS$facit_f_FS)
shapiro.test(shapiro_df_withoutPCS$facit_f_FS)
leveneTest(facit_f_FS~group,data = df_corr_frontal)# not significant
wilcox.test(facit_f_FS~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.0027
#effsize
df_corr_frontal%>%
  ungroup()%>%
  wilcox_effsize(facit_f_FS~group)# moderate

t.test(facit_f_FS~group, data = df_corr_frontal,
       alternative = "two.sided", paired = FALSE)# significant p = 0.0031

df_corr_frontal <- df_corr_frontal%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_frontal,facit_f_FS~group,
                            comparisons = list(c('withPCS','withoutPCS')),
                            alternative = 'two.sided',
                            distribution = 'exact',
                            conf.int = TRUE,
                            conf.level = 0.95)


# HADS
df_corr_frontal%>%
  group_by(group)%>%
  ggplot(aes(hads_d_total_score))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# withoutPCS is very skew, with PCS also a little

shapiro.test(shapiro_df_withPCS$hads_d_total_score)
shapiro.test(shapiro_df_withoutPCS$hads_d_total_score)
leveneTest(hads_d_total_score~group,data = df_corr_frontal)# not significant
t.test(hads_d_total_score~group, data = df_corr_frontal, alternative = "two.sided", paired = FALSE)# p = 0.095
wilcox.test(hads_d_total_score~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.0276
#effsize
df_corr_frontal%>%
  ungroup()%>%
  wilcox_effsize(hads_d_total_score~group)# moderate

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_frontal,hads_d_total_score~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')


# TMT A
df_corr_frontal%>%
  group_by(group)%>%
  ggplot(aes(tmt_a_time))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# looks okay

shapiro.test(shapiro_df_withPCS$tmt_a_time)
shapiro.test(shapiro_df_withoutPCS$tmt_a_time)
leveneTest(tmt_a_time~group,data = df_corr_frontal)# not significant
t.test(tmt_a_time~group, data = df_corr_frontal, alternative = "two.sided", paired = FALSE)# p = 0.059
wilcox.test(tmt_a_time~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.038
#effsize
df_corr_frontal%>%
  ungroup()%>%
  wilcox_effsize(tmt_a_time~group)# moderate

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_frontal,tmt_a_time~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')


# TMT B-A
df_corr_frontal%>%
  group_by(group)%>%
  ggplot(aes(tmt_b_minus_a))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# both skew

shapiro.test(shapiro_df_withPCS$tmt_b_minus_a)
shapiro.test(shapiro_df_withoutPCS$tmt_b_minus_a)
leveneTest(tmt_b_minus_a~group,data = df_corr_frontal)# not significant
t.test(tmt_b_minus_a~group, data = df_corr_frontal, alternative = "two.sided", paired = FALSE)# 0.2434
wilcox.test(tmt_b_minus_a~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.2768
#effsize
df_corr_frontal%>%
  ungroup()%>%
  wilcox_effsize(tmt_b_minus_a~group)# small
# in order to get the z value
result <- coin::wilcox_test(data = df_corr_frontal,tmt_b_minus_a~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')


# MOCA
df_corr_frontal%>%
  group_by(group)%>%
  ggplot(aes(moca))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# both skew

shapiro.test(shapiro_df_withPCS$moca)
shapiro.test(shapiro_df_withoutPCS$moca)
leveneTest(moca~group,data = df_corr_frontal)# not significant
t.test(moca~group, data = df_corr_frontal, alternative = "two.sided", paired = FALSE)
wilcox.test(moca~group, data = df_corr_frontal, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.6025
#effsize
df_corr_frontal%>%
  ungroup()%>%
  cohens_d(moca ~ group)

# number of epochs
df_corr_frontal%>%
  group_by(group)%>%
  summarise(mean_epoch = mean(number_epochs),
            sd_epoch = sd(number_epochs))

t.test(number_epochs~group, data = df_corr_frontal, alternative = "two.sided", paired = FALSE)# 0.4701
#effsize
df_corr_frontal%>%
  ungroup()%>%
  cohens_d(number_epochs~group)# small

df_corr_frontal%>%
  group_by(group)%>%
  summarise(max_epoch = max(number_epochs),
            min_epoch = min(number_epochs))

# do they correlate with the power?
df_corr_frontal%>%
  ggplot(aes(x = number_epochs, y = mean_delta_power, color = group))+
  geom_point()
cor.test(df_corr_frontal$mean_delta_power,df_corr_frontal$number_epochs)

df_corr_central%>%
  ggplot(aes(x = number_epochs, y = mean_beta_power, color = group))+
  geom_point()
cor.test(df_corr_central$mean_delta_power,df_corr_central$number_epochs)

# number of epochs correlates with fatigue score!
table_power_5%>%
  ggplot(aes(x = number_epochs, y = facit_f_FS))+
  geom_point()
cor.test(table_power_5$number_epochs, table_power_5$facit_f_FS)

# summarize values into one table
table_behav <- df_corr_frontal%>%
  group_by(group)%>%
  summarise(mean_facit = mean(facit_f_FS, na.rm = T),
            sd_facit = sd(facit_f_FS, na.rm = T),
            mean_hads = mean(hads_d_total_score, na.rm = T),
            sd_hads = sd(hads_d_total_score, na.rm = T),
            mean_tmta = mean(tmt_a_time),
            sd_tmta = sd(tmt_a_time),
            mean_tmtb_a = mean(tmt_b_minus_a),
            sd_tmtb_a = sd(tmt_b_minus_a),
            mean_y_o = mean(years_of_education),
            sd_y_o = sd(years_of_education),
            mean_epoc = mean(number_epochs),
            sd_epoc = sd(number_epochs),
            mean_moca = mean(moca, na.rm = T),
            sd_moca = sd(moca, na.rm = T))

# good channels
channel_artefacts <- table_power_5%>%
  group_by(group)%>%
  summarise(mean_channels_ica = mean(num_chan_ica),
            sd_channels_ica = sd(num_chan_ica),
            max_channels_ica = max(num_chan_ica),
            min_channels_ica = min(num_chan_ica),
            mean_channels_arte = mean(num_chan_artefact),
            sd_channels_arte = sd(num_chan_artefact),
            max_channels_arte = max(num_chan_artefact),
            min_channels_arte = min(num_chan_artefact))

#------ 5. exclude outliers--------
##-------- 5.1 delta ---------------
# relative delta power frontal
df_corr_frontal%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power, color = group))+
  geom_boxplot()# two outliers in the without PCS group

table_power_5%>%
  group_by(group)%>%
  count()

# outlier removal inside the participant with +/- 3SD
table_delta_filtered <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_rel_delta = mean(rel_delta),
         sd_rel_delta = sd(rel_delta),
         lower_bound = mean_rel_delta - 3 * sd_rel_delta,
         upper_bound = mean_rel_delta + 3 * sd_rel_delta) %>%
  filter(rel_delta >= lower_bound & rel_delta <= upper_bound) %>%
  ungroup()


table_delta_filtered%>%
  group_by(group)%>%
  count()#2830 in with PCS, 2844 in without PCS


# last step: remove negative values
table_delta_filtered$rel_delta <- ifelse(
  table_delta_filtered$rel_delta < 0, 0, table_delta_filtered$rel_delta)

table_delta_filtered%>%
  group_by(group)%>%
  filter(rel_delta == 0)%>%
  count()# 1123 in with PCS, 1082 in without PCS


# visualize the filtered data
table_delta_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = rel_delta, color = group))+
  geom_boxplot(outlier.colour = 'black')+
  geom_jitter()

# select only frontal channels
table_delta_frontal_filtered <- table_delta_filtered%>%
  filter(table_delta_filtered$channel %in% frontal_channels)

df_corr_frontal_filtered <- table_delta_frontal_filtered%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age,years_of_education,moca,sex)%>%
  summarise(mean_delta_power = mean(rel_delta))

df_corr_frontal_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power, color = group))+
  geom_boxplot(outlier.colour = 'black')+
  geom_jitter()# 4 outliers in without PCS

# additional filtering (across group)
table_delta_filtered_group <- table_delta_filtered%>%
  group_by(group)%>%
  mutate(mean_rel_delta = mean(rel_delta),
         sd_rel_delta = sd(rel_delta),
         lower_bound = mean_rel_delta - 3 * sd_rel_delta,
         upper_bound = mean_rel_delta + 3 * sd_rel_delta) %>%
  filter(rel_delta >= lower_bound & rel_delta <= upper_bound) %>%
  ungroup()

table_delta_filtered_group%>%
  group_by(group)%>%
  count()#2800 withPCS, 2785 without PCS

table_frontal_filtered_group <- table_delta_filtered_group%>%
  filter(table_delta_filtered_group$channel %in% frontal_channels)

df_corr_frontal_filtered_group <- table_frontal_filtered_group%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age,moca,hads_d_total_score,sex)%>%
  summarise(mean_delta_power = mean(rel_delta),
            mean_delta_power = mean(rel_delta),
            mean_aperiodic_exponent = mean(aperiodic_exponent))

df_corr_frontal_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power, color = group))+
  geom_boxplot(outlier.colour = 'black')+
  geom_jitter()# 2 outliers in withoutPCS and 1 in withPCS


# outliers delta absolute
table_delta_filtered_abs <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_abs_delta = mean(abs_delta),
         sd_abs_delta = sd(abs_delta),
         lower_bound = mean_abs_delta - 3 * sd_abs_delta,
         upper_bound = mean_abs_delta + 3 * sd_abs_delta) %>%
  filter(abs_delta >= lower_bound & abs_delta <= upper_bound) %>%
  ungroup()

table_delta_frontal_filtered_abs <- table_delta_filtered_abs%>%
  filter(table_delta_filtered_abs$channel %in% frontal_channels)

df_corr_frontal_filtered_abs <- table_delta_frontal_filtered_abs%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age,years_of_education)%>%
  summarise(mean_delta_power_abs = mean(abs_delta))

df_corr_frontal_filtered_abs%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power_abs, color = group))+
  geom_boxplot(outlier.colour = 'black')+
  geom_jitter()

t.test(mean_delta_power_abs~group, data = df_corr_frontal_filtered_abs, alternative = 'less')

## ----------- 5.2 relative beta power central ----------------------------
df_corr_central%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta_power, color = group))+
  geom_boxplot()# one outlier in without PCS, one in the withPCS group

# sd +- 3 for beta
table_beta_filtered <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_rel_beta = mean(rel_beta),
         sd_rel_beta = sd(rel_beta),
         lower_bound = mean_rel_beta - 3 * sd_rel_beta,
         upper_bound = mean_rel_beta + 3 * sd_rel_beta) %>%
  filter(rel_beta >= lower_bound & rel_beta <= upper_bound) %>%
  ungroup()

table_beta_filtered%>%
  group_by(group)%>%
  count()# 2847, 2847 that looks weird..

# last step: remove negative values
table_beta_filtered$rel_beta <- ifelse(
  table_beta_filtered$rel_beta < 0, 0, table_beta_filtered$rel_beta)

table_beta_filtered%>%
  group_by(group)%>%
  filter(rel_beta == 0)%>%
  count()# with PCS 58, without PCS 73.

table_central_filtered <- table_beta_filtered%>%
  filter(table_beta_filtered$channel %in% central_channels)

df_corr_central_filtered <- table_central_filtered%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age, moca,hads_d_total_score,sex)%>%
  summarise(mean_beta_power = mean(rel_beta))

df_corr_central_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta_power, color = group))+
  geom_boxplot(outlier.colour = 'black')+
  geom_jitter()# one outlier in withoutPCS group and one in withPCS group

# additional filtering (across group)
table_beta_filtered_group <- table_beta_filtered%>%
  group_by(group)%>%
  mutate(mean_rel_beta = mean(rel_beta),
         sd_rel_beta = sd(rel_beta),
         lower_bound = mean_rel_beta - 3 * sd_rel_beta,
         upper_bound = mean_rel_beta + 3 * sd_rel_beta) %>%
  filter(rel_beta >= lower_bound & rel_beta <= upper_bound) %>%
  ungroup()

table_beta_filtered_group%>%
  group_by(group)%>%
  count()

table_central_filtered_group <- table_beta_filtered_group%>%
  filter(table_beta_filtered_group$channel %in% central_channels)

df_corr_central_filtered_group <- table_central_filtered_group%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age, moca, hads_d_total_score,sex)%>%
  summarise(mean_beta_power = mean(rel_beta))

df_corr_central_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta_power, color = group))+
  geom_boxplot(outlier.colour = 'black')+
  geom_jitter()# 1 outlier in without PCS and 2 in with PCS

# sd +- 3 for beta 1
table_beta1_filtered <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_rel_beta1 = mean(rel_beta1),
         sd_rel_beta1 = sd(rel_beta1),
         lower_bound = mean_rel_beta1 - 3 * sd_rel_beta1,
         upper_bound = mean_rel_beta1 + 3 * sd_rel_beta1) %>%
  filter(rel_beta1 >= lower_bound & rel_beta1 <= upper_bound) %>%
  ungroup()

# last step: remove negative values
table_beta1_filtered$rel_beta1 <- ifelse(
  table_beta1_filtered$rel_beta1 < 0, 0, table_beta1_filtered$rel_beta1)

table_central1_filtered <- table_beta1_filtered%>%
  filter(table_beta1_filtered$channel %in% central_channels)

df_corr_central1_filtered <- table_central1_filtered%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age, moca,hads_d_total_score)%>%
  summarise(mean_beta1_power = mean(rel_beta1))

df_corr_central1_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta1_power, color = group))+
  geom_boxplot(outlier.colour = 'black')

# sd +- 3 for beta 2
table_beta2_filtered <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_rel_beta2 = mean(rel_beta2),
         sd_rel_beta2 = sd(rel_beta2),
         lower_bound = mean_rel_beta2 - 3 * sd_rel_beta2,
         upper_bound = mean_rel_beta2 + 3 * sd_rel_beta2) %>%
  filter(rel_beta2 >= lower_bound & rel_beta2 <= upper_bound) %>%
  ungroup()

# last step: remove negative values
table_beta2_filtered$rel_beta2 <- ifelse(
  table_beta2_filtered$rel_beta2 < 0, 0, table_beta2_filtered$rel_beta2)

table_central2_filtered <- table_beta2_filtered%>%
  filter(table_beta2_filtered$channel %in% central_channels)

df_corr_central2_filtered <- table_central2_filtered%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age, moca,hads_d_total_score)%>%
  summarise(mean_beta2_power = mean(rel_beta2))

df_corr_central2_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta2_power, color = group))+
  geom_boxplot(outlier.colour = 'black')

wilcox.test(mean_beta2_power~group, data = df_corr_central2_filtered, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)

## ----------- 5.3 aperiodic components -------------------------------------------
# +- 3 sd for aperiodic exponent
table_ape_filtered <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_ape = mean(aperiodic_exponent),
         sd_ape = sd(aperiodic_exponent),
         lower_bound = mean_ape - 3 * sd_ape,
         upper_bound = mean_ape + 3 * sd_ape) %>%
  filter(aperiodic_exponent >= lower_bound & aperiodic_exponent <= upper_bound) %>%
  ungroup()

table_ape_filtered%>%
  group_by(group)%>%
  count()

df_corr_ape <- table_ape_filtered%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age,moca,hads_d_total_score,sex)%>%
  summarise(mean_aperiodic_exponent = mean(aperiodic_exponent))

df_corr_ape%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_aperiodic_exponent, color = group))+
  geom_boxplot()+
  geom_jitter()

# +- 3 SD for aperiodic offset
table_apo_filtered <-table_power_5%>%
  group_by(participant_id) %>%
  mutate(mean_apo = mean(aperiodic_offset),
         sd_apo = sd(aperiodic_offset),
         lower_bound = mean_apo - 3 * sd_apo,
         upper_bound = mean_apo + 3 * sd_apo) %>%
  filter(aperiodic_offset >= lower_bound & aperiodic_offset <= upper_bound) %>%
  ungroup()

table_apo_filtered%>%
  group_by(group)%>%
  count()

df_corr_apo <- table_apo_filtered%>%
  group_by(participant_id,group,tmt_a_time,facit_f_FS, tmt_b_minus_a,age,moca,hads_d_total_score,sex)%>%
  summarise(mean_aperiodic_offset = mean(aperiodic_offset))

df_corr_apo%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_aperiodic_offset, color = group))+
  geom_boxplot()+
  geom_jitter()

#-------- 6. export tables for topoplots ---------------------
# beta power
export_beta_pcs <- table_beta_filtered_group%>%
  filter(group == 'withPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_rel_beta = mean(rel_beta))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

export_beta_c <- table_beta_filtered_group%>%
  filter(group == 'withoutPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_rel_beta = mean(rel_beta))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

# delta power
export_delta_pcs <- table_delta_filtered_group%>%
  filter(group == 'withPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_rel_delta = mean(rel_delta))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

export_delta_c <- table_delta_filtered_group%>%
  filter(group == 'withoutPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_rel_delta = mean(rel_delta))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

# save in folder
write.table(export_beta_pcs, file = "export_beta_pcs.txt", row.names = FALSE, col.names = FALSE)
write.table(export_beta_c, file = "export_beta_c.txt", row.names = FALSE, col.names = FALSE)
write.table(export_delta_pcs, file = "export_delta_pcs.txt", row.names = FALSE, col.names = FALSE)
write.table(export_delta_c, file = "export_delta_c.txt", row.names = FALSE, col.names = FALSE)

# have a look at min/max values for visualisation purposes
export_beta_c%>%
  summarise(min = min(mean_rel_beta),
            max = max(mean_rel_beta))

export_beta_pcs%>%
  summarise(min = min(mean_rel_beta),
            max = max(mean_rel_beta))

export_delta_c%>%
  summarise(min = min(mean_rel_delta),
            max = max(mean_rel_delta))

export_delta_pcs%>%
  summarise(min = min(mean_rel_delta),
            max = max(mean_rel_delta))

# now the same for aperiodic exponent
export_ape_pcs <- table_ape_filtered%>%
  filter(group == 'withPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_ape = mean(aperiodic_exponent))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

export_ape_c <- table_ape_filtered%>%
  filter(group == 'withoutPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_ape = mean(aperiodic_exponent))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

export_apo_pcs <- table_apo_filtered%>%
  filter(group == 'withPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_apo = mean(aperiodic_offset))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

export_apo_c <- table_apo_filtered%>%
  filter(group == 'withoutPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_apo = mean(aperiodic_offset))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

write.table(export_ape_pcs, file = "export_ape_pcs.txt", row.names = FALSE, col.names = FALSE)
write.table(export_ape_c, file = "export_ape_c.txt", row.names = FALSE, col.names = FALSE)
write.table(export_apo_pcs, file = "export_apo_pcs.txt", row.names = FALSE, col.names = FALSE)
write.table(export_apo_c, file = "export_apo_c.txt", row.names = FALSE, col.names = FALSE)


export_ape_pcs%>%
  summarise(min = min(mean_ape),
            max = max(mean_ape))

export_ape_c%>%
  summarise(min = min(mean_ape),
            max = max(mean_ape))

export_apo_pcs%>%
  summarise(min = min(mean_apo),
            max = max(mean_apo))

export_apo_c%>%
  summarise(min = min(mean_apo),
            max = max(mean_apo))

# now the same with the r squared
export_r_c <- table_power_5%>%
  filter(group == 'withoutPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_r = mean(r_squared, na.rm = T))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

export_r_pcs <- table_power_5%>%
  filter(group == 'withPCS')%>% 
  mutate(channel = as.numeric(channel)) %>%
  group_by(channel)%>%
  summarise(mean_r = mean(r_squared, na.rm = T))%>%
  arrange(channel)%>%
  mutate(channel = replace(channel, is.na(channel), "Gnd"))  

write.table(export_r_pcs, file = "export_r_pcs.txt", row.names = FALSE, col.names = FALSE)
write.table(export_r_c, file = "export_r_c.txt", row.names = FALSE, col.names = FALSE)

export_r_pcs%>%
  summarise(min = min(mean_r),
            max = max(mean_r))

export_r_c%>%
  summarise(min = min(mean_r),
            max = max(mean_r))

#--------- 7. check requirements-----------------------------
# I need data sets per group in order to check the normality requirement separately
shapiro_df_withPCS <- df_corr_frontal_filtered_group%>%
  filter(group == 'withPCS')

shapiro_df_withoutPCS <- df_corr_frontal_filtered_group%>%
  filter(group == 'withoutPCS')

# normality delta
df_corr_frontal_filtered_group%>%
  ggplot(aes(x = mean_delta_power))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# looks a bit weird but a similar kind of weird

shapiro.test(shapiro_df_withPCS$mean_delta_power)
shapiro.test(shapiro_df_withoutPCS$mean_delta_power)

# normality beta
shapiro_df_withPCS <- df_corr_central_filtered_group%>%
  filter(group == 'withPCS')
shapiro_df_withoutPCS <- df_corr_central_filtered_group%>%
  filter(group == 'withoutPCS')

df_corr_central_filtered%>%
  ggplot(aes(x = mean_beta_power))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')# looks really skew (bot equally skew in both groups, a bit worse in with PCS maybe) -> maybe use nonparametric stats

shapiro.test(shapiro_df_withPCS$mean_beta_power)# 8.953e-05
shapiro.test(shapiro_df_withoutPCS$mean_beta_power)# 0.001104

# beta 1
shapiro_df_withPCS <- df_corr_central1_filtered%>%
  filter(group == 'withPCS')
shapiro_df_withoutPCS <- df_corr_central1_filtered%>%
  filter(group == 'withoutPCS')

df_corr_central1_filtered%>%
  ggplot(aes(x = mean_beta1_power))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')# looks really skew (bot equally skew in both groups, a bit worse in with PCS maybe) -> maybe use nonparametric stats

shapiro.test(shapiro_df_withPCS$mean_beta1_power)# <.001
shapiro.test(shapiro_df_withoutPCS$mean_beta1_power)# <.001

# beta 2
shapiro_df_withPCS <- df_corr_central2_filtered%>%
  filter(group == 'withPCS')
shapiro_df_withoutPCS <- df_corr_central2_filtered%>%
  filter(group == 'withoutPCS')

df_corr_central2_filtered%>%
  ggplot(aes(x = mean_beta2_power))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')# looks really skew (bot equally skew in both groups, a bit worse in with PCS maybe) -> maybe use nonparametric stats

shapiro.test(shapiro_df_withPCS$mean_beta2_power)# <.001
shapiro.test(shapiro_df_withoutPCS$mean_beta2_power)# <.001

# normality aperiodic offset
shapiro_df_withPCS <- df_corr_apo%>%
  filter(group == 'withPCS')
shapiro_df_withoutPCS <- df_corr_apo%>%
  filter(group == 'withoutPCS')

df_corr_apo%>%
  ggplot(aes(x = mean_aperiodic_offset))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')+
  theme_classic()# looks quite normally distributed

shapiro.test(shapiro_df_withPCS$mean_aperiodic_offset)# 0.3911
shapiro.test(shapiro_df_withoutPCS$mean_aperiodic_offset)# 0.4375

# normality aperiodic exponent
shapiro_df_withPCS <- df_corr_ape%>%
  filter(group == 'withPCS')
shapiro_df_withoutPCS <- df_corr_ape%>%
  filter(group == 'withoutPCS')

df_corr_ape%>%
  ggplot(aes(x = mean_aperiodic_exponent))+
  geom_histogram(color = "black",
                 fill = "white", bins = sqrt(100))+
  facet_wrap(~group,scales = 'free')# looks different between the groups

shapiro.test(shapiro_df_withPCS$mean_aperiodic_exponent)# 0.01124
shapiro.test(shapiro_df_withoutPCS$mean_aperiodic_exponent)# 0.5865

# variance
leveneTest(mean_delta_power~group,data = df_corr_frontal_filtered_group)# not significant
leveneTest(mean_beta_power~group,data = df_corr_central_filtered_group)# not significant
leveneTest(mean_aperiodic_offset~group,data = df_corr_apo)# 0.0451
leveneTest(mean_aperiodic_exponent~group,data = df_corr_ape)# not significant
leveneTest(mean_beta1_power~group,data = df_corr_central1_filtered)# not significant
leveneTest(mean_beta2_power~group,data = df_corr_central2_filtered)# not significant


# conclusion: variances are not that big of a problem, normality is though! with the aperiodic offset we have normality in both groups
# but then there is no equal variances in that case
# => use NONPARAMETRIC Tests for beta, delta and the exponent + offset?

# ----- 8. boxplots and stats -------------------
# Define custom colors
color_palette <- c("without PCS" = '#02CAF5',
                   "with PCS" = "#F59541")

##---- 8.1 aperiodic exponent general ------------
# mean
df_corr_ape%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_aperiodic_exponent, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = function(p) sprintf("p = %.2g", p), test = 'wilcox.test', color = 'black')+
  labs(y = 'mean aperiodic exponent')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )


wilcox.test(mean_aperiodic_exponent~group, data = df_corr_ape, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)# 0.6057

df_corr_ape <- df_corr_ape%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_ape,mean_aperiodic_exponent~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')

# get the effsize
df_corr_ape%>%
  ungroup()%>% # apparently you have to ungroup here, otherwise, wilcox_effsize does not work
  wilcox_effsize(mean_aperiodic_exponent~group)

##-------8.2 aperiodic offset general-----------
# mean
df_corr_apo%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_aperiodic_offset, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = function(p) sprintf("p = %.2g", p), test = 'wilcox.test', color = 'black')+
  labs(y = 'mean aperiodic offset')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )

wilcox.test(mean_aperiodic_offset~group, data = df_corr_apo, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)
df_corr_apo <- df_corr_apo%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work
result <- coin::wilcox_test(data = df_corr_apo,mean_aperiodic_offset~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')
# get the effsize
df_corr_apo%>%
  ungroup()%>% # apparently you have to ungroup here, otherwise, wilcox_effsize does not work
  wilcox_effsize(mean_aperiodic_offset~group)

# but be aware that the without PCS group is a bit younger
# exponent
df_corr_ape%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = age, y = mean_aperiodic_exponent, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE)+
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.20, 0.15))+
  stat_cor(aes(color = "Correlation: "),method = "pearson", label.x = 60, label.y = 1.1,hjust=0)+
  labs(y = 'mean aperiodic exponent')

cor.test(df_corr_ape$age,df_corr_ape$mean_aperiodic_exponent)# r -.61, p 6.83e-06

# offset
df_corr_apo%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = age, y = mean_aperiodic_offset, color = group))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE)+
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.20, 0.15))+
  stat_cor(aes(color = "Correlation: "),method = "pearson", label.x = 60, label.y = 0.8,hjust=0)+
  labs(y = 'mean aperiodic offset')

cor.test(df_corr_apo$age,df_corr_apo$mean_aperiodic_offset)# significant (p = 0.036) r = -0.31


##------ 8.3 rel and absolute delta frontal---------------
df_corr_frontal_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power, color = group))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = TRUE, color = 'black')+
  labs(y = 'mean delta power [μV^2]')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 18)  # Adjust the size here
  )

g1 <- df_corr_frontal_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power, color = sex))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2)+
  labs(y = 'mean delta power [μV^2]')+
  theme_classic()+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )

# function for exporting figure: https://r-coder.com/save-plot-r/?utm_content=cmp-true

# in order to get the W statistics
wilcox.test(mean_delta_power~group, data = df_corr_frontal_filtered_group, 
            alternative = 'less',
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)

df_corr_frontal_filtered_group <-df_corr_frontal_filtered_group%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_frontal_filtered_group,mean_delta_power~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'less')

# get the effsize
df_corr_frontal_filtered_group%>%
  ungroup()%>%
  wilcox_effsize(mean_delta_power~group)

# absolute delta
df_corr_frontal_filtered_abs%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_delta_power_abs))+
  geom_boxplot()+
  geom_jitter(width = 0.3, height = 0, alpha = 0.1)

wilcox.test(mean_delta_power_abs~group, data = df_corr_frontal_filtered_abs, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)
t.test(mean_delta_power_abs~group, data = df_corr_frontal_filtered_abs, alternative = "less", paired = FALSE)

##----- 8.4 rel beta -----------------
df_corr_central_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta_power, color = group))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+                                         # Add p-value to plot
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = TRUE, color = 'black')+
  labs(y = 'mean beta power [μV^2]')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 18)  # Adjust the size here
  )

g2 <- df_corr_central_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta_power, color = sex))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2)+
  labs(y = 'mean beta power [μV^2]')+
  theme_classic()+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )



wilcox.test(mean_beta_power~group, data = df_corr_central_filtered_group, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)

df_corr_central_filtered_group <-df_corr_central_filtered_group%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_central_filtered_group,mean_beta_power~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')

# get effsize
df_corr_central_filtered_group%>%
  ungroup()%>%
  wilcox_effsize(mean_beta_power~group)

# exploratory: beta 1
df_corr_central1_filtered%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta1_power, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+                                         # Add p-value to plot
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = TRUE, color = 'black')+
  labs(y = 'mean beta power [μV^2] in central ROI')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )

wilcox.test(mean_beta1_power~group, data = df_corr_central1_filtered, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)

df_corr_central1_filtered <-df_corr_central1_filtered%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_central1_filtered,mean_beta1_power~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')

# get effsize
df_corr_central1_filtered%>%
  ungroup()%>%
  wilcox_effsize(mean_beta1_power~group)
# exploratory: beta 2
df_corr_central2_filtered%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = mean_beta2_power, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+                                         # Add p-value to plot
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = TRUE, color = 'black')+
  labs(y = 'mean beta power [μV^2] in central ROI')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 14)  # Adjust the size here
  )

wilcox.test(mean_beta2_power~group, data = df_corr_central2_filtered, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)

df_corr_central2_filtered <-df_corr_central2_filtered%>%
  mutate(group = as.factor(group))# otherwise wilcoxon_test from coin does not work

# in order to get the z value
result <- coin::wilcox_test(data = df_corr_central2_filtered,mean_beta2_power~group, comparisons = list(c('withPCS','withoutPCS')), alternative = 'two.sided')

# get effsize
df_corr_central2_filtered%>%
  ungroup()%>%
  wilcox_effsize(mean_beta2_power~group)

## --------- 8.5 tables of all EEG values ----
df_corr_frontal_filtered_group%>%
  group_by(group)%>%
  summarise(mean_delta = mean(mean_delta_power),
            sd_delta = sd(mean_delta_power))

df_corr_central_filtered_group%>%
  group_by(group)%>%
  summarise(mean_beta = mean(mean_beta_power),
            sd_beta = sd(mean_beta_power))

df_corr_ape%>%
  group_by(group)%>%
  summarise(mean_ape = mean(mean_aperiodic_exponent),
            sd_ape = sd(mean_aperiodic_exponent))

df_corr_apo%>%
  group_by(group)%>%
  summarise(mean_apo = mean(mean_aperiodic_offset),
            sd_apo = sd(mean_aperiodic_offset))

#------ 9. plot behavioral data and corr test ------
## ------- 9.1 just behavioral data ---------------
# TMT A
df_corr_frontal_filtered%>%
  ggplot(aes(x = group, y = tmt_a_time))+
  geom_boxplot()+
  geom_jitter(width = 0.3, height = 0, alpha = 0.1)

# FACIT
df_corr_frontal_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = facit_f_FS, color = group))+
  geom_boxplot(size = 0.75,outlier.colour = 'black', width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+                                         # Add p-value to plot
  geom_signif(comparisons = list(c("with PCS","without PCS")),map_signif_level = function(p) sprintf("p = %.2g", p),test = "t.test", color = 'black')+
  labs(y = 'FACIT Fatigue Scale [Range: 0-52]')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 15)  # Adjust the size here
  )

# TMT B-A
df_corr_frontal_filtered%>%
  ggplot(aes(x = group, y = tmt_b_minus_a))+
  geom_boxplot()+
  geom_jitter(width = 0.3, height = 0, alpha = 0.1)

# moca
df_corr_frontal_filtered%>%
  ggplot(aes(x = group, y = moca))+
  geom_boxplot()+
  geom_jitter(width = 0.3, height = 0, alpha = 0.1)

wilcox.test(moca~group, data = df_corr_frontal_filtered, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)

# FACIT and HADS-D
p8<- df_corr_frontal_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  ggplot(aes(x = hads_d_total_score,y = facit_f_FS, color = group))+
  geom_point(size = 2)+
  labs(y = 'FACIT Fatigue Scale [Range: 0-52]',
       x = 'HADS-D Score [Range: 0-21]')+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.20, 0.15))+
  stat_cor(aes(color = "Correlation: "),method = "spearman", label.x = 12, label.y = 40,hjust=0)
ggMarginal(p8, type = "densigram")

cor.test(df_corr_frontal_filtered_group$facit_f_FS,df_corr_frontal_filtered_group$hads_d_total_score, method = 'spearman', exact = FALSE)

# TMT with FACIT
df_corr_frontal_filtered%>%
  ggplot(aes(x = facit_f_FS, y = tmt_a_time))+
  geom_point()

cor.test(df_corr_frontal_filtered$tmt_a_time,df_corr_frontal_filtered$facit_f_FS, method = 'spearman', exact = FALSE)

df_corr_frontal_filtered%>%
  ggplot(aes(x = facit_f_FS, y = tmt_b_minus_a))+
  geom_point()

cor.test(df_corr_frontal_filtered$tmt_b_minus_a,df_corr_frontal_filtered$facit_f_FS, method = 'spearman', exact = FALSE)
#moca
cor.test(df_corr_frontal_filtered$moca,df_corr_frontal_filtered$facit_f_FS, method = 'spearman', exact = FALSE)

# HADS with TMT-A/B-A and MoCA
cor.test(df_corr_frontal_filtered_group$tmt_a_time,df_corr_frontal_filtered_group$hads_d_total_score, method = 'spearman', exact = FALSE)
cor.test(df_corr_frontal_filtered_group$tmt_b_minus_a,df_corr_frontal_filtered_group$hads_d_total_score, method = 'spearman', exact = FALSE)
cor.test(df_corr_frontal_filtered_group$moca,df_corr_frontal_filtered_group$hads_d_total_score, method = 'spearman', exact = FALSE)

# TMT-A with TMT B-A and MoCA
cor.test(df_corr_frontal_filtered_group$tmt_a_time,df_corr_frontal_filtered_group$tmt_b_minus_a, method = 'spearman', exact = FALSE)
cor.test(df_corr_frontal_filtered_group$tmt_a_time,df_corr_frontal_filtered_group$moca, method = 'spearman', exact = FALSE)

# TMT B-A and MoCA
cor.test(df_corr_frontal_filtered_group$tmt_b_minus_a,df_corr_frontal_filtered_group$moca, method = 'spearman', exact = FALSE)


# very high correlation
## --------- 9.2 corr tests ---------------------------
### ---- 9.2.1 relative delta power with TMT-A and TMT-B-A------------
df_corr_frontal_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = mean_delta_power,y = tmt_a_time, color = group))+
  geom_point()

cor.test(df_corr_frontal_filtered_group$mean_delta_power,df_corr_frontal_filtered_group$tmt_a_time, method = 'spearman', exact = FALSE)

df_corr_frontal_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = tmt_b_minus_a,y = mean_delta_power,color = group))+
  geom_point()

cor.test(df_corr_frontal_filtered_group$mean_delta_power,df_corr_frontal_filtered_group$tmt_b_minus_a, method = 'spearman', exact = FALSE)


### ------ 9.2.2 rel delta and moca ---------------
df_corr_frontal_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = moca,y = mean_delta_power,color = group))+
  geom_point()

cor.test(df_corr_frontal_filtered$mean_delta_power,df_corr_frontal_filtered$moca, method = 'spearman', exact = FALSE)

###-------- 9.2.3 relative and absolute delta and FACIT score ------------------------
p1<- df_corr_frontal_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  ggplot(aes(x = mean_delta_power,y = facit_f_FS, color = group))+
  geom_point(size = 2.5)+
  labs(y = 'FACIT Fatigue Scale [Range: 0-52]',
       x = 'mean delta power [μV^2]')+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.14, 0.15))+
  stat_cor(aes(color = "Correlation: "),method = "spearman", label.x = 2, label.y = 41,hjust=0)+
  theme(
    text = element_text(size = 15)  # Adjust the size here
  )

ggMarginal(p1, type = "densigram")


# one NA
cor.test(df_corr_frontal_filtered_group$mean_delta_power,df_corr_frontal_filtered_group$facit_f_FS, method = 'spearman', exact = FALSE)

# just curious = > divide into the two groups and test separately
p5<- df_corr_frontal_filtered_group%>%
  filter(group == 'withPCS')%>%
  ggplot(aes(x = mean_delta_power,y = facit_f_FS))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme_classic() +
  theme(legend.position = c(0.25, 0.15))

ggMarginal(p5, type = "densigram")


test <- df_corr_frontal_filtered_group%>%
  filter(group == 'withPCS')

cor.test(test$mean_delta_power,test$facit_f_FS) # 

p6 <- df_corr_frontal_filtered_group%>%
  filter(group == 'withoutPCS')%>%
  ggplot(aes(x = mean_delta_power,y = facit_f_FS))+
  geom_point()+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme_classic() +
  theme(legend.position = c(0.25, 0.15))

ggMarginal(p6, type = "densigram")

test2 <- df_corr_frontal_filtered_group%>%
  filter(group == 'withoutPCS')

cor.test(test2$mean_delta_power,test2$facit_f_FS) 


# absolute delta power
p1<- df_corr_frontal_filtered_abs%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  ggplot(aes(x = mean_delta_power_abs,y = facit_f_FS, color = group))+
  geom_point(size = 2.5)+
  labs(y = 'FACIT Fatigue Scale [Range: 0-52]',
       x = 'mean delta power [μV^2]')+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.14, 0.15))+
  stat_cor(aes(color = "Correlation: "),method = "spearman", label.x = 2, label.y = 41,hjust=0)+
  theme(
    text = element_text(size = 15)  # Adjust the size here
  )

ggMarginal(p1, type = "densigram")

###------------ 9.2.4 delta and hads d--------------------
p7<- df_corr_frontal_filtered_group%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = mean_delta_power,y = hads_d_total_score, color = group))+
  geom_point(size = 2.5)+
  labs(y = 'HADS-D Score [Range: 0-21]',
       x = 'mean delta power [μV^2]')+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.15, 0.85))+
  stat_cor(aes(color = "Correlation: "),method = "spearman", label.x = 2, label.y = 17,hjust=0)+
  theme(
    text = element_text(size = 15)  # Adjust the size here
  )

ggMarginal(p7, type = "densigram")
cor.test(df_corr_frontal_filtered_group$mean_delta_power,df_corr_frontal_filtered_group$hads_d_total_score, method = 'spearman', exact = FALSE)


### ------ 9.2.5 relative beta power and with TMT-A and TMT-B-A--------------------
df_corr_central_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta_power,y = tmt_a_time, color = group))+
  geom_point()

cor.test(df_corr_central_filtered_group$mean_beta_power,df_corr_central_filtered_group$tmt_a_time, method = 'spearman', exact = FALSE)


df_corr_central_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta_power,y = tmt_b_minus_a,color = group))+
  geom_point()

cor.test(df_corr_central_filtered_group$mean_beta_power,df_corr_central_filtered_group$tmt_b_minus_a, method = 'spearman', exact = FALSE)

# beta 1
df_corr_central1_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta1_power,y = tmt_a_time, color = group))+
  geom_point()

cor.test(df_corr_central1_filtered$mean_beta1_power,df_corr_central1_filtered$tmt_a_time, method = 'spearman', exact = FALSE)


df_corr_central1_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta1_power,y = tmt_b_minus_a,color = group))+
  geom_point()

cor.test(df_corr_central1_filtered$mean_beta1_power,df_corr_central1_filtered$tmt_b_minus_a, method = 'spearman', exact = FALSE)

# beta 2
df_corr_central2_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta2_power,y = tmt_a_time, color = group))+
  geom_point()

cor.test(df_corr_central2_filtered$mean_beta2_power,df_corr_central2_filtered$tmt_a_time, method = 'spearman', exact = FALSE)


df_corr_central2_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta2_power,y = tmt_b_minus_a,color = group))+
  geom_point()

cor.test(df_corr_central2_filtered$mean_beta2_power,df_corr_central2_filtered$tmt_b_minus_a, method = 'spearman', exact = FALSE)

###---- 9.2.6 relative beta power and with FACIT score --------------------
df_corr_central_filtered_group%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta_power,y = facit_f_FS,color = group))+
  geom_point()

cor.test(df_corr_central_filtered_group$mean_beta_power,df_corr_central_filtered_group$facit_f_FS, method = 'spearman', exact = FALSE)
#hads
cor.test(df_corr_central_filtered_group$mean_beta_power,df_corr_central_filtered_group$hads_d_total_score, method = 'spearman', exact = FALSE)
#moca
cor.test(df_corr_central_filtered_group$mean_beta_power,df_corr_central_filtered_group$moca, method = 'spearman', exact = FALSE)

# beta 1
df_corr_central1_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta1_power,y = facit_f_FS,color = group))+
  geom_point()

cor.test(df_corr_central1_filtered$mean_beta1_power,df_corr_central1_filtered$facit_f_FS, method = 'spearman', exact = FALSE)
cor.test(df_corr_central1_filtered$mean_beta1_power,df_corr_central1_filtered$hads_d_total_score, method = 'spearman', exact = FALSE)
cor.test(df_corr_central1_filtered$mean_beta1_power,df_corr_central1_filtered$moca, method = 'spearman', exact = FALSE)

# beta 2
df_corr_central2_filtered%>%
  group_by(group)%>%
  ggplot(aes(x = mean_beta2_power,y = facit_f_FS,color = group))+
  geom_point()

cor.test(df_corr_central2_filtered$mean_beta2_power,df_corr_central2_filtered$facit_f_FS, method = 'spearman', exact = FALSE)
cor.test(df_corr_central2_filtered$mean_beta2_power,df_corr_central2_filtered$hads_d_total_score, method = 'spearman', exact = FALSE)
cor.test(df_corr_central2_filtered$mean_beta2_power,df_corr_central2_filtered$moca, method = 'spearman', exact = FALSE)

# beta and delta
corr_power <- cbind(df_corr_central_filtered_group,df_corr_frontal_filtered_group)
cor.test(corr_power$mean_beta_power,corr_power$mean_delta_power, method = 'spearman', exact = FALSE)

### --------9.2.7 aperiodic exponent with everything --------
corr_ape_apo <- cbind(df_corr_central_filtered_group,df_corr_frontal_filtered_group,df_corr_ape,df_corr_apo)
cor.test(corr_ape_apo$mean_aperiodic_exponent...21,corr_ape_apo$mean_delta_power, method = 'spearman', exact = FALSE)
cor.test(corr_ape_apo$mean_aperiodic_exponent...21,corr_ape_apo$mean_beta_power, method = 'spearman', exact = FALSE)
cor.test(df_corr_ape$mean_aperiodic_exponent,df_corr_ape$facit_f_FS, method = 'spearman', exact = FALSE)
cor.test(df_corr_ape$mean_aperiodic_exponent,df_corr_ape$hads_d_total_score, method = 'spearman', exact = FALSE)
cor.test(df_corr_ape$mean_aperiodic_exponent,df_corr_ape$tmt_a_time, method = 'spearman', exact = FALSE)
cor.test(df_corr_ape$mean_aperiodic_exponent,df_corr_ape$tmt_b_minus_a, method = 'spearman', exact = FALSE)
cor.test(df_corr_ape$mean_aperiodic_exponent,df_corr_ape$moca, method = 'spearman', exact = FALSE)
cor.test(corr_ape_apo$mean_aperiodic_exponent...21,corr_ape_apo$mean_aperiodic_offset, method = 'spearman', exact = FALSE)


corr_ape_apo%>%
  group_by(group...2)%>%
  ggplot(aes(x = mean_aperiodic_exponent...19,y = mean_delta_power,color = group...2))+
  geom_point()

p1<- corr_ape_apo%>%
  mutate(group = fct_recode(group...2,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"))%>%
  ggplot(aes(x = mean_delta_power,y = mean_aperiodic_exponent...19, color = group))+
  geom_point(size = 2.5)+
  labs(y = 'mean aperiodic exponent',
       x = 'mean delta power [μV^2]')+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(legend.position = c(0.75, 0.77))+
  stat_cor(aes(color = "Correlation: "),method = "spearman", label.x = 0.1, label.y = 0.3,hjust=0)+
  theme(
    text = element_text(size = 15)  # Adjust the size here
  )

ggMarginal(p1, type = "densigram")
### ---------9.2.8 aperiodic offset with everything -------------
cor.test(corr_ape_apo$mean_aperiodic_offset,corr_ape_apo$mean_delta_power, method = 'spearman', exact = FALSE)
cor.test(corr_ape_apo$mean_aperiodic_offset,corr_ape_apo$mean_beta_power, method = 'spearman', exact = FALSE)
cor.test(df_corr_apo$mean_aperiodic_offset,df_corr_ape$facit_f_FS, method = 'spearman', exact = FALSE)
cor.test(df_corr_apo$mean_aperiodic_offset,df_corr_ape$hads_d_total_score, method = 'spearman', exact = FALSE)
cor.test(df_corr_apo$mean_aperiodic_offset,df_corr_ape$tmt_a_time, method = 'spearman', exact = FALSE)
cor.test(df_corr_apo$mean_aperiodic_offset,df_corr_ape$tmt_b_minus_a, method = 'spearman', exact = FALSE)
cor.test(df_corr_apo$mean_aperiodic_offset,df_corr_ape$moca, method = 'spearman', exact = FALSE)

# ---------- 10. r squared ----------------
table_power_5%>%
  ggplot(aes(x = r_squared, y = aperiodic_exponent))+
  geom_point()
cor.test(table_power_5$aperiodic_exponent, table_power_5$r_squared)

table_power_5%>%
  ggplot(aes(x = r_squared, y = aperiodic_offset))+
  geom_point()
cor.test(table_power_5$aperiodic_offset, table_power_5$r_squared)

# ---------- 11. permutation tests -------------------------
# comparing apo and ape at every channel
table_apo_filtered <- table_apo_filtered%>%
  mutate(group = as.factor(group))

table_ape_filtered <- table_ape_filtered%>%
  mutate(group = as.factor(group))

# Create permutation test function
permutation_function <- function(nsim,df,chan,offset_or_exponent){
  res <- numeric(nsim) ## set aside space for results
  df <- df%>%
    filter(channel == chan)
  for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df))
    bdat <- transform(df,group = group[perm])
    ## compute & store difference in means; store the value
    tt <- t.test(bdat[[offset_or_exponent]]~group,data=bdat,var.equal=FALSE)
    res[i] <- tt$statistic
  }
  obs <- t.test(df[[offset_or_exponent]]~group,data=df,var.equal=FALSE)
  ## append the observed value to the list of results
  res <- c(res,obs$statistic)
  res <<- res
  obs <<- obs$statistic
}

# Create permutation test function for wilcox test
permutation_function_wilcox <- function(nsim,df,chan,offset_or_exponent){
  res <- numeric(nsim) ## set aside space for results
  df <- df%>%
    filter(channel == chan)
  for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df))
    bdat <- transform(df,group = group[perm])
    ## compute & store difference in means; store the value
    w <- wilcox.test(bdat[[offset_or_exponent]]~group,data=bdat,var.equal=FALSE)
    res[i] <- w$statistic
  }
  obs <- wilcox.test(df[[offset_or_exponent]]~group,data=df,var.equal=FALSE)
  ## append the observed value to the list of results
  res <- c(res,obs$statistic)
  res <<- res
  obs <<- obs$statistic
}

# create a function for visualisation
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

visualize_results_wilcox <- function(res,obs){
  result_sim <- sum(res > obs)/(length(res)-1)
  Histogramm_pre <- hist(res,
                         ylab = 'Anzahl',
                         xlab = substitute(paste(italic('W'), ' values')),
                         main = substitute(paste(italic('W'), ' values with observed value')))
  Histogramm <- abline(v=obs,col="red")
  return(result_sim)
  return(Histogramm)
}


# now do this for every channel: offset 
set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,1,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.298

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,2,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.241

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,3,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.467

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,4,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.554

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,5,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.342

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,6,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.385

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,7,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.255

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,8,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.396

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,9,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.296

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,10,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.654

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,11,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.342

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,12,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.319

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,13,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.431

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,14,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.401

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,15,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.608

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,16,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.611

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,17,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.394

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,18,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.643

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,19,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.614

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,21,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.79

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,22,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.101

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,23,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.537

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,24,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.765

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,25,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.333

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,26,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.617

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,27,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.373

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,28,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.638

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,29,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.818

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,30,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.579

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,33,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.385

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,34,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.685

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,35,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.337

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,36,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.417

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,37,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.357

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,38,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.505

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,39,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.583

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,40,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.606

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,41,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.219

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,42,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.188

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,43,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.474

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,44,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.389

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,45,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.273

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,46,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.452

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,47,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.211

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,48,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.309

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,49,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.351

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,50,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.784

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,51,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.235

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,52,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.479

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,53,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.54

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,54,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.545

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,55,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.392

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,56,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.505

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,57,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.662

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,58,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.755

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,59,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.472

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,60,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.821

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,61,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.279

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,62,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.436

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,63,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.52

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,64,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.587

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,65,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.456

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,66,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.455

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,67,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.268

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,68,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.337

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,69,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.393

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,70,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.422

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,71,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.415

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,72,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.281

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,73,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.397

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,74,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.335

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,75,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.561

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,76,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.359

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,77,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.447

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,78,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.258

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,79,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.28

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,80,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.28

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,81,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.29

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,82,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.397

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,83,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.66

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,84,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.329

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,85,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.243

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,86,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.238

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,87,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.393

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,88,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.593

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,89,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.519

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,90,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.265

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,91,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.542

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,92,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.349

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,93,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.721

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,94,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.395

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,95,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.559

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,96,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.655

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,97,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.124

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,98,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.179

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,99,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.444

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,100,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.406

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,101,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.533

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,102,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.361

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,103,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.428

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,104,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.433

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,105,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.412

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,106,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.44

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,107,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.748

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,108,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.688

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,109,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.499

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,110,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.481

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,111,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.636

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,112,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.616

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,113,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.317

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,114,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.32

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,115,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.319

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,116,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.284

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,117,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.723

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,118,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.416

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,119,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.686

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,120,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.755

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,121,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.499

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,122,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.675

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,123,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.251

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,124,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.914

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,125,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.233

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,126,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.329

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,127,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.507

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,128,'aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.418

set.seed(2)
permutation_function_wilcox(1000,table_apo_filtered,'Gnd','aperiodic_offset')
visualize_results_wilcox(res,obs)# 0.672


# aperiodic exponent
set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,1,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.062

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,2,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.062

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,3,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.137

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,4,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.507

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,5,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.076

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,6,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.071

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,7,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.057

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,8,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.13

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,9,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.141

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,10,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.13

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,11,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.319

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,12,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.366

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,13,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.313

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,14,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.439

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,15,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.534

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,16,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.357

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,17,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.435

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,18,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.578

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,19,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.308

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,21,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.847

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,22,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.239

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,23,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.399

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,24,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.762

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,25,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.349

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,26,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.634

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,27,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.198

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,28,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.306

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,29,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.46

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,30,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.257

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,33,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.159

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,34,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.149

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,35,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.063

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,36,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.122

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,37,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.106

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,38,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.129

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,39,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.126

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,40,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.436

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,41,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.496

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,42,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.24

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,43,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.082

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,44,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.156

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,45,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.16

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,46,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.481

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,47,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.22

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,48,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.214

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,49,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.113

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,50,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.538

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,51,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.253

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,52,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.34

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,53,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.413

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,54,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.429

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,55,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.265

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,56,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.32

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,57,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.577

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,58,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.501

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,59,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.407

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,60,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.249

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,61,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.337

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,62,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.246

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,63,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.223

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,64,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.26

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,65,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.067

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,66,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.121

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,67,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.118

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,68,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.046

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,69,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.072

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,70,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.102

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,71,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.062

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,72,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.08

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,73,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.1

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,74,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.11

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,75,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.293

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,76,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.225

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,77,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.198

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,78,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.082

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,79,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.053

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,80,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.078

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,81,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.164

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,82,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.162

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,83,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.312

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,84,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.249

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,85,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.487

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,86,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.046

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,87,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.133

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,88,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.252

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,89,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.372

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,90,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.442

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,91,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.102

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,92,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.119

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,93,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.569

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,94,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.49

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,95,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.708

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,96,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.549

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,97,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.114

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,98,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.332

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,99,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.278

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,100,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.277

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,101,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.371

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,102,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.447

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,103,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.561

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,104,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.235

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,105,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.389

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,106,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.469

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,107,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.582

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,108,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.413

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,109,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.212

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,110,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.365

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,111,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.602

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,112,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.459

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,113,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.447

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,114,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.195

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,115,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.206

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,116,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.144

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,117,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.419

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,118,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.255

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,119,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.263

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,120,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.594

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,121,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.131

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,122,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.336

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,123,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.475

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,124,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.249

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,125,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.267

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,126,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.207

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,127,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.379

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,128,'aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.212

set.seed(2)
permutation_function_wilcox(1000,table_ape_filtered,'Gnd','aperiodic_exponent')
visualize_results_wilcox(res,obs)# 0.702


# visualize significant permutations
table_ape_filtered%>%
  filter(channel == 86)%>%
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"),
         group = fct_relevel(group, "with PCS", "without PCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = aperiodic_exponent, color = group))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  labs(y = 'mean AE at channel 86')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 18)  # Adjust the size here
  )

table_ape_filtered%>%
  filter(channel == 68)%>%#68
  mutate(group = fct_recode(group,
                            "with PCS" = "withPCS",
                            "without PCS" = "withoutPCS"),
         group = fct_relevel(group, "with PCS", "without PCS"))%>%
  group_by(group)%>%
  ggplot(aes(x = group, y = aperiodic_exponent, color = group))+
  geom_boxplot(size = 0.75,outlier.shape = NA, width=0.5)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2)+
  labs(y = 'mean AE at channel 68')+
  scale_color_manual(values = color_palette) +
  theme_classic()+
  guides(color = FALSE)+
  theme(
    text = element_text(size = 18)  # Adjust the size here
  )


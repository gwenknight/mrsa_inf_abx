#### Data analysis - antibiotic use and MRSA incidence
## MRSA incidence from Cassini et al 2018 (see data_manip.R for aggregation)
## Antibiotic use data from https://www.ecdc.europa.eu/en/antimicrobial-consumption/surveillance-and-disease-data/database

## Needed libraries
library(tidyverse)
library(naniar)
library(ggpmisc)
library(broom)
library(nlme)
library(here)
setwd(here())
theme_set(theme_bw(base_size = 14))

##******************************** Data ***********************************#####################
## Read in data - combined infection incidence and antibiotic use 
data <- read_csv("data/data_amr_abx_use_macotra.csv")
data <- data %>% replace_with_na_all(condition = ~.x == "-") # remove all "-" and replace with NA
data[,2:128] <- apply(data[,2:128], 2, as.numeric)

# Extract setting
data_abx <- data %>% pivot_longer(cols = "J01 C antibacterials_for_systemic":"J01XE CH nitrofuran") %>% 
  filter(!is.na(value)) %>% 
  mutate(ATC_code = gsub( " .*$", "", name ),
         setting_drug = sub( ".*? ", "", name),
         setting = gsub( " .*$", "", setting_drug),
         drug = gsub(".*? ", "", name)) %>% select(-c("setting_drug"))

data_abx[which(data_abx$setting == "C"), "setting"] <- "Community"
data_abx[which(data_abx$setting == "H"), "setting"] <- "Hospital"
data_abx[which(data_abx$setting == "CH"), "setting"] <- "Community & Hospital"

# Label 3 key countries for MACOTRA
w<-which(data_abx$country %in% c("UK","Netherlands","France") ) 
data_abx$macotra <- ""
data_abx[w,"macotra"] <- data_abx[w,"country"]

# Make numeric
data_abx$value <- as.numeric(data_abx$value)
data_abx$infection_per_100000 <- as.numeric(data_abx$infection_per_100000)

# Remove NA values 
data_abx <- data_abx %>% filter(!is.na(value),!is.na(infection_per_100000)) 

# Make additional columns
data_abx$log_data <- log(data_abx$infection_per_100000) # log of infection data
data_abx$sett_drug <- interaction(data_abx$setting,data_abx$drug,data_abx$ATC_code) # combined flag of setting and drug

##******************************** Relationships between drug use and infection incidence ***********************************#####################
## Plot everything by drug 
theme_set(theme_bw(base_size = 10))
ggplot(data_abx, aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_wrap(~drug, scales = "free") + 
  geom_smooth(method='lm', aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France")), pch = 3) + 
  scale_x_continuous("Antibiotic usage (DDD per 1000 inhabitants and per day)") + 
  scale_y_continuous("Infection incidence per 100,000 inhabitants")
ggsave("plots/data_by_drug.pdf", width = 25, height = 20)

## Plot everything by drug and setting with linear regression line
formula <- y ~ x   
ggplot(data_abx, aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(drug~setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=0,vjust=0.1) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula, 
               label.x.npc = "right", label.y.npc = 0.15) 
ggsave("data_by_drug&setting.pdf", width = 12, height = 20)

## Explore specific drug classes
ggplot(data_abx %>% filter(drug %in% c("penicillins_with_extended_spectrum","b-lactam",
                                       "blactamase_sensitive_penicillins","blactamase_resistant_pencillins",
                                       "blactamase_inhibitors","combination_penicillins_incl_blactamase_inhibitors",
                                       "other_blactams")), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(Country %in% c("UK","Netherlands","France"))%>% filter(drug %in% c("pen","b-lact")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)")
ggsave("fig_blact.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(drug %in% c("Quinolones", "1G_quinolones", "2G_quinolones", "3G_quinolones", "fluoroquinolones", "other_quinolones")), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(drug~setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(Country %in% c("UK","Netherlands","France"))%>% filter(drug %in% c("streptogramins")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)")
ggsave("fig_quinolones.pdf")

ggplot(data_abx %>% filter(drug %in% c("macrolides_lincosamides_streptogramins", 
                                       "short_acting_macrolides", "intermediate_acting_macrolides", 
                                       "long_acting_macrolides", "macrolides", "lincosamides", "streptogramins")), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(Country %in% c("UK","Netherlands","France"))%>% filter(drug %in% c("streptogramins")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)")
ggsave("fig_macrolide.pdf", width = 20, height = 10)

##******************************** Model fitting analysis ***********************************#####################
#### Linear model 
## Is abx use linked to infection across all data?
all_data_model <- lm(data_abx$infection_per_100000 ~ data_abx$value)

ggplot(data_abx, aes(x=value, y = infection_per_100000, label = macotra)) + geom_point() + 
  geom_smooth(method = "lm", formula = formula) + 
  ggtitle("All data") + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)")
ggsave("all_data.pdf")

summary(all_data_model)

#### Fit linear model to each drug across all settings
u <- unique(data_abx$sett_drug) # which drugs and settings? 
table_store<- c() # drug / intercept / value / p value / Rsquared

# linear
fit <- lmList(infection_per_100000 ~ value | sett_drug, data=data_abx)
r2 <- sapply(fit,function(x) summary(x)$r.squared)
s <- sapply(fit,function(x) tidy(x))
A <- sapply(fit,function(x) AIC(x))

for(i in u){
  w_col = which(colnames(s) == i)
  table_store <- rbind(table_store, 
                       c(unlist(strsplit(i, ".", fixed = TRUE)), round(s[,w_col][[2]],4),signif(s[,w_col][[5]][2],3),round(r2[which(names(r2)==i)],4),signif(A[which(names(A)==i)],6),"linear"))
}

# exponential
fit <- lmList(log(infection_per_100000) ~ value | sett_drug, data=data_abx)
#fit <- lmList(infection_per_100000 ~ exp(value) | sett_drug, data=data_abx.na)
r2 <- sapply(fit,function(x) summary(x)$r.squared)
s <- sapply(fit,function(x) tidy(x))
A <- sapply(fit,function(x) AIC(x))

for(i in u){
  w_col = which(colnames(s) == i)
  table_store <- rbind(table_store, 
                       c(unlist(strsplit(i, ".", fixed = TRUE)), round(s[,w_col][[2]],4),signif(s[,w_col][[5]][2],3),round(r2[which(names(r2)==i)],4),signif(A[which(names(A)==i)],6),"exponential"))
}
colnames(table_store) <- c("setting","drug","ATC code","intercept","value","pval_value","Rsquared","AIC","model")
table_store <- as.data.frame(table_store)

# Analyse stored data - should be min of AIC or Rsquared? 
table_stored <- table_store %>% group_by(setting,drug) %>% 
  mutate(max_R = max(Rsquared), min_AIC = min(AIC), min_pval_value = min(pval_value)) %>% 
  ungroup() %>% 
  mutate(significant = ifelse(pval_value < 0.05, ifelse(pval_value  < 0.01,"**","*"),""))

best_aic <- table_stored %>% group_by(setting,drug) %>% filter(AIC == min_AIC) 
best_rsq <- table_stored %>% group_by(setting,drug) %>% filter(Rsquared == max_R) 

# Check each has a "best fit" model
nrow(best_rsq) 
nrow(best_aic) 
length(u)  

# Different dependening on prediction (AIC) / current fit (Rsquared) distinction
best_aic %>% filter(pval_value < 0.01) %>% print(n=Inf)
best_rsq %>% filter(pval_value < 0.01) %>% print(n=Inf)

### Export table
write.table(table_stored[order(table_stored$drug),c("setting","drug","model","intercept","value","Rsquared","AIC","pval_value","significant")], 
            file = "table_supp.txt", sep = ",", quote = FALSE, row.names = F)

### plot those significant 
best_rsq$pval_value <- as.numeric(best_rsq$pval_value)
df_linear <- best_rsq %>% filter(pval_value < 0.01, model == "linear")
df_exponential <- best_rsq %>% filter(pval_value < 0.01, model == "exponential")

new_data_linear <- c()
for(i in 1:dim(df_linear)[1]){
  new_data_linear <- rbind(new_data_linear, data_abx %>% filter(drug == df_linear$drug[i] & setting == df_linear$setting[i]))
}

ggplot(new_data_linear, aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = new_data_linear %>% filter(Country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  guides(fill="none", col = "none")
ggsave("linear_just_significant_0.01.pdf", width = 22, height = 8)

new_data_exponential <- c()
for(i in 1:dim(df_exponential)[1]){
  new_data_exponential <- rbind(new_data_exponential, data_abx %>% filter(drug == df_exponential$drug[i] & setting == df_exponential$setting[i]))
}

formula = y ~ x
ggplot(new_data_exponential, aes(x=value, y = log_data, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(setting ~ drug, scales = "free",labeller = labeller(facet_category = label_wrap_gen(width = 16))) + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = new_data_exponential %>% filter(Country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(new_data_exponential$log_data + 10)), "Log of MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  guides(fill="none", col = "none")
ggsave("exponential_just_significant_0.01.pdf", width = 22, height = 8)



##### Figure 4
formula = y ~ x
ggplot(data_abx %>% filter(drug %in% c("b-lactam","blactamase_sensitive_penicillins")), 
       aes(x=value, y = log_data, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(Country %in% c("UK","Netherlands","France"))%>% filter(drug %in% c("pen","b-lact")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$log_data + 1)), "Log (base 10) MRSA incidence per 100,000") + 
  scale_x_continuous("Total usage (DDD 1000 inhabitants and per day)")
ggsave("fig_blact.pdf", width = 20, height = 10)

ggplot(data_abx %>% filter(drug %in% c("b-lactam","other_blactams"), setting %in% c("Community","Hospital")), 
       aes(x=value, y = log_data, label = macotra, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(drug ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(Country %in% c("UK","Netherlands","France"))%>% filter(drug %in% c("pen","b-lact")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0, size =5) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous("Log (base 10) MRSA incidence per 100,000") + 
  scale_x_continuous("Total usage (DDD 1000 inhabitants and per day)")
ggsave("fig_blact_rotate.pdf", width = 15, height = 10)

# Ery ciprofloxacin trimethoprim important in THD

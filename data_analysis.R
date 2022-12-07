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
library(patchwork)
library(ggh4x) # for facet_nested
setwd(here())
theme_set(theme_bw(base_size = 14))

##******************************** Data ***********************************#####################
## Read in data - combined infection incidence and antibiotic use 
data <- read_csv("data/data_amr_abx_use_macotra.csv")
data <- data %>% replace_with_na_all(condition = ~.x == "-") # remove all "-" and replace with NA
data[,2:131] <- apply(data[,2:131], 2, as.numeric)

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
data_abx$macotra <- factor(data_abx$macotra, levels = c("France","Netherlands","UK", ""))

# Make numeric
data_abx$value <- as.numeric(data_abx$value)
data_abx$infection_per_100000 <- as.numeric(data_abx$infection_per_100000)

# Remove NA values 
data_abx <- data_abx %>% filter(!is.na(value),!is.na(infection_per_100000)) 

# Make additional columns
data_abx$log_data <- log(data_abx$infection_per_100000) # log of infection data
data_abx$sett_drug <- interaction(data_abx$setting,data_abx$drug,data_abx$ATC_code) # combined flag of setting and drug

## Additional codes for specific drug groupings 
data_abx[which(data_abx$drug == "1G_quinolones"),"ATC_code"] = "J01M_frst"
data_abx[which(data_abx$drug == "2G_quinolones"),"ATC_code"] = "J01M_scnd"
data_abx[which(data_abx$drug == "3G_quinolones"),"ATC_code"] = "J01M_thrd"

data_abx[which(data_abx$drug == "short_acting_macrolides"),"ATC_code"] = "J01FA_shrt"
data_abx[which(data_abx$drug == "intermediate_acting_macrolides"),"ATC_code"] = "J01FA_intr"
data_abx[which(data_abx$drug == "long_acting_macrolides"),"ATC_code"] = "J01FA_long"

# Extract ATC code family and assign summary classes 
data_abx <- data_abx %>% rowwise() %>% mutate(ATC_code_family = substring(ATC_code,1,4)) %>% 
  mutate(summary_class = ifelse(drug %in% c("all_b_lactam","b-lactam","other_blactams",
                                            "sulfonamides_and_trimethoprim","macrolides_lincosamides_streptogramins",
                                            "Quinolones","macrolides","fluoroquinolones"),
                                "Yes","No"))



data_abx$drug <- factor(data_abx$drug, levels = c("antibacterials_for_systemic", 
                                                  "tetracyclines", 
                                                  "amphenicols", 
                                                  "all_b_lactam", 
                                                  "b-lactam",
                                                  "penicillins_with_extended_spectrum", 
                                                  "blactamase_resistant_penicillins", "blactamase_inhibitors", 
                                                  "combination_penicillins_incl_blactamase_inhibitors", 
                                                  "other_blactams", 
                                                  "1G_cephalosporins", "2G_cephalosporins", "3G_cephalosporins", 
                                                  "blactamase_sensitive_penicillins", 
                                                  "4G_cephalosporin", "monobactams", "carbapenems", "other_cephalosporins_and_penems", 
                                                  "Quinolones", 
                                                  "fluoroquinolones", "other_quinolones", 
                                                  "1G_quinolones", "2G_quinolones", "3G_quinolones", 
                                                  "sulfonamides_and_trimethoprim", 
                                                  "trimethoprim_and_derivatives", "short_acting_sulfonamides", "intermediate_acting_sulfonamides", 
                                                  "sulfonamides_and_trimethoprim_combs",
                                                  "macrolides_lincosamides_streptogramins", 
                                                  "macrolides", "lincosamides", "streptogramins", 
                                                  "short_acting_macrolides", "intermediate_acting_macrolides", "long_acting_macrolides", 
                                                  "aminoglycosides", 
                                                  "other_aminoglycosides", 
                                                  "glycopeptides", "polymyxins", "steroid_antibacterials", "imadazole", "nitrofuran"))

data_abx$ATC_code <- factor(data_abx$ATC_code, c("J01", "J01A", "J01B", "J01C_D", "J01C", "J01CA", "J01CE", 
                                                 "J01CF", "J01CG", "J01CR", "J01D", "J01DB", "J01DC", "J01DD", 
                                                 "J01DE", "J01DF", "J01DH", "J01DI", "J01E", "J01EA", "J01EB", 
                                                 "J01EC", "J01EE", "J01F", 
                                                 "J01FA", "J01FA_shrt", "J01FA_intr", "J01FA_long", 
                                                 "J01FF", "J01FG", "J01G", "J01GB", "J01M", "J01MA", "J01M_frst", 
                                                 "J01M_scnd", "J01M_thrd", "J01MB", "J01XA", "J01XB", 
                                                 "J01XC", "J01XD", "J01XE"))
##******************************** Relationships between drug use and infection incidence ***********************************#####################
## Plot everything by drug 
theme_set(theme_bw(base_size = 8))
ggplot(data_abx, aes(x=value, y = infection_per_100000, group = setting)) + geom_point(aes(col = setting)) + 
  facet_wrap(~drug, scales = "free") + 
  geom_smooth(method='lm', aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France")), pch = 3) + 
  scale_x_continuous("Antibiotic usage (DDD per 1000 inhabitants and per day)") + 
  scale_y_continuous("Infection incidence per 100,000 inhabitants")
ggsave("plots/data_by_drug.pdf", width = 20, height = 7)


## Plot everything by drug and setting with linear regression line
formula <- y ~ x   
ggplot(data_abx, aes(x=value, y = infection_per_100000, group = setting)) + geom_point(aes(col = setting)) + 
  facet_grid(drug~setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = setting)) + 
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France")), pch = 3) + 
  scale_x_continuous("Antibiotic usage (DDD per 1000 inhabitants and per day)") + 
  scale_y_continuous("Infection incidence per 100,000 inhabitants") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  theme(legend.position="bottom")
ggsave("plots/data_by_drug&setting.pdf", width = 10, height = 14)


## Just beta-lactams
ggplot(data_abx %>% filter(ATC_code_family %in% c("J01C", "J01D")), 
       aes(x=value, y = infection_per_100000,  group = setting)) + 
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) +
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01C", "J01D")), pch = 3) + 
  geom_text(aes(label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Beta-lactams, all data, linear models") + 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_blact.pdf", width = 25, height = 10)


# blactamase tenuous - driven by only two data points
data_abx %>% filter(drug %in% c("blactamase_inhibitors")) %>% filter(value > 0)
# Carbapenems? 
data_abx %>% filter(drug %in% c("carbapenems")) %>% filter(value > 0) %>% arrange(desc(value)) %>% print(n=Inf)


ggplot(data_abx %>% filter(ATC_code_family %in% c("J01C", "J01D"),!country %in% c("Portugal")), 
       aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), 
                    values = c("white","grey90")) +
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula) +
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, 
                                            "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01C", "J01D"),!country %in% c("Portugal")), pch = 3) + 
  geom_text(aes(label = macotra),hjust=-0.17,vjust=0) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Beta-lactams, without Portugal data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_blact_remove_top1.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(ATC_code_family %in% c("J01C", "J01D"),
                           !country %in% c("Italy","Portugal")), 
       aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = factor(summary_class))) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01C", "J01D")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") 
ggsave("plots/fig_blact_remove_top2.pdf", width = 25, height = 10)

data_abx %>% filter(drug == "blactamase_inhibitors") %>% arrange(desc(value)) # Portugal outlier in blactamase_inhibitors use
data_abx %>% filter(drug == "1G_cephalosporins") %>% arrange(desc(value)) # Finland outlier in first gen ceph use
data_abx %>% filter(drug == "other_cephalosporins_and_penems") %>% arrange(desc(value)) # Slovenia outlier in other_cephalosporins_and_penems use

##### Exponential 
ggplot(data_abx %>% filter(ATC_code_family %in% c("J01C", "J01D")), 
       aes(x=value, y = log(infection_per_100000), group = setting)) + 
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01C", "J01D")), pch = 3) + 
  geom_text(aes(label = macotra), hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(log(data_abx$infection_per_100000 + 10))), "Log. MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Beta-lactams, all data, exponential model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_blact_exponential.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(ATC_code_family %in% c("J01C", "J01D"), 
                           !country == "Portugal"), 
       aes(x=value, y = log(infection_per_100000), group = setting)) + 
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01C", "J01D")), pch = 3) + 
  geom_text(aes(label=macotra), hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(log(data_abx$infection_per_100000 + 10))), "Log. MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Beta-lactams, without Portugal data, exponential model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_blact_exponential_noPortugal.pdf", width = 25, height = 10)


####### Other drugs
### Quinolones
ggplot(data_abx %>% filter(ATC_code_family == "J01M"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code == "J01M"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Quinolones, all data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_quin_linear.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(ATC_code_family == "J01M"), 
       aes(x=value, y = log(infection_per_100000), group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code == "J01M"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(log(data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Quinolones, all data, exponential model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_quin_exponential.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(!country == "Portugal", ATC_code_family == "J01M"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(!country == "Portugal", ATC_code == "J01M"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Quinolones, without Portugal data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_quin_linear_noPortugal.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(!country == "Portugal",ATC_code_family == "J01M"), 
       aes(x=value, y = log(infection_per_100000), group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(!country == "Portugal",ATC_code == "J01M"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(log(data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Quinolones, without Portugal data, exponential model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_quin_exponential_noPortugal.pdf", width = 25, height = 10)

### Macrolides
ggplot(data_abx %>% filter(ATC_code_family == "J01F"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code == "J01F"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Macrolides, Lincosamides and Streptogramins, all data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_mac_linear.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(ATC_code_family == "J01F"), 
       aes(x=value, y = log(infection_per_100000), group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code == "J01F"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(log(data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Macrolides, Lincosamides and Streptogramins, all data, exponential model") + 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_mac_exponential.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(!country == "Portugal", ATC_code_family == "J01F"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(!country == "Portugal", ATC_code == "J01F"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Macrolides, Lincosamides and Streptogramins, without Portugal data, linear model") + 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_mac_linear_noPortugal.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(!country == "Portugal",ATC_code_family == "J01F"), 
       aes(x=value, y = log(infection_per_100000), group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(!country == "Portugal",ATC_code == "J01F"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(log(data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Macrolides, Lincosamides and Streptogramins, without Portugal data, exponential model") + 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),         
        panel.grid.major=element_line(colour="black"),         
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_mac_exponential_noPortugal.pdf", width = 15, height = 10)

### Tetracyclines
ggplot(data_abx %>% filter(ATC_code_family == "J01A"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code == "J01A"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Tetracyclines, all data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_tet_linear.pdf", width = 25, height = 10)

### Trimethoprim
ggplot(data_abx %>% filter(ATC_code_family == "J01E"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code == "J01E"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Sulfonamides and trimethoprim, all data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_trim_linear.pdf", width = 25, height = 10)

ggplot(data_abx %>% filter(!country == "Portugal", ATC_code_family == "J01E"), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.1) + 
  scale_fill_manual("Summary\nclass",breaks = c("No","Yes"), values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(!country == "Portugal", ATC_code == "J01E"), pch = 3) + 
  geom_text(aes( label = macotra),hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Sulfonamides and trimethoprim, no Portugal data, linear model")+ 
  theme(strip.text = element_text(size = 12),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"), legend.position = "none")
ggsave("plots/fig_trim_noPort_linear.pdf", width = 25, height = 10)

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
ggsave("plots/all_data.pdf")

summary(all_data_model)

#######******** All data **********###########
# List of blactams
blacts = c("b-lactam",
           "penicillins_with_extended_spectrum",
           "blactamase_sensitive_penicillins","blactamase_resistant_penicillins",
           "blactamase_inhibitors","combination_penicillins_incl_blactamase_inhibitors",
           "other_blactams", 
           "1G_cephalosporins", "2G_cephalosporins", "3G_cephalosporins", 
           "4G_cephalosporin", "monobactams", "carbapenems", "other_cephalosporins_and_penems")

#######*#### Fit linear model to each drug across all settings
u <- unique(data_abx$sett_drug) # which drugs and settings? 
table_store<- c() # drug / intercept / value / p value / Rsquared

# linear
fit <- lmList(infection_per_100000 ~ value | sett_drug, data=data_abx)
r2 <- sapply(fit,function(x) summary(x)$r.squared)
s <- sapply(fit,function(x) tidy(x))
A <- sapply(fit,function(x) AIC(x))

for(i in u){
  if(i %in% colnames(s)){
    w_col = which(colnames(s) == i)
    table_store <- rbind(table_store, 
                         c(unlist(strsplit(i, ".", fixed = TRUE)), round(s[,w_col][[2]],4),signif(s[,w_col][[5]][2],3),round(r2[which(names(r2)==i)],4),signif(A[which(names(A)==i)],6),"linear"))
  }
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

# Different depending on prediction (AIC) / current fit (Rsquared) distinction
# Does an exponential or linear model fit better? 
best_aic %>% filter(pval_value < 0.05) %>% print(n=Inf) %>% ungroup() %>% summarise(unique(model))
best_rsq %>% filter(pval_value < 0.05) %>% print(n=Inf) %>% ungroup() %>% summarise(unique(model))
tt <- table_stored %>% arrange(drug, setting) %>% filter(pval_value < 0.05, setting == "Community") %>% 
  group_by(drug, setting) %>% mutate(n = n()) %>% filter(n > 1) %>% print(n=Inf) %>% 
  mutate(top = ifelse(Rsquared == max_R, 0, 1)) %>%
  select(setting, drug, Rsquared, model, top)

# Examples
eg_exp_lin <- left_join(data_abx %>% ungroup() %>% filter(setting == "Community", drug %in% c("2G_quinolones", "3G_cephalosporins", "b-lactam")),
          tt)

g1 <- ggplot(eg_exp_lin %>% filter(model == "linear"),
       aes(x=value, y = infection_per_100000)) + 
  geom_point() + 
  geom_label(y = 95, aes(color = factor(top), x = Inf, label = paste0("R^2 = ", Rsquared)),hjust = 1) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = y ~ x) +
  scale_x_continuous("Antibiotic usage (DDD per 1000 inhabitants and per day)") + 
  scale_y_continuous("Infection incidence per 100,000 inhabitants") + 
  scale_color_manual("Highest R^2",breaks = c(1,0), labels = c("yes","no"), values = c("red", "black")) + 
  ggtitle("Linear - trend line")
ggsave("plots/data_by_drug_eg_exponential_linear_g1.pdf", width = 5, height = 10)

g2 <- ggplot(eg_exp_lin %>% filter(model == "exponential"),
       aes(x=value, y = log(infection_per_100000))) + 
  geom_point() + 
  geom_label(y = 5.5, aes(color = factor(top), x = Inf, label = paste0("R^2 = ", Rsquared)),hjust = 1) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = y ~ x) +
  scale_x_continuous("Antibiotic usage (DDD per 1000 inhabitants and per day)") + 
  scale_y_continuous("Logarithm of Infection incidence per 100,000 inhabitants", lim = c(0,6)) + 
  expand_limits(x = 0, y = 0) + 
  scale_color_manual("Highest R^2",breaks = c(1,0), labels = c("yes","no"), values = c("red", "black")) + 
  ggtitle("Exponential - linear fit in algorithm to log(y)~x")
ggsave("plots/data_by_drug_eg_exponential_linear_g2.pdf", width = 5, height = 10)

g3 <- ggplot(eg_exp_lin %>% filter(model == "exponential"),
             aes(x=value, y = infection_per_100000)) + 
  geom_point() + 
  geom_label(y = 95, aes(color = factor(top), x = Inf, label = paste0("R^2 = ", Rsquared)),hjust = 1) + 
  facet_grid(setting~drug, scales = "free") + 
  geom_smooth(method='lm', formula = (y ~ exp(x))) +
  scale_x_continuous("Antibiotic usage (DDD per 1000 inhabitants and per day)") + 
  scale_y_continuous("Infection incidence per 100,000 inhabitants") + 
  expand_limits(x = 0, y = 0) + 
  scale_color_manual("Highest R^2",breaks = c(1,0), labels = c("yes","no"), values = c("red", "black")) + 
  ggtitle("Exponential - trend line")
ggsave("plots/data_by_drug_eg_exponential_linear_g3.pdf", width = 5, height = 10)

g1 / g2 / g3  + plot_layout(guides = 'collect') 
#g1 / g2  + plot_layout(guides = 'collect') 
ggsave("plots/data_by_drug_eg_exponential_linear.pdf", width = 10, height = 10)

### Export table
write.table(table_stored[order(table_stored$drug),c("setting","drug","model","intercept","value","Rsquared","AIC","pval_value","significant")], 
            file = "table_supp.txt", sep = ",", quote = FALSE, row.names = F)

#######******************* plot those significant ***************************########################
best_rsq$pval_value <- as.numeric(best_rsq$pval_value)
df_linear <- best_rsq %>% filter(pval_value < 0.05, model == "linear")
df_exponential <- best_rsq %>% filter(pval_value < 0.05, model == "exponential")

new_data_linear <- c()
for(i in 1:dim(df_linear)[1]){
  new_data_linear <- rbind(new_data_linear, data_abx %>% filter(drug == df_linear$drug[i] & setting == df_linear$setting[i]))
}

new_data_linear <- as.data.frame(new_data_linear)

g1 <- ggplot(new_data_linear %>% filter(!setting == "Community & Hospital"), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",sett_drug %in% unique(new_data_linear$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")
ggsave("plots/linear_just_significant_0.01.pdf", width = 22, height = 8)

g1_nb <- ggplot(new_data_linear %>% filter(!drug %in% blacts, !setting == "Community & Hospital"), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!drug %in% blacts, !setting == "Community & Hospital",sett_drug %in% unique(new_data_linear$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")

new_data_exponential <- c()
for(i in 1:dim(df_exponential)[1]){
  new_data_exponential <- rbind(new_data_exponential, data_abx %>% filter(drug == df_exponential$drug[i] & setting == df_exponential$setting[i]))
}

formula = y ~ x
new_data_exponential <- as.data.frame(new_data_exponential)
g2 <- ggplot(new_data_exponential %>% filter(!setting == "Community & Hospital"), aes(x=value, y = log_data, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",sett_drug %in% unique(new_data_exponential$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(new_data_exponential$log_data + 2)), "Log of MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")
ggsave("plots/exponential_just_significant_0.01.pdf", width = 22, height = 8)

g2_nb <- ggplot(new_data_exponential %>% filter(!drug %in% blacts, !setting == "Community & Hospital"), aes(x=value, y = log_data, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!drug %in% blacts, !setting == "Community & Hospital",sett_drug %in% unique(new_data_exponential$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(new_data_exponential$log_data + 2)), "Log of MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")

##### Remove outlier = Portugal ##########
#### Fit linear model to each drug across all settings
table_store<- c() # drug / intercept / value / p value / Rsquared

# linear
fit <- lmList(infection_per_100000 ~ value | sett_drug, data=data_abx %>% filter(!country == "Portugal"))
r2 <- sapply(fit,function(x) summary(x)$r.squared)
s <- sapply(fit,function(x) tidy(x))
A <- sapply(fit,function(x) AIC(x))

for(i in u){
  w_col = which(colnames(s) == i)
  table_store <- rbind(table_store, 
                       c(unlist(strsplit(i, ".", fixed = TRUE)), round(s[,w_col][[2]],4),signif(s[,w_col][[5]][2],3),round(r2[which(names(r2)==i)],4),signif(A[which(names(A)==i)],6),"linear"))
}

# exponential
fit <- lmList(log(infection_per_100000) ~ value | sett_drug, data=data_abx %>% filter(!country == "Portugal"))
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
            file = "table_supp_no_portugal.txt", sep = ",", quote = FALSE, row.names = F)

### plot those significant 
best_rsq$pval_value <- as.numeric(best_rsq$pval_value)
df_linear <- best_rsq %>% filter(pval_value < 0.05, model == "linear")
df_exponential <- best_rsq %>% filter(pval_value < 0.05, model == "exponential")

new_data_linear <- c()
for(i in 1:dim(df_linear)[1]){
  new_data_linear <- rbind(new_data_linear, data_abx %>% filter(drug == df_linear$drug[i] & setting == df_linear$setting[i]))
}
new_data_linear <- as.data.frame(new_data_linear)
g3 <- ggplot(new_data_linear %>% filter(!setting == "Community & Hospital",!country == "Portugal"), 
             aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",
                                        sett_drug %in% unique(new_data_linear$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")
ggsave("plots/linear_no_portugal_just_significant_0.01.pdf", width = 22, height = 8)

g3_nb <- ggplot(new_data_linear %>% filter(!country == "Portugal",!drug %in% blacts, !setting == "Community & Hospital"), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!drug %in% blacts, !setting == "Community & Hospital",sett_drug %in% unique(new_data_linear$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")

new_data_exponential <- c()
for(i in 1:dim(df_exponential)[1]){
  new_data_exponential <- rbind(new_data_exponential, data_abx %>% filter(drug == df_exponential$drug[i] & setting == df_exponential$setting[i]))
}

formula = y ~ x
new_data_exponential <- as.data.frame(new_data_exponential)


g4 <- ggplot(new_data_exponential %>% filter(!country == "Portugal", !setting == "Community & Hospital"), aes(x=value, y = log_data, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",sett_drug %in% unique(new_data_exponential$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(new_data_exponential$log_data + 2)), "Log of MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")
ggsave("plots/exponential_no_portugal_just_significant_0.01.pdf", width = 22, height = 8)

g4_nb <- ggplot(new_data_exponential %>% filter(!country == "Portugal", !drug %in% blacts, !setting == "Community & Hospital"), aes(x=value, y = log_data, label = macotra, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(size = summary_class)) +
  geom_point(data = data_abx %>% filter(!drug %in% blacts, !setting == "Community & Hospital",sett_drug %in% unique(new_data_exponential$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(new_data_exponential$log_data + 2)), "Log of MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class") + 
  guides(fill="none", col = "none")



### Combine for supplementary
g1 / g3 + plot_annotation(tag_levels = "A")
ggsave("plots/linear_both_just_significant_0.01_all.pdf", width = 12, height = 8)

g2 / g4 + plot_annotation(tag_levels = "A")
ggsave("plots/exponential_both_just_significant_0.01_all.pdf", width = 22, height = 16)

g1_nb / g3_nb + plot_annotation(tag_levels = "A")
ggsave("plots/linear_both_just_significant_0.01_nb.pdf", width = 12, height = 8)

g2_nb / g4_nb + plot_annotation(tag_levels = "A")
ggsave("plots/exponential_both_just_significant_0.01_nb.pdf", width = 22, height = 16)

####### FIGURE 
g3 <- ggplot(new_data_linear %>% filter(!setting == "Community & Hospital"), aes(x=value, y = infection_per_100000, label = macotra, group = setting)) + 
  geom_point(aes(col = factor(summary_class))) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes(fill = factor(summary_class))) + 
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",sett_drug %in% unique(new_data_linear$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000 + 10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  guides(fill="none", col = "none")

g4 <- ggplot(new_data_exponential %>% filter(!setting == "Community & Hospital"), aes(x=value, y = log_data, label = macotra, group = setting)) + 
  geom_point(aes(col = factor(summary_class))) + 
  facet_nested(setting~ATC_code_family + ATC_code, scales = "free") + #facet_grid(setting ~ drug, scales = "free",labeller = labeller(facet_category = label_wrap_gen(width = 16))) + 
  geom_smooth(method='lm', formula = formula, aes(fill = factor(summary_class))) + 
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",sett_drug %in% unique(new_data_exponential$sett_drug), country %in% c("UK","Netherlands","France")), pch = 3) + 
  geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label),
                                  after_stat(p.value.label),
                                  sep = "*\", \"*")), formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) + 
  scale_y_continuous(lim = c(0,max(new_data_exponential$log_data + 10)), "Log of MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  guides(fill="none", col = "none")

g3 / g4 + plot_annotation(tag_levels = "A")
ggsave("plots/both_just_significant_0.01_linear_and_exp.pdf", width = 22, height = 16)

# ##### Figure 5 for paper 
formula = y ~ x
data_abx$ATC_code <- factor(data_abx$ATC_code, c("J01", "J01A", "J01B", "J01C_D", "J01C", "J01CA", 
                                                 "J01CF", "J01CG", "J01CR", "J01D", "J01DB", "J01DC", "J01DD", 
                                                 "J01CE", 
                                                 "J01DE", "J01DF", "J01DH", "J01DI", "J01E", "J01EA", "J01EB", 
                                                 "J01EC", "J01EE", "J01F", 
                                                 "J01FA", "J01FA_shrt", "J01FA_intr", "J01FA_long", 
                                                 "J01FF", "J01FG", "J01G", "J01GB", "J01M", "J01MA", "J01M_frst", 
                                                 "J01M_scnd", "J01M_thrd", "J01MB", "J01XA", "J01XB", 
                                                 "J01XC", "J01XD", "J01XE"))

# New facet label names for dose variable
dose.labs <- c("Total\nbeta-lactams",  
               "Combination\npenicillins\nincluding\nbeta-lactamase\ninhibitors","3rd generation\ncepahlosporins",
               "Beta-lactamase\nsensitive\npenicillins",
               "Fluoroquinolones","Macrolides","Aminoglycosides")
names(dose.labs) <- c("all_b_lactam", 
                      "combination_penicillins_incl_blactamase_inhibitors","3G_cephalosporins",
                      "blactamase_sensitive_penicillins",
                      "fluoroquinolones","macrolides","other_aminoglycosides")

g5a<- ggplot(data_abx %>% filter(ATC_code %in% c("J01C_D","J01DD","J01CR","J01CE")), 
       aes(x=value, y = infection_per_100000, group = setting)) +
  geom_rect(aes(fill = summary_class),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha=0.1) + 
  scale_fill_manual("",values = c("white","grey90")) + 
  geom_point(aes()) + 
  facet_grid(drug~setting, scales = "free",labeller = labeller(drug = dose.labs)) + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code %in% c("J01C_D","J01DD","J01CR","J01CE")), pch = 3) + 
  geom_text(aes(label = macotra),nudge_x = 0.4,nudge_y=0.1) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000+10)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Beta-lactams") + 
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, angle = 0),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
          plot.title = element_text(size = 12), 
        panel.grid.major=element_line(colour="black"),
        panel.grid.minor=element_line(colour="black"),
        legend.position = "none")
ggsave("plots/fig5a.pdf")

g5b<- ggplot(data_abx %>% filter(!country == "Portugal",ATC_code %in% c("J01MA","J01FA","J01GB")), 
             aes(x=value, y = infection_per_100000, group = setting)) +
  geom_point(aes()) + 
  facet_grid(drug~setting, scales = "free",labeller = labeller(drug = dose.labs)) + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(!country == "Portugal",ATC_code %in% c("J01MA","J01FA","J01GB")), pch = 3) + 
  geom_text(aes(label = macotra),hjust=0.27,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max(data_abx$infection_per_100000+5)), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  ggtitle("Other antibiotics") + 
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, angle = 0),
        axis.text =element_text(size=12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12), 
        legend.position = "none")
ggsave("plots/fig5b.pdf")

g5a + g5b +
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths = c(2, 2))
ggsave("plots/Figure5.pdf", width = 16, height = 8)

ggplot(data_abx %>% filter(ATC_code_family %in% c("J01E"), !country == "Portugal"), 
       aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01E"),!country == "Portugal"), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")

g5a<- ggplot(data_abx %>% filter(ATC_code_family %in% c("J01C","J01D"),
                                 !ATC_code == "J01CG"), 
             aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01C","J01D"),
                      !ATC_code == "J01CG"), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")
ggsave("plots/fig5a_allblactams.pdf")

g5a_egs<- ggplot(data_abx %>% filter(!setting == "Community & Hospital", 
                 ATC_code_family %in% c("J01C","J01D"),
                                 ATC_code %in% c("J01CE","J01CR",
                                                 "J01DD")), 
             aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(!setting == "Community & Hospital",
             country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family == "J01C",
                      ATC_code %in% c(
                                      "J01CE","J01CR",
                                      "J01DD")), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")
ggsave("plots/fig5a_selective.pdf")

g5b_all<- ggplot(data_abx %>% filter(ATC_code_family %in% c("J01F","J01M")), 
             aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01F","J01M")), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")
ggsave("plots/fig5b_all.pdf")

g5b_quin<- ggplot(data_abx %>% filter(ATC_code_family %in% c("J01F","J01M"),
                                 ATC_code %in% c("J01M","J01MA","J01M_scnd")), 
             aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01F","J01M"),
                      ATC_code %in% c("J01M","J01MA","J01M_scnd")), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")
ggsave("plots/fig5b_quin_selective.pdf")

g5b_mac_all<- ggplot(data_abx %>% filter(ATC_code_family %in% c("J01F"),
                                 !country %in% c("Portugal")), 
             aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01F"), !country %in% c("Portugal")), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")
ggsave("plots/fig5b_mac_all.pdf")

g5b_mac_selective<- ggplot(data_abx %>% filter(ATC_code_family %in% c("J01F"),
                                 !country %in% c("Portugal"), 
                                 ATC_code %in% c("J01F", "J01FA","J01FA_intr",
                                                 "J01FA_long","J01FA_shrt")), 
             aes(x=value, y = infection_per_100000, group = setting)) + 
  geom_point(aes()) + 
  facet_nested(ATC_code_family + ATC_code ~ setting, scales = "free") + 
  geom_smooth(method='lm', formula = formula, aes()) +
  geom_point(data = data_abx %>% filter(country %in% c("UK","Netherlands","France"))%>% 
               filter(ATC_code_family %in% c("J01F"),
                      ATC_code %in% c("J01F", "J01FA","J01FA_intr",
                                      "J01FA_long","J01FA_shrt")), col ="blue",pch = 3) + 
  #geom_text(hjust=-0.17,vjust=0) + 
  stat_poly_eq(mapping = aes(label =  paste(after_stat(rr.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*"),
                             color = ifelse(after_stat(p.value) < 0.05, "red", "black")), 
               formula = formula,label.x.npc = 0.9, label.y.npc = 0.9) +
  scale_color_identity() + 
  scale_y_continuous(lim = c(0,max((data_abx$infection_per_100000 + 10))), "MRSA incidence per 100,000") + 
  scale_x_continuous(lim = c(0,NA), "Total usage (DDD 1000 inhabitants and per day)") + 
  scale_size_manual(values = c(1, 1.5), "Summary class")
ggsave("plots/fig5b_mac_selective.pdf")


g5a_egs + (g5b_quin / g5b_mac_selective) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths = c(1, 1))
ggsave("plots/fig5.pdf", width = 25, height = 10)

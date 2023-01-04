### Data manipulation and exploration of EARS-Net data on antibitoic usage
# https://www.ecdc.europa.eu/en/antimicrobial-consumption/surveillance-and-disease-data/database

## Needed libraries
library(tidyverse)
library(naniar)
library(here)
library(patchwork)
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
  mutate(summary_class = ifelse(drug %in% c("b-lactam","other_blactams",
                                            "sulfonamides_and_trimethoprim","macrolides_lincosamides_streptogramins",
                                            "Quinolones","macrolides","fluoroquinolones"),
                                "Yes","No"))



data_abx$drug <- factor(data_abx$drug, levels = c("antibacterials_for_systemic", 
                                                  "tetracyclines", 
                                                  "amphenicols", 
                                                  "b-lactam",
                                                  "penicillins_with_extended_spectrum", "blactamase_sensitive_penicillins", 
                                                  "blactamase_resistant_penicillins", "blactamase_inhibitors", 
                                                  "combination_penicillins_incl_blactamase_inhibitors", 
                                                  "other_blactams", 
                                                  "1G_cephalosporins", "2G_cephalosporins", "3G_cephalosporins", 
                                                  "4G_cephalosporin", "monobactams", "carbapenems", "other_cephalosporins_and_penems", 
                                                  "Quinolones", 
                                                  "1G_quinolones", "2G_quinolones", "3G_quinolones", 
                                                  "fluoroquinolones", "other_quinolones", 
                                                  "sulfonamides_and_trimethoprim", 
                                                  "trimethoprim_and_derivatives", "short_acting_sulfonamides", "intermediate_acting_sulfonamides", 
                                                  "sulfonamides_and_trimethoprim_combs",
                                                  "macrolides_lincosamides_streptogramins", 
                                                  "macrolides", "lincosamides", "streptogramins", 
                                                  "short_acting_macrolides", "intermediate_acting_macrolides", "long_acting_macrolides", 
                                                  "aminoglycosides", 
                                                  "other_aminoglycosides", 
                                                  "glycopeptides", "polymyxins", "steroid_antibacterials", "imadazole", "nitrofuran"))

##******************************** Exploration of trends in antibiotic use ***********************************#####################
## Variation in drug use 
ggplot(data_abx %>% filter(!drug == "antibacterials_for_systemic"), 
       aes(x = country, y = value)) + geom_point(aes(col = drug)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")
ggsave("plots/drug_variation_setting_country.pdf", width = 11, height = 6)

totals_abx = data_abx %>% filter(!drug == "antibacterials_for_systemic") %>% group_by(country,setting) %>% summarise(total = sum(value))
data_abx_t <- left_join(data_abx, totals_abx %>% filter(setting == "Community & Hospital") %>% select(country, total))
write.csv(data_abx_t, "data/data_abx_with_totals.csv")

ggplot(data_abx_t %>% filter(!drug == "antibacterials_for_systemic", setting == "Community & Hospital"), 
       aes(x = reorder(country,total), y = value)) + 
  geom_bar(stat = "identity", position = "stack", aes(x = reorder(country,total), fill = ATC_code_family)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_fill_manual("ATC code family", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')) + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")
ggsave("plots/drug_variation_setting_country_bar.pdf", width = 10, height = 7)

ggplot(data_abx %>% filter(!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors"), 
       aes(x = country, y = value)) + geom_point(aes(col = ATC_code)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")
ggsave("plots/ATC_code_variation_setting_country.pdf", width = 10, height = 6)


ggplot(data_abx %>% filter(summary_class == "No",!ATC_code == "JJ", !drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors"), 
       aes(x = country, y = value)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = drug)) + 
  facet_grid(setting ~ ATC_code_family, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() 
ggsave("plots/Comm_hosp_variation_fillbar_country.pdf", width = 25, height = 15)

ggplot(data_abx %>% filter(summary_class == "No",!ATC_code == "JJ", !drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors") %>% 
         group_by(country, setting) %>% summarise(total = sum(value)), 
       aes(x = country, y = total)) + 
  geom_point(aes()) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Total DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() 
ggsave("plots/Comm_hosp_variation_total_country.pdf", width = 10, height = 6)


ggplot(data_abx %>% filter(summary_class == "No",!ATC_code == "JJ", !macotra == "",!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors"), 
       aes(x = country, y = value)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = drug)) + 
  facet_grid(setting ~ ATC_code_family, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() 
ggsave("plots/Comm_hosp_variation_fillbar_country_macotra.pdf", width = 20, height = 6)

g2 <- ggplot(data_abx %>% filter(summary_class == "No",!ATC_code == "JJ", !macotra == "", !drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors") %>% 
         group_by(country, setting) %>% summarise(total = sum(value)), 
       aes(x = country, y = total)) + 
  geom_point(aes()) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Total DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() 
ggsave("plots/Comm_hosp_variation_total_country_mactora.pdf", width = 10, height = 3)



# g1 / g2
# ggsave("plots/Comm_hosp_variation_total_country_mactora_both.pdf", width = 10, height = 3)

## Highest use 
data_abx %>% filter(!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors") %>% 
  arrange(desc(value)) %>% print(n = 200)

data_abx %>% filter(!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors", setting == "Community") %>% 
  arrange(desc(value)) %>% print(n = 20)

data_abx %>% filter(!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors", setting == "Hospital") %>% 
  arrange(desc(value)) %>% print(n = 40)

## Beta-lactams
# J01C BETA-LACTAM ANTIBACTERIALS, PENICILLINS
# J01CA Penicillins with extended spectrum
# J01CE Beta-lactamase sensitive penicillins
# J01CF Beta-lactamase resistant penicillins
# J01CG Beta-lactamase inhibitors
# J01CR Combinations of penicillins, incl. beta-lactamase inhibitors

g1 <- ggplot(data_abx %>% filter(!drug == "antibacterials_for_systemic", 
                                 ATC_code %in% c("J01CA", "J01CE", "J01CF", "J01CG","J01CR")), 
             aes(x = country, y = value)) + geom_point(aes(col = drug)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom") + 
  ggtitle("J01C Beta-lactams")

## Other beta-lactams
# J01DB First-generation cephalosporins
# J01DC Second-generation cephalosporins
# J01DD Third-generation cephalosporins
# J01DE Fourth-generation cephalosporins
# J01DF Monobactams
# J01DH Carbapenems
# J01DI Other cephalosporins and penems
g2 <- ggplot(data_abx %>% filter(!drug == "antibacterials_for_systemic", 
                                 ATC_code %in% c("J01DB","J01DC","J01DD","J01DE","J01DF","J01DH","J01DI")), 
             aes(x = country, y = value)) + geom_point(aes(col = drug)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom") + 
  ggtitle("J01D Other Beta-lactams")
g1 / g2 
ggsave("plots/beta_lactam_variation_setting_country.pdf", width = 15, height = 10)

## Quinolones
g1 <- ggplot(data_abx %>% filter(drug %in% c("Quinolones", "1G_quinolones", "2G_quinolones", "3G_quinolones", "fluoroquinolones", "other_quinolones")), 
             aes(x = country, y = value)) + geom_point(aes(col = drug)) + 
  geom_point(data = data_abx %>% filter(drug == "Quinolones"), pch = "x", size = 3) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom") + 
  ggtitle("J01M Quinolones")
ggsave("plots/quinolones_variation_setting_country.pdf", width = 10, height = 8)


#### Explore ranges
data_abx %>% group_by(drug,setting) %>% summarise(min = min(value), max = max(value)) %>% arrange(desc(max))
data_abx %>% group_by(drug,setting) %>% summarise(min = min(value), max = max(value)) %>% arrange(desc(-min)) %>% print(n=30)

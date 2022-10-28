### Data manipulation and exploration of EARS-Net data on antibitoic usage
# https://www.ecdc.europa.eu/en/antimicrobial-consumption/surveillance-and-disease-data/database

## Needed libraries
library(tidyverse)
library(naniar)
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

##******************************** Exploration of trends in antibiotic use ***********************************#####################
## Variation in drug use 
ggplot(data_abx %>% filter(!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors"), 
       aes(x = country, y = value)) + geom_point(aes(col = drug)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("Drug") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")
ggsave("plots/drug_variation_setting_country.pdf", width = 15, height = 15)

ggplot(data_abx %>% filter(!drug == "antibacterials_for_systemic", !drug == "combination_penicillins_incl_blactamase_inhibitors"), 
       aes(x = country, y = value)) + geom_point(aes(col = ATC_code)) + 
  facet_wrap(~setting, scales = "free") + 
  scale_colour_discrete("ATC_code") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")
ggsave("plots/ATC_code_variation_setting_country.pdf", width = 15, height = 15)

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

### Data manipulation and exploration of Cassini data on infection 
# Taken from Country Specific Results in Cassini supplementary 

library(tidyverse)
theme_set(theme_bw(base_size = 11))
library(here)
library(patchwork)
setwd(here())

## Data
data_cass_orig <- read_csv("data/cassini_supp_mrsa_incidence.csv")

##******************************** Exploration of trends in infection ***********************************#####################
g1 <- ggplot(data_cass_orig, aes(x=country, y = inc_median, group = syndrome)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = syndrome)) + 
  scale_fill_discrete("Infection type") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Median incidence per 100,000") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() 

g2 <- ggplot(data_cass_orig, aes(x=country, y = inc_median, group = syndrome)) + 
  geom_bar(stat = "identity", position = "fill", aes(fill = syndrome)) + 
  scale_fill_discrete("Infection type") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Proportion of median incidence per 100,000") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() 

g1 + g2 + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
ggsave("plots/mrsa_incidence_by_syndrome_props.pdf")

##******************************** Country totals ***********************************#####################
data_cass <- data_cass_orig %>% group_by(country) %>% summarise(tot_inc_low = sum(inc_low),
                                              total_inc_median = sum(inc_median),
                                              total_inc_high = sum(inc_high))

ggplot(data_cass, aes(x=country,y = total_inc_median)) + geom_point() + 
  geom_errorbar(aes(ymin = tot_inc_low, ymax = total_inc_high)) + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Incidence per 100,000")
ggsave("plots/mrsa_incidence_by_country.pdf")

## Data to use 
write.csv(data_cass %>% select(country, total_inc_median) %>% rename("infection_per_100000" = "total_inc_median"), 
          "data/cass_infection_output.csv")
## cut and paste into file with antibiotic use - first column of data_amr_abx_use_macotra.csv


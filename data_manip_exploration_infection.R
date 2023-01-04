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

## Figure for paper: 
data_cass_all <- left_join(data_cass_orig, data_cass)
g1 <- ggplot(data_cass_all %>% filter(!country == "Iceland"), aes(x=reorder(country,total_inc_median), y = inc_median)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = syndrome)) + 
  geom_point(aes(x = reorder(country,total_inc_median), y = total_inc_median)) + 
  geom_errorbar(aes(x = reorder(country,total_inc_median), ymin = tot_inc_low, ymax = total_inc_high)) + 
  scale_fill_discrete("Infection type") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Median incidence per 100,000 inhabitants") + 
  theme(strip.text.y = element_text(angle = 0),legend.position="bottom") + 
  coord_flip() 
ggsave("plots/Fig4_mrsa_incidence_by_country.pdf", width = 10, height = 7)

# Add in antibiotic usage
data_abx_t <- read_csv("data/data_abx_with_totals.csv")

g2 <- ggplot(data_abx_t %>% filter(!drug == "antibacterials_for_systemic", setting == "Community & Hospital"), 
       aes(x = reorder(country,total), y = value)) + 
  geom_bar(stat = "identity", position = "stack", aes(x = reorder(country,total), fill = ATC_code_family)) + 
  scale_fill_manual("ATC code family", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')) + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")
ggsave("plots/drug_variation_setting_country_bar.pdf", width = 10, height = 7)

g1 + g2 + plot_layout(widths = c(1.5,1))
ggsave("plots/Fig4_proposed2.pdf", width = 20, height = 8)

fig4_data <- left_join(data_abx_t, data_cass_all)
g3 <- ggplot(fig4_data %>% filter(!drug == "antibacterials_for_systemic", setting == "Community & Hospital"), 
             aes(x = reorder(country,total_inc_median), y = value)) + 
  geom_bar(stat = "identity", position = "stack", aes(x = reorder(country,total_inc_median), fill = ATC_code_family)) + 
  scale_fill_manual("ATC code family", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')) + 
  scale_x_discrete("Country") + 
  scale_y_continuous("DDD per 1000 inhabitants and per day") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_flip() + theme(legend.position="bottom")

g1 + g3 + plot_layout(widths = c(1.5,1))
ggsave("plots/Fig4_proposed2.pdf", width = 20, height = 8)

# Split out low incidence 
g1a <- ggplot(data_cass_all %>% filter(!country == "Iceland", total_inc_median > 7), aes(x=reorder(country,total_inc_median), y = inc_median)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = syndrome)) + 
  geom_point(aes(x = reorder(country,total_inc_median), y = total_inc_median)) + 
  geom_errorbar(aes(x = reorder(country,total_inc_median), ymin = tot_inc_low, ymax = total_inc_high)) + 
  scale_fill_discrete("Infection type") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Median incidence per 100,000 inhabitants") + 
  theme(strip.text.y = element_text(angle = 0),legend.position="bottom") + 
  coord_flip() 

g1b <- ggplot(data_cass_all %>% filter(!country == "Iceland", total_inc_median <= 7), aes(x=reorder(country,total_inc_median), y = inc_median)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = syndrome)) + 
  geom_point(aes(x = reorder(country,total_inc_median), y = total_inc_median)) + 
  geom_errorbar(aes(x = reorder(country,total_inc_median), ymin = tot_inc_low, ymax = total_inc_high)) + 
  scale_fill_discrete("Infection type") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Median incidence per 100,000 inhabitants") + 
  theme(strip.text.y = element_text(angle = 0),legend.position="bottom") + 
  coord_flip() 

((g1a + g1b + plot_layout(heights =  c(5,1), guides = "collect") & theme(legend.position = 'bottom'))| g3) + plot_layout(widths = c(1.5,1)) & theme(legend.position = 'bottom')
ggsave("plots/Fig4_proposed3.pdf", width = 20, height = 8)

 ## Data to use 
write.csv(data_cass %>% select(country, total_inc_median) %>% rename("infection_per_100000" = "total_inc_median"), 
          "data/cass_infection_output.csv")
## cut and paste into file with antibiotic use - first column of data_amr_abx_use_macotra.csv

# Ranges on data
data_cass %>% filter(total_inc_median == min(data_cass$total_inc_median, na.rm=TRUE))
data_cass %>% filter(total_inc_median == max(data_cass$total_inc_median, na.rm=TRUE))


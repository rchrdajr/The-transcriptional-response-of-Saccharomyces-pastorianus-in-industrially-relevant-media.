library(tidyverse)
#library(forcats)

ethGO <- read.csv("gProfiler_eth_intersections-Copy.csv")
leuGO <- read.csv("gProfiler_leu_intersections-Copy.csv")
wortGO <- read.csv("gProfiler_wort_intersections-Copy.csv")

head(ethGO)
str(ethGO)

as.double(ethGO$negative_log10_of_adjusted_p_value)
#------- DATA WRANGLING -------

ethGO_clean <- ethGO %>%
  filter(source == "GO:MF" | source == "GO:CC" | source == "GO:BP")

ethGO_clean <- ethGO_clean %>%
  mutate(source = recode(source, "GO:MF" = "Molecular Function",
                                 "GO:CC" = "Cellular "))

leuGO_clean <- leuGO %>%
  filter(source == "GO:MF" | source == "GO:CC" | source == "GO:BP")

wortGO_clean <- wortGO %>%
  filter(source == "GO:MF" | source == "GO:CC" | source == "GO:BP")


#------- HIST PLOTS -----------
ethGO_clean %>%
  ggplot(aes(x = negative_log10_of_adjusted_p_value, y = fct_reorder(term_name, source))) +
  geom_col(aes(fill = source)) +
  theme_classic() +
  labs(x = "-log 10 (p-value)", 
       y = "",
       title = "SD vs SD ethanol")

leuGO_clean %>%
  ggplot(aes(x = negative_log10_of_adjusted_p_value, y = fct_reorder(term_name, source))) +
  geom_col(aes(fill = source)) +
  theme_classic() +
  labs(x = "-log 10 (p-value)", 
       y = "",
       title = "SD vs SD - Leucine")

wortGO_clean %>%
  ggplot(aes(x = negative_log10_of_adjusted_p_value, y = fct_reorder(term_name, source))) +
  geom_col(aes(fill = source)) +
  theme_classic() +
  labs(x = "-log 10 (p-value)", 
       y = "",
       title = "SD vs Wort")

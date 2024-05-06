library(tidyverse)
#library(forcats)

ethGO <- read.csv("gProfiler_eth_new.csv")
leuGO <- read.csv("gProfiler_leu_new.csv")
wortGO <- read.csv("gProfiler_wort_new.csv")

head(ethGO)
str(ethGO)

as.double(ethGO$negative_log10_of_adjusted_p_value)
#------- DATA WRANGLING -------

ethGO_clean <- ethGO %>%
  filter(source == "GO:MF" | source == "GO:CC" | source == "GO:BP")

ethGO_clean <- ethGO_clean %>%
  mutate(source = recode(source, "GO:MF" = "Molecular Function",
                         "GO:CC" = "Cellular Component",
                         "GO:BP" = "Biological Process"))

leuGO_clean <- leuGO %>%
  filter(source == "GO:MF" | source == "GO:CC" | source == "GO:BP")

leuGO_clean <- leuGO_clean %>%
  mutate(source = recode(source, "GO:MF" = "Molecular Function",
                         "GO:CC" = "Cellular Component",
                         "GO:BP" = "Biological Process"))

wortGO_clean <- wortGO %>%
  filter(source == "GO:MF" | source == "GO:CC" | source == "GO:BP")

wortGO_clean <- wortGO_clean %>%
  mutate(source = recode(source, "GO:MF" = "Molecular Function",
                         "GO:CC" = "Cellular Component",
                         "GO:BP" = "Biological Process"))


#------- HIST PLOTS -----------

#png(file = "ethGO_hist_new.png", width = 500, height = 550)
ethhist <- ethGO_clean %>%
  ggplot(aes(x = negative_log10_of_adjusted_p_value, y = fct_reorder(term_name, source))) +
  geom_col(aes(fill = source)) +
  theme_classic() +
  labs(x = "-log 10 (p-value)", 
       y = "",
       title = "SD vs SD ethanol")
#dev.off()



png(file = "leuGO_hist_new.png", width = 500, height = 550)
leuhist <- leuGO_clean %>%
  ggplot(aes(x = negative_log10_of_adjusted_p_value, y = fct_reorder(term_name, source))) +
  geom_col(aes(fill = source)) +
  theme_classic() +
  labs(x = "-log 10 (p-value)", 
       y = "",
       title = "SD vs SD - Leucine",)
dev.off()

png(file = "wortGO_hist_new.png", width = 500, height = 550)
worthist <- wortGO_clean %>%
  ggplot(aes(x = negative_log10_of_adjusted_p_value, y = fct_reorder(term_name, source))) +
  geom_col(aes(fill = source)) +
  theme_classic() +
  labs(x = "-log 10 (p-value)", 
       y = "",
       title = "SD vs Wort")
dev.off()

allGO_hist <- (ethhist + leuhist + worthist)
ggsave("allGO_hist.png", plot = allGO_hist, width = 15, height = 5)

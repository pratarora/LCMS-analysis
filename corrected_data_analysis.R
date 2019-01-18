metadata <- read.csv("sample_metadata.csv", header = TRUE)
compound_rt_database <-
  read.csv("compound_rt_database.csv", header = TRUE)
maven_data <-
  read.csv("Maven_corrected_for_analysis.csv", header = TRUE)

library(dplyr)
library(ggplot2)
library(ggbiplot)
library(tidyr)
library(reshape2)


pool_total <- maven_data %>% select(-c(note)) %>%
  group_by(compound) %>% summarize_all(funs(sum))
metadata$Sample <- paste0("X", metadata$Sample)
df2 <- maven_data %>%
  group_by(compound, note)

df3 <- melt(df2, id = c("compound", "note"))

names(df3)[names(df3) == "variable"] <- "Sample"

df4 <- full_join(df3, metadata)

df5 <- pool_total %>% ungroup() %>%
  group_by(compound)
df6 <- melt(df5, id = c("compound"))
names(df6)[names(df6) == "variable"] <- "Sample"
names(df6)[names(df6) == "value"] <- "Total"
df7 <- full_join(df6, metadata)


fraction_enrichment <- full_join(df4, df7)
fraction_enrichment <-
  fraction_enrichment %>% mutate("fraction_enrich" = (value / Total) * 100)


totalgraph <- ggplot(
  fraction_enrichment,
  aes(
    x = Time..mins.,
    y = Total,
    group = interaction(compound, Phenotype),
    color = Phenotype
  )
) +
  geom_line() +
  labs(title = "Pool Total",
       x = "Time Points (in mins)",
       y = "Metabolite Levels") +
  facet_wrap(. ~ compound)
print(totalgraph)
dev.copy(pdf,
         file = "correctedtotalgraph.pdf",
         width = 5,
         height = 5)
dev.off ()


fractiongraph <- ggplot(
  fraction_enrichment,
  aes(
    x = Time..mins.,
    y = fraction_enrich,
    group = interaction(note, Phenotype),
    color = note
  )
) +
  geom_line() +
  labs(title = "Fraction Enrichment",
       x = "Time Points (in mins)",
       y = "Fraction Enrichment") +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ compound + Phenotype)

print(fractiongraph)
dev.copy(pdf,
         file = "correctedfractiongraph.pdf",
         width = 8,
         height = 8)
dev.off ()

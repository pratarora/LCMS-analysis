rm(list=ls())
options(digits=9)

setwd("C:/Users/prata/Desktop/Elucidata/LCMS analysis/LCMS-analysis/")
getwd()
metadata <- read.csv("sample_metadata.csv", header = TRUE)
maven_data <- read.csv("Maven_processed.csv", header=TRUE)

library(dplyr)
library(ggplot2)
library(ggbiplot)
library(tidyr)
library(reshape2)


pool_total <- maven_data %>% select(-c(label,metaGroupId,groupId,goodPeakCount,
                                       medMz,medRt,maxQuality,expectedRtDiff,
                                       ppmDiff,parent,note,formula,compoundId)) %>%
  group_by(compound) %>% summarize_all(funs(sum))

# Correction for R reading of numerical as factors
metadata$Sample <- paste0("X",metadata$Sample)

df2 <- maven_data %>% select(-c(label,metaGroupId,groupId,goodPeakCount,
                                medMz,medRt,maxQuality,expectedRtDiff,
                                ppmDiff,parent, compoundId, formula)) %>% 
  group_by(compound,note)

df3 <- melt(df2,id=c("compound","note"))

names(df3)[names(df3)=="variable"] <- "Sample"

df4 <- full_join(df3,metadata)

df5 <- pool_total %>% ungroup() %>% 
  group_by(compound)
df6 <- melt(df5,id=c("compound"))
names(df6)[names(df6)=="variable"] <- "Sample"
names(df6)[names(df6)=="value"] <- "Total"
df7 <- full_join(df6,metadata)


fraction_enrichment <- full_join(df4,df7)
fraction_enrichment <- fraction_enrichment %>% mutate("fraction_enrich"=(value/Total)*100)
write.csv(fraction_enrichment,
          file="Data with Pool and fractions.csv")




# PCA non transformed raw data
for (i in 2:ncol(pool_total)) {
  p <- ggplot(pool_total,
              aes(x=unlist(pool_total[,i],use.names=FALSE)))+
    geom_histogram(bins=10)
  print(p)
}

pca <- prcomp(pool_total[,c(2:ncol(pool_total))])
summary(pca)
ggbiplot(pca)



# PCA normalized data
pool_total_normalized <- pool_total
for (i in 2:ncol(pool_total)) {
  pool_total_normalized[,i] <- scale(pool_total[,i],center = TRUE) %>% as.vector
  
}


pca.normalized <- prcomp(pool_total_normalized[,c(2:ncol(pool_total))])
summary(pca.normalized)
ggbiplot(pca.normalized)

# PCA log transformed data
pool_total_log <- pool_total
for (i in 2:ncol(pool_total)) {
  pool_total_log[,i] <- log10(pool_total[,i]+1) %>% as.vector
  
}

for (i in 4:ncol(pool_total)) {
  p <- ggplot(pool_total_log,
              aes(x=unlist(pool_total_log[,i],use.names=FALSE)))+geom_histogram(bins=10)
  print(p)
}

pca.labels <- pool_total_log$compound
pca.log <- prcomp(pool_total_log[,c(2:ncol(pool_total))])

pca.log$sdev
screeplot(pca, type="lines",col=3)
pca$rotation
summary(pca.log)

ggbiplot(pca.log, ellipse=TRUE, labels=pca.labels)
dev.copy(pdf,
         file="PCAlogPC1PC2.pdf",
         height= 8,
         width=12)
dev.off()

ggbiplot(pca.log,choices = c(2,3), labels=pca.labels)
dev.copy(pdf,
         file="PCAlogPC2PC3.pdf",
         height= 8,
         width=12)
dev.off()

ggbiplot(pca.log,choices = c(3,4), labels=pca.labels)
dev.copy(pdf,
         file="PCAlogPC3PC4.pdf",
         height= 8,
         width=12)
dev.off()

pool_total_T <- add_rownames(pool_total) %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) 

names(pool_total_T)[names(pool_total_T)=="var"] <- "Sample"
colnames(pool_total_T) <- unlist(pool_total_T[1,],use.names=FALSE)
pool_total_T <- pool_total_T[-1,]

pool_total_T[,2:ncol(pool_total_T)] <- sapply(pool_total_T[,2:ncol(pool_total_T)], as.double)
str(pool_total_T)



for (i in 2:ncol(pool_total_T)) {
  p <- ggplot(pool_total_T,
              aes(x=unlist(pool_total_T[,i],use.names=FALSE)))+
    geom_histogram(bins=10)
  print(p)
}

pool_total_T_log <-pool_total_T 

for (i in 2:ncol(pool_total_T_log)) {
  pool_total_T_log[,i] <- log10(pool_total_T[,i]+1) %>% as.vector
  
}

for (i in 2:ncol(pool_total_T_log)) {
  p <- ggplot(pool_total_T_log,
              aes(x=unlist(pool_total_T_log[,i],use.names=FALSE)))+
    geom_histogram(bins=10)
  print(p)
}

pca.labels_T <- pool_total_T_log$compound
pca.log_T <- prcomp(pool_total_T_log[,c(2:ncol(pool_total))])

pca.log_T$sdev
screeplot(pca.log_T, type="lines",col=3)
pca.log_T$rotation
summary(pca.log_T)

ggbiplot(pca.log_T, ellipse=TRUE, labels=pca.labels_T)
dev.copy(pdf,
         file="PCAlogTPC1PC2.pdf",
         height= 8,
         width=12)
dev.off()


totalgraph <- ggplot(fraction_enrichment,
                     aes(x=Time..mins.,
                         y=Total,
                         group=interaction(compound,Phenotype),
                         color=Phenotype
                     ))+
  geom_line()+
  labs(
    title="Pool Total",
    x="Time Points (in mins)",
    y="Metabolite Levels"
  )+
  facet_wrap(.~compound, scales="free_y")
print(totalgraph)
dev.copy(
  pdf,
  file = "totalgraph.pdf",
  width = 25,
  height = 25
)
dev.off ()


fractiongraph <- ggplot(fraction_enrichment,
                        aes(x=Time..mins.,
                            y=fraction_enrich,
                            group=interaction(note,Phenotype),
                            color=note
                        ))+
  geom_line()+
  labs(
    title="Fraction Enrichment",
    x="Time Points (in mins)",
    y="Fraction Enrichment"
  )+
  theme(legend.position="bottom")+
  facet_wrap(.~compound+Phenotype, scales="free_y")

print(fractiongraph)
dev.copy(
  pdf,
  file = "fractiongraph.pdf",
  width = 25,
  height = 25
)
dev.off ()

library(accucor)

corrected_maven_data <- natural_abundance_correction(
  path="Maven_processed_2.csv", resolution = 100000,
  output_filetype='csv',
  columns_to_skip = c('label','metaGroupId','groupId','goodPeakCount','medMz', 'medRt','maxQuality', 'compoundId', 'expectedRtDiff','ppmDiff','parent'),
  report_pool_size_before_df=TRUE
  
)

##### Data analysis asthma - Bioproject PRJNA446042  ######
library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)
library(picante)
library(TSA)
library(nortest)
library(multcomp)
library(car)
source('http://bioconductor.org/biocLite.R')
library(devtools)
library(mctoolsr)
library(microbiome)
library(mvabund)
library(MASS)
library(nlme)
library(lmerTest)
library(geepack)c
library(doBy)
library(lattice)
library(MuMIn)
library(DESeq2)

setwd("path")

#read phyloseq
readRDS("path/ps4.RDS") -> ps4

########## ANALYSIS OF ALL SAMPLES #############
otuDA<-as.data.frame((otu_table(ps4)))
phylodiversityRAREF_QA <- pd(otuDA, phy_tree(ps4), include.root=F) ### filogenetic diversity. Include root=FALSE tree generated not a root tree
diversityRAREF_QA <- estimate_richness(ps4, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #Fisher error measure.
diversityRAREF_QA1<-cbind(sample_data(ps4),diversityRAREF_QA,phylodiversityRAREF_QA)

# Beta-Diversity comparisons with adonis
uniunA <-phyloseq::distance(ps4, method="unifrac")
uniweighA<-phyloseq::distance(ps4, method="wunifrac")
braydA <-phyloseq::distance(ps4, method="bray")
jaccdA <-phyloseq::distance(ps4, method="jaccard")

# Adonis to compare the matrix between oral and nasal samples, Aplicar strata
t1A<-adonis(uniunA~diversityRAREF_QA1$host_tissue_sampled,perm=10000, parallel = 4) ; t1A
t1A<-adonis(uniweighA~diversityRAREF_QA1$host_tissue_sampled,perm=10000, parallel = 4) ; t1A
t1A<-adonis(braydA~diversityRAREF_QA1$host_tissue_sampled,perm=10000, parallel = 4); t1A
t1A<-adonis(jaccdA~diversityRAREF_QA1$host_tissue_sampled,perm=10000, parallel = 4) ; t1A

# PCoA plots
colors <- c("nasal mucosa" = "dark green", "oral mucosa" = "orange")
p1A = plot_ordination(ps4, ordinate(ps4, method="PCoA", dist="unifrac", weighted=TRUE), type = "samples", color = "host_tissue_sampled")
p1A + geom_point(size = 2, alpha = .8) +
  theme(title = element_text(face="bold", colour="#225ea8", size=12)) +
  theme(axis.title = element_text(face="bold", colour="#225ea8", size=14)) + 
  theme(axis.text.x  = element_text(size=8, angle = 45, hjust = 1, face = "bold")) +
  theme(axis.text.y  = element_text(size=8)) +
  guides(color = guide_legend(title = "Host tissue sampled",title.position = "top")) + 
  labs(title = "PCoA UNIFRAC", caption = "Distance measure: unifrac") +
  scale_color_manual(values = colors) 

p1A = plot_ordination(ps4, ordinate(ps4, method="PCoA", dist="wunifrac"), type = "samples", color = "host_tissue_sampled") 
p1A + geom_point(size = 2, alpha = .8) +
  theme(title = element_text(face="bold", colour="#225ea8", size=12)) +
  theme(axis.title = element_text(face="bold", colour="#225ea8", size=14)) + 
  theme(axis.text.x  = element_text(size=8, angle = 45, hjust = 1, face = "bold")) +
  theme(axis.text.y  = element_text(size=8)) +
  guides(color = guide_legend(title = "Host tissue sampled",title.position = "top")) + 
  labs(title = "PCoA Unweigthed WUNIFRAC", caption = "Distance measure: wunifrac") +
  scale_color_manual(values = colors) 


########## ANALYSIS OF ALL NASAL SAMPLES #############

subset_samples(ps4, host_tissue_sampled%in%c("nasal mucosa"))-> ps4.nasal
sample_data(ps4.nasal)$host_disease <- relevel(sample_data(ps4.nasal)$host_disease, "Healthy", ref = "Healthy")
otuD<-as.data.frame((otu_table(ps4.nasal)))
phylodiversityRAREF_Q <- pd(otuD, phy_tree(ps4.nasal), include.root=F) ### filogenetic diversity. Include root=FALSE tree generated not a root tree
diversityRAREF_Q <- estimate_richness(ps4.nasal, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #Fisher error measure.
diversityRAREF_Q1<-cbind(sample_data(ps4.nasal),diversityRAREF_Q,phylodiversityRAREF_Q)

par(mfrow = c(2,2))

######## funcion multiplot ######## 

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##### graphs boxplots with jitter the outlaiers, if there are any, they are in black, the white circles are the values of your samples
#diversityRAREF_Q <- estimate_richness(phyloseq, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #Fisher error measure.
observed <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), Observed))
observed2 <-observed + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("ASV richness")+labs(y = "ASV richness")
chao <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), Chao1))
chao2<-chao + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Chao1 richness")+labs(y = "Chao1 richness")
shan <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), Shannon))
shan2<-shan + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Shannon diversity")+labs(y = "Shannon diversity")
invsimp <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), InvSimpson))
invsimp2<- invsimp + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("InvSimpson diversity")+labs(y = "InvSimpson diversity")
ACE <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), ACE))
ACE2<- ACE + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("ACE diversity")+labs(y = "ACE diversity")
phyl <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), PD))
phyl2<-phyl + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Phylogenetic diversity")+labs(y = "Faith's diversity")


#multiplot(observed2, chao2, shan2, fish2, ACE2, phyl2, cols=2)  #### to run this fucntion I have to run function multiplot first
multiplot(observed2, chao2, invsimp2, shan2, ACE2, phyl2, cols=3)  #### to run this fucntion I have to run function multiplot first


#Apply test LMER in data
##### LMER TEST ####
#test <-lmer(Observed~host_disease+collection_date+host_disease+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
#summary(test)`
#Observed
t1<-lmer(Observed~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#PD
t1<-lmer(PD~host_disease+collection_date+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#Shannon
t1<-lmer(Shannon~host_disease+collection_date+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#Chao1
t1<-lmer(Chao1~host_disease+collection_date+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#ACE
t1<-lmer(ACE~host_disease+collection_date+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#Simpson
t1<-lmer(Simpson~host_disease+collection_date+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")

#Plot glm (https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_lm.html)
#library(ggfortify)
#autoplot(t1, which = 1:6, label.size = 3)

# Beta-Diversity comparisons with adonis
uniun <-phyloseq::distance(ps4.nasal, method="unifrac")

uniweigh<-phyloseq::distance(ps4.nasal, method="wunifrac")

brayd <-phyloseq::distance(ps4.nasal, method="bray")

jaccd <-phyloseq::distance(ps4.nasal, method="jaccard")

# Referee suggests for adonis (no strata)
t1<-adonis(uniun~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4) ; t1
t1<-adonis(uniweigh~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4) ; t1
t1<-adonis(jaccd~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4) ; t1
t1<-adonis(brayd~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4); t1


colors_nasal <- c("Asthma" = "Red", "Healthy" = "Blue")
# PCoA plots
p1 = plot_ordination(ps4.nasal, ordinate(ps4.nasal, method="PCoA", dist="unifrac", weighted=TRUE), type = "samples", color = "host_disease")
p1 + geom_point(size = 2, alpha = .8) +
  theme(title = element_text(face="bold", colour="#225ea8", size=12)) +
  theme(axis.title = element_text(face="bold", colour="#225ea8", size=14)) + 
  theme(axis.text.x  = element_text(size=8, angle = 45, hjust = 1, face = "bold")) +
  theme(axis.text.y  = element_text(size=8)) +
  guides(color = guide_legend(title = "Host tissue sampled",title.position = "top")) + 
  labs(title = "PCoA Unweigthed UNIFRAC", caption = "Distance measure: unifrac") +
  scale_color_manual(values = colors_nasal)
p1 = plot_ordination(ps4.nasal, ordinate(ps4.nasal, method="PCoA", dist="wunifrac"), type = "samples", color = "host_disease") 
p1 + geom_point(size = 2) + ggtitle("PCoA weighted-UniFrac") +
  theme(title = element_text(face="bold", colour="#225ea8", size=12)) +
  theme(axis.title = element_text(face="bold", colour="#225ea8", size=14)) + 
  theme(axis.text.x  = element_text(size=8, angle = 45, hjust = 1, face = "bold")) +
  theme(axis.text.y  = element_text(size=8)) +
  guides(color = guide_legend(title = "Host disease",title.position = "top")) + 
  labs(title = "PCoA weighted-UniFrac nasal samples", caption = "Distance measure: weighted-UniFrac") +
  scale_color_manual(values = colors_nasal)
p1 = plot_ordination(ps4.nasal, ordinate(ps4.nasal, method="PCoA", dist="brayd"), type = "samples", color = "host_disease") 
p1 + geom_point(size = 2) + ggtitle("PCoA Bray-Curtis") 
p1 = plot_ordination(ps4.nasal, ordinate(ps4.nasal, method="PCoA", dist="jaccd"), type = "samples", color = "host_disease") 
p1 + geom_point(size = 2) + ggtitle("PCoA Jaccard")

######### Estimate proportions for taxonomic ranks

# Phyla

rank_names(ps4.nasal, errorIfNULL=TRUE) # view taxonomic ranks

taxonP <- tax_glom(ps4.nasal, taxrank = "Phylum")
taxon.tr.P <- transform_sample_counts(taxonP, function (x) x/sum(x))
taxon.tr.Phyla <- taxon.tr.P 
taxon.tr.Phyla

#Export table of phyla with taxonomy

table_all_phyla<-cbind(sample_data(ps4.nasal),otu_table(taxon.tr.Phyla))
tax_table(taxon.tr.Phyla)
table_all_phyla <- rename(table_all_phyla,c("ASV144" = "Fusobacteria", "ASV283" = "Firmicutes","ASV1002" = "Actinobacteria","ASV1851" = "Epsilonbacteraeota",
                                "ASV2213" = "Bacteroidetes","ASV2771" = "Proteobacteria"))
table_all_phyla
write.csv(table_all_phyla,file="table_all_phyla_nasal.csv")
# Genera
rank_names(ps4.nasal, errorIfNULL=TRUE) # view taxonomic ranks

taxonG <- tax_glom(ps4.nasal, taxrank = "Genus")
taxon.tr.G <- transform_sample_counts(taxonG, function (x) x/sum(x))
taxon.tr.Genus <- taxon.tr.G
taxon.tr.Genus

#Export table of genera with taxonomy

table_all_Genus<-cbind(sample_data(ps4.nasal),otu_table(taxon.tr.Genus))
tax_table(taxon.tr.Genus)
table_all_Genus <- rename(table_all_Genus,c(c("ASV106"="Leptotrichia",
                                              "ASV114"="Streptobacillus",
                                              "ASV144"="Fusobacterium",
                                              "ASV203"="Streptococcus",
                                              "ASV261"="Abiotrophia",
                                              "ASV283"="Dolosigranulum",
                                              "ASV296"="Granulicatella",
                                              "ASV335"="Gemella",
                                              "ASV352"="Staphylococcus",
                                              "ASV433"="Megasphaera",
                                              "ASV450"="Veillonella",
                                              "ASV494"="Ruminococcaceae_UCG-014",
                                              "ASV531"="Stomatobaculum",
                                              "ASV536"="Lachnoanaerobaculum",
                                              "ASV588"="Oribacterium",
                                              "ASV656"="Finegoldia",
                                              "ASV840"="Kocuria",
                                              "ASV848"="Pseudarthrobacter",
                                              "ASV878"="Rothia",
                                              "ASV882"="F0332",
                                              "ASV952"="Actinomyces",
                                              "ASV991"="Lawsonella",
                                              "ASV1002"="Corynebacterium_1",
                                              "ASV1015"="Corynebacterium",
                                              "ASV1596"="Paracoccus",
                                              "ASV1851"="Campylobacter",
                                              "ASV2007"="Porphyromonas",
                                              "ASV2020"="Tannerella",
                                              "ASV2059"="Prevotella_2",
                                              "ASV2081"="Prevotella_6",
                                              "ASV2125"="Prevotella_7",
                                              "ASV2179"="Prevotella",
                                              "ASV2213"="Alloprevotella",
                                              "ASV2285"="F0058",
                                              "ASV2372"="Capnocytophaga",
                                              "ASV2478"="Bergeyella",
                                              "ASV2715"="Acinetobacter",
                                              "ASV2730"="Cardiobacterium",
                                              "ASV2771"="Moraxella",
                                              "ASV2812"="Aggregatibacter",
                                              "ASV2828"="Actinobacillus",
                                              "ASV2853"="Haemophilus",
                                              "ASV2868"="Escherichia/Shigella",
                                              "ASV2879"="Salmonella",
                                              "ASV2991"="Lautropia",
                                              "ASV3042"="Neisseria",
                                              "ASV3054"="Kingella",
                                              "ASV3078"="Massilia")))

write.csv(table_all_Genus,file="table_all_Genus_nasal.csv")
# Phyla
table_all_phyla

#Proteobacteria
t1<-lmer(Proteobacteria~host_disease + collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Firmicutes
t1<-lmer(Firmicutes~host_disease + collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Actinobacteria
t1<-lmer(Actinobacteria~host_disease + collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Bacteroidetes
t1<-lmer(Bacteroidetes~host_disease + collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE) 
summary(t1)
anova(t1, test="F")

### GENERA
table_all_Genus

#Moraxella
t1<-lmer(Moraxella~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Haemophilus
t1<-lmer(Haemophilus~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Streptococcus
t1<-lmer(Streptococcus~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Dolosigranulum
t1<-lmer(Dolosigranulum~host_disease + collection_date  + (1|host_subject_id), table_all_Genus, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Corynebacterium_1
t1<-lmer(Corynebacterium_1~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE) 
summary(t1)
anova(t1, test="F")

#Staphylococcus
t1<-lmer(Staphylococcus~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE) 
summary(t1)
anova(t1, test="F")

######################## I. DATA NORMALIZATION ######################## 
library(DESeq2)
# Using Negative Binomial distribution

# IMP: if samples are going to be rarefaccionated skip this section and move to Rarefacction of samples to the minimum sample size
#Relevel so that fold-change is compared to the Healthy state
sample_data(ps4.nasal)$host_disease <- relevel(sample_data(ps4.nasal)$host_disease, "Asthma")

diagdds = phyloseq_to_deseq2(ps4.nasal, ~collection_date + host_disease) # Any variable of the metadata would work. You need one to create the DESeq object
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate size factors
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res = res[order(res$padj, na.last=NA), ]
#View(res)
#alpha = 0.05
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps4.nasal)[rownames(sigtab), ], "matrix"))
write.csv(sigtab,file="DESeqTaxa_Signif_nasal.csv")

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


########## ANALYSIS OF ALL ORAL SAMPLES #############
subset_samples(ps4, host_tissue_sampled%in%c("oral mucosa"))-> ps4.oral
sample_data(ps4.oral)$host_disease <- relevel(sample_data(ps4.oral)$host_disease, "Asthma")
otuD<-as.data.frame((otu_table(ps4.oral)))
phylodiversityRAREF_Q <- pd(otuD, phy_tree(ps4.oral), include.root=F) ### filogenetic diversity. Include root=FALSE tree generated not a root tree
diversityRAREF_Q <- estimate_richness(ps4.oral, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #Fisher error measure.
diversityRAREF_Q1<-cbind(sample_data(ps4.oral),diversityRAREF_Q,phylodiversityRAREF_Q)
#write.table(diversityRAREF_Q1, "diversityRAREF_Q1_oral.tsv", sep = "\t")

par(mfrow = c(2,2))

##### graphs boxplots with jitter the outlaiers, if there are any, they are in black, the white circles are the values of your samples
#diversityRAREF_Q <- estimate_richness(phyloseq, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #Fisher error measure.
observed <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), Observed))
observed2 <-observed + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("ASV richness")+labs(y = "ASV richness")
chao <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), Chao1))
chao2<-chao + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Chao1 richness")+labs(y = "Chao1 richness")
shan <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), Shannon))
shan2<-shan + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Shannon diversity")+labs(y = "Shannon diversity")
invsimp <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), InvSimpson))
invsimp2<- invsimp + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("InvSimpson diversity")+labs(y = "InvSimpson diversity")
ACE <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), ACE))
ACE2<- ACE + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("ACE diversity")+labs(y = "ACE diversity")
phyl <- ggplot(diversityRAREF_Q1, aes(factor(host_disease), PD))
phyl2<-phyl + geom_boxplot(aes(fill = factor(host_disease)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Phylogenetic diversity")+labs(y = "Faith's diversity")


#multiplot(observed2, chao2, shan2, fish2, ACE2, phyl2, cols=2)  #### to run this fucntion I have to run function multiplot first
multiplot(observed2, chao2, invsimp2, shan2, ACE2, phyl2, cols=3)  #### to run this fucntion I have to run function multiplot first


#Apply test LMER in data
##### LMER TEST ####
#test <-lmer(Observed~host_disease+collection_date+host_disease+T_or_C+host_sex+Host_Age+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
#summary(test)
#Observed
t1<-lmer(Observed~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#PD
t1<-lmer(PD~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#Shannon
t1<-lmer(Shannon~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#Chao1
t1<-lmer(Chao1~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#ACE
t1<-lmer(ACE~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")
#Simpson
t1<-lmer(Simpson~host_disease+collection_date+T_or_C+host_sex+drug_usage+(1|host_subject_id), diversityRAREF_Q1)
summary(t1)
anova(t1, test="F")

#Plot glm (https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_lm.html)
#library(ggfortify)
#autoplot(t1, which = 1:6, label.size = 3)

# Beta-Diversity comparisons with adonis
uniun <-phyloseq::distance(ps4.oral, method="unifrac")

uniweigh<-phyloseq::distance(ps4.oral, method="wunifrac")

brayd <-phyloseq::distance(ps4.oral, method="bray")

jaccd <-phyloseq::distance(ps4.oral, method="jaccard")

# Referee suggests for adonis (no strata)
t1<-adonis(uniweigh~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4) ; t1
t1<-adonis(brayd~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4); t1
t1<-adonis(jaccd~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4) ; t1
t1<-adonis(uniun~diversityRAREF_Q1$host_disease+diversityRAREF_Q1$collection_date,perm=10000,parallel = 4) ; t1

# PCoA plots

p1 = plot_ordination(ps4.oral, ordinate(ps4.oral, method="PCoA", dist="unifrac", weighted=TRUE), type = "samples", color = "host_disease")
p1 + geom_point(size = 5) + ggtitle("PCoA Weigthed UNIFRAC") 
p1 = plot_ordination(ps4.oral, ordinate(ps4.oral, method="PCoA", dist="wunifrac"), type = "samples", color = "host_disease") 
p1 + geom_point(size = 5) + ggtitle("PCoA Unweigthed UNIFRAC") 
p1 = plot_ordination(ps4.oral, ordinate(ps4.oral, method="PCoA", dist="brayd"), type = "samples", color = "host_disease") 
p1 + geom_point(size = 5) + ggtitle("PCoA Bray-Curtis") 
p1 = plot_ordination(ps4.oral, ordinate(ps4.oral, method="PCoA", dist="jaccd"), type = "samples", color = "host_disease") 
p1 + geom_point(size = 5) + ggtitle("PCoA Jaccard")

######### Estimate proportions for taxonomic ranks

# Phyla
rank_names(ps4.oral, errorIfNULL=TRUE) # view taxonomic ranks

taxonP <- tax_glom(ps4.oral, taxrank = "Phylum")
taxon.tr.P <- transform_sample_counts(taxonP, function (x) x/sum(x))
taxon.tr.Phyla <- taxon.tr.P 
taxon.tr.Phyla 
tax_table(taxon.tr.Phyla)
#Export table of phyla with taxonomy

table_all_phyla <-cbind(sample_data(ps4.oral),otu_table(taxon.tr.Phyla))
table_all_phyla <- rename(table_all_phyla,c("ASV144" = "Fusobacteria", "ASV194" = "Firmicutes",
                                            "ASV952" = "Actinobacteria","ASV1851" = "Epsilonbacteraeota",
                                            "ASV2213" = "Bacteroidetes","ASV2839" = "Proteobacteria"))
table_all_phyla
write.csv(table_all_phyla,file="table_all_phyla_oral.csv")

# Genera
rank_names(ps4.oral, errorIfNULL=TRUE) # view taxonomic ranks

taxonG <- tax_glom(ps4.oral, taxrank = "Genus")
taxon.tr.G <- transform_sample_counts(taxonG, function (x) x/sum(x))
taxon.tr.Genus <- taxon.tr.G
taxon.tr.Genus
tax_table(taxon.tr.Genus)
#Export table of genera with taxonomy

table_all_Genus<-cbind(sample_data(ps4.oral),otu_table(taxon.tr.Genus))
tax_table(taxon.tr.Genus)
table_all_Genus <- rename(table_all_Genus,c(c("ASV49" = "Leptotrichia",
                                              "ASV114" = "Streptobacillus",
                                              "ASV144" = "Fusobacterium",
                                              "ASV194" = "Streptococcus",
                                              "ASV261" = "Abiotrophia",
                                              "ASV283" = "Dolosigranulum",
                                              "ASV296" = "Granulicatella",
                                              "ASV335" = "Gemella",
                                              "ASV352" = "Staphylococcus",
                                              "ASV433" = "Megasphaera",
                                              "ASV441" = "Veillonella",
                                              "ASV494" = "Ruminococcaceae_UCG-014",
                                              "ASV531" = "Stomatobaculum",
                                              "ASV536" = "Lachnoanaerobaculum",
                                              "ASV588" = "Oribacterium",
                                              "ASV656" = "Finegoldia",
                                              "ASV840" = "Kocuria",
                                              "ASV848" = "Pseudarthrobacter",
                                              "ASV871" = "Rothia",
                                              "ASV882" = "F0332",
                                              "ASV952" = "Actinomyces",
                                              "ASV991" = "Lawsonella",
                                              "ASV1002" = "Corynebacterium_1",
                                              "ASV1033" = "Corynebacterium",
                                              "ASV1596" = "Paracoccus",
                                              "ASV1851" = "Campylobacter",
                                              "ASV2007" = "Porphyromonas",
                                              "ASV2020" = "Tannerella",
                                              "ASV2059" = "Prevotella_2",
                                              "ASV2081" = "Prevotella_6",
                                              "ASV2125" = "Prevotella_7",
                                              "ASV2179" = "Prevotella",
                                              "ASV2213" = "Alloprevotella",
                                              "ASV2285" = "F0058",
                                              "ASV2372" = "Capnocytophaga",
                                              "ASV2478" = "Bergeyella",
                                              "ASV2702" = "Acinetobacter",
                                              "ASV2730" = "Cardiobacterium",
                                              "ASV2771" = "Moraxella",
                                              "ASV2812" = "Aggregatibacter",
                                              "ASV2828" = "Actinobacillus",
                                              "ASV2839" = "Haemophilus",
                                              "ASV2868" = "Escherichia/Shigella",
                                              "ASV2879" = "Salmonella",
                                              "ASV2991" = "Lautropia",
                                              "ASV3047" = "Neisseria",
                                              "ASV3054" = "Kingella",
                                              "ASV3078" = "Massilia")))

write.csv(table_all_Genus,file="table_all_Genus_oral.csv")

# Phyla
table_all_phyla

#Proteobacteria
t1<-lmer(Proteobacteria~host_disease+ collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Firmicutes
t1<-lmer(Firmicutes~host_disease+ collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Actinobacteria
t1<-lmer(Actinobacteria~host_disease+ collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Bacteroidetes
t1<-lmer(Bacteroidetes~host_disease+ collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Fusobacteria
t1<-lmer(Fusobacteria~host_disease+ collection_date + (1|host_subject_id), table_all_phyla, REML = FALSE)
summary(t1)
anova(t1, test="F")

# Genera
table_all_Genus

#Streptococcus
t1<-lmer(Streptococcus~host_disease + collection_date+ (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Haemophilus
t1<-lmer(Haemophilus~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Gemella
t1<-lmer(Gemella~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Moraxella
t1<-lmer(Moraxella~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Neisseria
t1<-lmer(Neisseria~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Veillonella
t1<-lmer(Veillonella~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

#Porphyromonas
t1<-lmer(Porphyromonas~host_disease + collection_date + (1|host_subject_id), table_all_Genus, REML = FALSE)
summary(t1)
anova(t1, test="F")

######################## I. DATA NORMALIZATION ######################## 
library(DESeq2)
# Using Negative Binomial distribution

# IMP: if samples are going to be rarefaccionated skip this section and move to Rarefacction of samples to the minimum sample size
#Relevel so that fold-change is compared to the Healthy state
sample_data(ps4.oral)$host_disease <- relevel(sample_data(ps4.oral)$host_disease, "Asthma")

diagdds = phyloseq_to_deseq2(ps4.oral, ~collection_date + host_disease) # Any variable of the metadata would work. You need one to create the DESeq object
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate size factors
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res = res[order(res$padj, na.last=NA), ]
#View(res)
#alpha = 0.05
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps4.oral)[rownames(sigtab), ], "matrix"))
write.csv(sigtab,file="DESeqTaxa_Signif_oral.csv")

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


###ampvis2 all samples

###### AMPVIS2 ####
library(phytools)
library(ampvis2)
ps4@otu_table@.Data <- t(ps4@otu_table@.Data)
## extraer las tablas
otutable <- data.frame(OTU = rownames(phyloseq::otu_table(ps4)@.Data),
                       phyloseq::otu_table(ps4)@.Data,
                       phyloseq::tax_table(ps4)@.Data,
                       check.names = FALSE)
### 'otutable' contains the ASVs, the accounts for each ASV, and the corresponding lineage
## extract the metadata
metadata <- data.frame(phyloseq::sample_data(ps4), 
                       check.names = FALSE)
row.names(metadata) -> patients
col -> colnames_metada
metadata <- data.frame(metadata,patients)
col_metadata <- as.vector(colnames(sample_data(ps4)))
metadata <- metadata[,c("patients",col_metadata)]

## extract the tree and root it

rooted_tree <- midpoint.root(ps4@phy_tree)
## Ampvis2 requires that
### 1) the taxonomic ranges are seven and go from Kingdom to Species
### 2) the first column of the metadata is the identifier of each sample
#We transform the collection_date to a factor to be able to make shapes
metadata$collection_date <- as.factor(metadata$collection_date)
## Finalmente, generamos el objeto ampvis
av <- amp_load(otutable = otutable, metadata = metadata, tree = rooted_tree)
av
##### Microbiota ####
amp_heatmap(av,
            plot_values = T,
            color_vector = c("white", "dark red"),
            tax_show = 10,
            group_by = "host_tissue_sampled",
            tax_aggregate = "Genus",
            tax_add = "Kingdom",
            plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = -45, size=10, vjust = 1,hjust = 0),
        axis.text.y = element_text(size=10, face = "italic"), legend.position="right") +
  ggtitle("The most abundant taxa") +
  theme(plot.title = element_text(face="bold", colour="#225ea8", size=12)) -> heatmap
heatmap

av_oral <- amp_subset_samples(av, host_tissue_sampled %in% c("oral mucosa"))

##### Microbiota ####
amp_heatmap(av_oral,
            plot_values = T,
            color_vector = c("white", "dark red"),
            tax_show = 10,
            group_by = "host_disease",
            tax_aggregate = "Phylum",
            tax_add = "Kingdom",
            plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = -45, size=10, vjust = 1,hjust = 0),
        axis.text.y = element_text(size=10, face = "italic"), legend.position="right") +
  ggtitle("The most abundant taxa") +
  theme(plot.title = element_text(face="bold", colour="#225ea8", size=12))

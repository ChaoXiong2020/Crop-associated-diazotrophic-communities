
### Bioinformatic analysis for amplicon sequencing based on USEARCH (v10.0) ###
### More information about the script can be found at https://drive5.com/usearch/manual/ex_miseq.html ###

###
mkdir -p fq
mkdir -p output
cd fq

# Merge paired reads
$usearch -fastq_mergepairs *R1*.fq -relabel @ \
-fastq_trunctail 35 -fastq_minlen 120 \
-fastq_maxdiffs 8 -fastq_pctid 90 \
-fastq_minmergelen 340 -fastq_maxmergelen 380 \
-fastqout ../output/mergedraw.fq

# Strip primers
$usearch -fastx_truncate ../output/mergedraw.fq -stripleft 20 -stripright 20 -fastqout ../output/merged.fq

# Quality filter
$usearch -fastq_filter ../output/merged.fq -fastq_maxee 0.5 \
-fastaout ../output/filtered.fa -relabel Filt

# Find unique read sequences and abundances
$usearch -fastx_uniques ../output/filtered.fa -sizeout -relabel Uniq -fastaout ../output/uniques.fa

# Make 95% OTUs and filter chimeras
$usearch -uchime2_ref ../output/uniques.fa -db ../output_ref/nif_database.fasta -uchimeout ../output_ref/out.txt -strand plus -mode sensitive -chimeras ../output_ref/chimeras.fa -notmatched ../output_ref/nochimeras.fa -threads 20
$usearch -cluster_fast ../output_ref/nochimeras.fa -id 0.95 -centroids ../output_ref/otus_demo.fa -uc ../output_ref/clusters.uc -minsize 8 -relabel OTU -sort size -threads 20

# Make OTU table
$usearch -otutab ../output/merged.fq -otus ../output_ref/otus_demo.fa -otutabout ../output_ref/otutab_demo.txt -threads 20


### Qiime process QIIME v1.91 ### 
### More information about the script can be found at http://qiime.org/scripts/ ###

### CSS normalize process ###
normalize_table.py -i otu_table.biom -a CSS -o otu_table_norm.biom

### Normalize to 1380 reads / sample ###
single_rarefaction.py -i otu_table.biom -d 1380 -o otu_table_even.biom

### Alpha_diversity ###
alpha_diversity.py -i otu_table_even.biom -m chao1,shannon,simpson,simpson_e,heip_e,goods_coverage,observed_species,PD_whole_tree -t otus.tree -o alpha.txt

### Beta_diversity ###
beta_diversity.py -i ../qiime_output/otu_table_norm.biom -m weighted_unifrac -t ../output/otus.tree -o ../qiime_output/beta_diversity

### Beta_diversity_nmds ###
nmds.py -i ../qiime_output/beta_diversity -o ../qiime_output/nmds

### Summarize_Community composition ###
summarize_taxa.py -i ../qiime_output/otu_table.biom -o ../qiime_output/tax_summary_a -L 1,2,3,4,5,6,7 -a 
summarize_taxa.py -i ../qiime_output/otu_table.biom -o ../qiime_output/tax_summary_r -L 1,2,3,4,5,6,7



### R process ###
rm(list=ls()) 
#########
library(vegan)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(grid)
library(rfPermute)

####### Alpha diversity #######

#### Import data ####
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/alpha/")
alpha=read.csv(file="alpha.csv",header=T,check.names=FALSE,sep=",")
head(alpha)
#### Plot Figure  ###
p=ggplot(alpha, aes(x=niche, y=shannon, fill=niche,col=niche))+ 
  geom_boxplot(width=0.5,position=position_dodge(0.6),alpha=0.8)+geom_jitter(width = 0.2,alpha=0.8,size=1)
p=p+scale_colour_manual(values =c('deeppink',"#009E73","#56B4E9","#E69F00","#D55E00",'darkmagenta',"#999999","gray33"))
p=p+scale_fill_manual(values =c('deeppink',"#009E73","#56B4E9","#E69F00","#D55E00",'darkmagenta',"#999999","gray33"))
p=p+theme(axis.title.x = element_text(face="bold", colour="black", size=6),
            axis.title.y = element_text(face="bold", colour="black", size=6),
            axis.text.x = element_text(face="bold", colour="black", size=6),
            axis.text.y  = element_text(face="bold", colour="black", size=6))
p
#####

####### Beta diversity #######

### Import data ###
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/beta/")
beta=read.csv(file="alpha.csv",header=T,check.names=FALSE,sep=",")
head(beta)
### Plot Figure  ###
p1=ggplot(sp, aes(x=xvar, y=yvar, color=nich,shape=crop)) + geom_point(size=0.8,alpha=0.8)
p1=p1+scale_shape_manual(name="Crop",breaks = c('maizex','maizeq','wheat','barley'),values = c(17,24,16,22))
p1=p1+scale_colour_manual(values =c('deeppink',"#009E73",'darkmagenta',"#56B4E9",'darkmagenta',"#E69F00","#D55E00", "#999999","gray33"))
p1=p1 + theme_bw()+theme(panel.grid = element_blank())
p=p+theme(axis.title.x = element_text(face="bold", colour="black", size=6),
            axis.title.y = element_text(face="bold", colour="black", size=6),
            axis.text.x = element_text(face="bold", colour="black", size=6),
            axis.text.y = element_text(face="bold", colour="black", size=6))
p1


####### PERMANOVA #######

### Nested PERMANOVA ###
### Import data ###
setwd("F:/Desktop/r/chaoxiong/PERMANOVA/")
sp.dist=read.table(file="weighted_unifrac_zotu_table_norm.txt",header=T,check.names=FALSE,sep="")
map=read.table(file="map.txt",header=T,row.names=1,check.names=FALSE,sep="")
head(sp.dist)
head(map)
### Analysis ###
adonis_all=adonis(sp.dist~Crop_species/Niche/Stage+Niche/Treatment+Crop_species/Site,data =map,permutations=999,strata = map$Crop_species)

### PERMANOVA ###
### Import data ###
setwd("F:/Desktop/r/chaoxiong/PERMANOVA/")
sp.dist=read.table(file="weighted_unifrac_zotu_table_norm_leaf_epi.txt",header=T,check.names=FALSE,sep="")
map=read.table(file="map_leaf_epi.txt",header=T,row.names=1,check.names=FALSE,sep="")
head(sp.dist)
head(map)
### Analysis ###
adonis_leaf_epi=adonis(sp.dist~Site+Stage+Treatment,data =map,permutations=999)


### Random forest analyses ### 
##############randomForest
data=read.csv(file="data.csv",header=T,check.names=FALSE,sep=",")
set.seed(315)
rf = randomForest(y ~ ., data=data, importance=TRUE, proximity=TRUE, ntree = 19999)
print(rf)
round(importance(rf), 2)
varImpPlot(rf)  
im= as.data.frame(importance(rf))
im = im[order(im[,1],decreasing = T),]
head(im,10)
##############rfPermute
data=read.csv(file="data.csv",header=T,check.names=FALSE,sep=",")
set.seed(315)
en.rfP <- rfPermute(y ~ ., data =data, ntree = 19999, na.action = na.omit, nrep = 100, num.cores = 1)
print(en.rfP)
plot(rp.importance(en.rfP, scale = TRUE))
imp.score <- rp.importance(en.rfP, scale = TRUE)
head(imp.score,10)
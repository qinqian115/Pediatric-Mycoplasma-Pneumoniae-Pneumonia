library(microeco)
library(dplyr)
library(ggplot2)
library(magrittr)
library(agricolae)
library(Rmisc)
library(Hmisc)
library(tidyr)
library(ggbeeswarm)
library(ggpubr)
library(reshape2)
otu_ibd<- read.table("G:/metagenomics/UC_bracken_count.tsv",sep = "\t",header = T)

colnames(otu_ibd)<- gsub("GY.","",colnames(otu_ibd))
colnames(otu_ibd)<- gsub("_bracken_report","",colnames(otu_ibd))

otu_ibd<- otu_ibd[grep("Bacteria",otu_ibd$taxonomy),]
taxonomy_table<- otu_ibd[,c(1,31)]
taxonomy_table$OTU.ID<- paste("OTU",taxonomy_table$OTU.ID,sep = "_")
rownames(taxonomy_table)<- taxonomy_table$OTU.ID
taxonomy_table<- separate(taxonomy_table,"taxonomy"
                          ,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                          sep = ";")
taxonomy_table<- taxonomy_table[,-1]
taxonomy_table<- tidy_taxonomy(taxonomy_table)
otu_ibd$OTU.ID<- paste("OTU",otu_ibd$OTU.ID,sep = "_")
rownames(otu_ibd)<- otu_ibd$OTU.ID
otu_ibd<- otu_ibd[,-c(1,31)]
sample_info<- data.frame(ID=colnames(otu_ibd))
sample_info<- separate(sample_info,ID,c("Group","Replicate"),sep = "_",remove = F)
rownames(sample_info)<- sample_info$ID


dataset <- microtable$new(sample_table=sample_info,
                          otu_table =otu_ibd,
                          tax_table=taxonomy_table)
dataset$sample_sums()%>% range
dataset$rarefy_samples(sample.size= 3236000)
dataset$cal_abund()
dataset$cal_alphadiv()
dataset$cal_betadiv(Curtis = TRUE)
t1<- trans_abund$new(dataset= dataset,taxrank = "Phylum" ,ntaxa = 12)
t1$plot_bar(others_color='gray70',facet = "Group")
species_all<- dataset$taxa_abund$Species
cor_species<- cor(species_all)
pheatmap::pheatmap(cor_species)

alpha_diversity<- dataset$alpha_diversity
t2<- trans_alpha$new(dataset = dataset,group = "Group")
t2$cal_diff(method = "anova")
t2$plot_alpha(measure="Shannon")
alpha_diversity_plot<- alpha_diversity
alpha_diversity_plot$Group<- rownames(alpha_diversity_plot)
alpha_diversity_plot<- alpha_diversity_plot[,c(2,4,6,7)]
alpha_diversity_plot$Group<- rownames(alpha_diversity)
alpha_diversity_plot<- melt(alpha_diversity_plot, id.vars = c("Group"), variable.name = "Index")
alpha_diversity_plot<- separate(alpha_diversity_plot,"Group",c("Class",NA),sep = "_",remove = F)
ggplot(alpha_diversity_plot,aes(Class,value,color=Class))+
  geom_violin()+
  geom_signif(comparisons = list(c("Control","MPP")),
              test = wilcox.test,
              map_signif_level = F,
              color = "black",
              textsize = 4)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))+labs(x="")+facet_wrap(~Index,scales = "free")
ggplot(alpha_diversity_plot[c(1:18),],aes(Class,value,color=Class))+
  geom_boxplot()+stat_boxplot(geom = "errorbar",
                              width=0.3)+
  geom_signif(comparisons = list(c("Control","MPP")),
              test = wilcox.test,
              map_signif_level = T,
              color = "black",
              textsize = 4,y=c(9100,9100,9500))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))+labs(x="",y="Chao1")+ylim(6000,9700)
ggplot(alpha_diversity_plot[c(19:36),],aes(Class,value,color=Class))+
  geom_boxplot()+stat_boxplot(geom = "errorbar",
                              width=0.3)+
  geom_signif(comparisons = list(c("Control","MPP")),
              test = wilcox.test,
              map_signif_level = T,
              color = "black",
              textsize = 4,y_position = c(9100,9100,9400))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))+labs(x="",y="ACE")+ylim(6000,9600)

ggplot(alpha_diversity_plot[c(37:54),],aes(Class,value,color=Class))+
  geom_boxplot()+stat_boxplot(geom = "errorbar",
                              width=0.3)+
  geom_signif(comparisons = list(c("Control","MPP")),
              test = wilcox.test,
              map_signif_level = T,
              color = "black",
              textsize = 4,y_position = c(5.8,5.8,6.2))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))+labs(x="",y="Shannon")+ylim(2,6.5)
ggplot(alpha_diversity_plot[c(55:72),],aes(Class,value,color=Class))+
  geom_boxplot()+stat_boxplot(geom = "errorbar",
                              width=0.3)+
  geom_signif(comparisons = list(c("Control","MPP")),
              test = wilcox.test,
              map_signif_level = T,
              color = "black",
              textsize = 4,y_position = c(1,1,1.05))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))+labs(x="",y="Simpson")+ylim(0.65,1.1)

beta_diversity<- dataset$beta_diversity$jaccard
t3 <- trans_beta$new(dataset= dataset,group="Group",measure= "bray")

t3$cal_ordination(ordination ="PCoA")
class(t3$cal_ordination)
t3$plot_ordination(plot_color ="Group")
t3$plot_ordination(plot_color ="Group",plot_shape = "Group",plot_type= c("point","ellipse"))
t3$cal_group_distance()
t3$plot_group_distance(distance_pair_stat = TRUE)
distance<- t3$res_group_distance
ggplot(distance,aes(Group,Value,color=Group))+
  geom_boxplot()+stat_boxplot(geom = "errorbar",
                              width=0.3)+
  geom_signif(comparisons = list(c("Control","MPP")),
              test =t.test,
              map_signif_level = F,
              color = "black",
              textsize = 4,y_position = c(0.65,0.65,0.7))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))+labs(x="",y="Bray-Curtis distance")+ylim(0,0.8)

PCoA<- data.frame(t3$res_ordination$scores)
ggplot(PCoA,aes(PCo1,PCo2,color=Group))+
  geom_point(size=3)+
  labs(x="PC1", y = "PC2")+theme_classic()+
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=15),legend.position = "none")+
  stat_ellipse (type = "norm",level = 0.7,size=1.2)
ggplot(PCoA,aes(Group,PCo1,color=Group))+
  geom_boxplot()+stat_boxplot(geom = "errorbar",
                              width=0.3)+
  geom_signif(comparisons = list(c("Control","MPP")),
              test = wilcox.test,
              map_signif_level = T,
              color = "black",
              textsize = 4)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15))




t4<- trans_diff$new(dataset = dataset,method ="lefse",group = "Group",alpha = 0.05,
                    lefse_subgroup =   NULL)
t4$plot_diff_bar(use_number=1:30,width=0.6)
t4$plot_diff_cladogram(use_taxa_num=500,use_feature_num = 200,
                       clade_label_level=1)

# I use this script to take tcr distances calculated in tcr-dist/ph-bradley in python to customise and re-analyse PCA plots.

library(tidyverse)
library(plyr)
library(dplyr)
library(readr)
library(circlize)
require(useful)
require(ggplot2)
require(stringr)
library(plotly)
require(viridis)
library(viridis)
library(ggpubr)
library(igraph)
library(ggrepel)
library(IDPmisc)

#####getting ready#####

#Import the TCRdist matrix file: naive sort
#epitope NP or PA
f <- read_delim("JZLCK12_naive_NP_PA_subset_parsed_seqs_probs_mq20_clones_AB.dist",
                delim = " ",
                col_names = F)
f <- read_delim("JZLCK12_naive_NP_PA_parsed_seqs_probs_mq20_clones_AB_PA.dist",
                delim = " ",
                col_names = F)

f <- read_delim("JZLCK01_WTfilter_parsed_seqs_probs_mq20_clones_AB.dist",
                delim = " ",
                col_names = F)

f <- read_delim("JZLCK01_PAd10_fullWT_parsed_seqs_probs_mq20_clones_AB.dist",
                delim = " ",
                col_names = F)

f<-read_delim("JZLCK24_221205_parsed_seqs_probs_mq20_clones_AB.dist",
              delim = " ",
              col_names = F)

f<-read_delim("JZLCK28_all_CD4NP311_parsed_seqs_probs_mq20_clones_AB.dist",
              delim = " ",
              col_names = F)
f<-read_delim("JZLCK01_230428_PA_all_parsed_seqs_probs_mq20_clones_AB.dist",
              delim = " ",
              col_names = F)


nnseq<-read_tsv("JZLCK24_221205.tsv")
nnseq$clone_id <-str_c(nnseq$id, ".clone")
                
#Import the tsv file containing all the annotations for each clone: PA naive sort 1
notations <- read_tsv("JZLCK12_naive_NP_PA_subset_parsed_seqs_probs_mq20_clones_nbrdists.tsv")
notations <- read_tsv("JZLCK12_naive_NP_PA_parsed_seqs_probs_mq20_clones_nbrdists.tsv")

#or
notations <- read_tsv("JZLCK01_WTfilter_parsed_seqs_probs_mq20_clones_nbrdists.tsv")
notations <- read_tsv("JZLCK01_PAd10_fullWT_parsed_seqs_probs_mq20_clones_nbrdists.tsv")
notations <- read_tsv("JZLCK24_221205_parsed_seqs_probs_mq20_clones_nbrdists.tsv")
notations <- read_tsv("JZLCK28_all_CD4NP311_parsed_seqs_probs_mq20_clones_nbrdists.tsv")
notations <- read_tsv("JZLCK01_230428_PA_all_parsed_seqs_probs_mq20_clones_nbrdists.tsv")


## Now need to create the matrix from the distances and calculate the X and Y co-ordinates for MDS plot
b <- as.matrix(f)
rownames(b) <- b[,1]
b <- b[,-1]
colnames(b) <- rownames(b)
b <- as.dist(b)
mds <- as.data.frame(cmdscale(b))

mds$clone_id = row.names(mds) 

#Now that we have X and Y for MDS plot, we can link these to the annotation data containing other info such as gene, etc.
mds_combined = left_join(mds,notations,by="clone_id")
#mds_combined0 = left_join(mds,notations,by="clone_id")

#or use f for UMAP
c <- as.matrix(f)
c <- c[,-1]
d<- matrix(as.numeric(c), ncol=ncol(c))


c<-sapply(c, as.numeric)
umap_d <-umap(d)
plot(umap_d$layout)
mds_combined <-cbind(mds_combined, umap_d$layout)


#For JZLCK12
mds_combined = filter(mds_combined0, grepl("NP", mds_combined0$epitope))
mds_combined = left_join(mds_combined,nnseq,by="clone_id")

mds_combined$epitope = str_sub(mds_combined$epitope,1,2)

## Optional, get MFI too from Merge_IND_and_seq.R
index_updated_IND = mutate(index_updated_IND, clone_id=str_c(id, ".clone", sep="", collapse=NULL))
mds_combined = left_join(mds_combined, index_updated_IND, by="clone_id")

## Be carefule with duplicate columns

# extract genotype, need to change case by case (for JZLCK01)
mds_combined$genotype <- substr(mds_combined$subject,1,3)
mds_combined$genotype <- sub("mou", "WT", mds_combined$genotype)
mds_combined$genotype <- sub("lck", "LCKfree", mds_combined$genotype)
mds_combined$genotype <- factor(mds_combined$genotype, c("WT","LCKfree"))

# for JZLCK12

mds_combined <- mutate (
  mds_combined, genotype= str_extract(subject, "[A-Z]*")
)
mds_combined$genotype <- sub("LCK", "LCKfree", mds_combined$genotype)
mds_combined$genotype <- factor(mds_combined$genotype, c("WT","LCKfree"))

#for JZLCK28
mds_combined$genotype <- sub(".*-", "", mds_combined$epitope)
mds_combined$state <- sub("-.*", "", mds_combined$epitope)
#mds_combined<-filter(mds_combined, state!="iCD8"&state!="nCD8")

  ## Use this code for shortening the Va gene alleles to gene only

mds_combined$va_gene_short <- as.factor(str_sub(mds_combined$va_gene, 1, -4))
mds_combined$va_gene_short <- as.factor(str_replace(mds_combined$va_gene_short, "TR", ""))
mds_combined$va_gene_short <- as.factor(str_replace(mds_combined$va_gene_short, "-.*", ""))

  ## Use this code for shortening the Vb gene alleles to gene only
mds_combined$vb_gene_short <- as.factor(str_sub(mds_combined$vb_gene, 1, -4))
mds_combined$vb_gene_short <- as.factor(str_replace(mds_combined$vb_gene_short, "TR", ""))


mds_combined$va_vb <- str_c(mds_combined$va_gene_short, "_",mds_combined$vb_gene_short)

#for LCK01 rerun
mds_combined$epitope <- factor(mds_combined$epitope, levels = c("iCD8-WT", "iCD8-LCK", "nCD8-WT", "nCD8-LCK"))


####Additional stats####
## Added function to categorize TCRs based on hydrophobicity refer to Standinsky 2016
## First import p6p7hydrochart.csv as reference matrix(note: p7 is row, p6 is column)
hydrochart <- read_csv("p6p7hydrochart.csv")
hydrochart <- as.matrix(hydrochart)
rownames(hydrochart) <- hydrochart[,1]
hydrochart <- hydrochart[,-1]

# extract p7 and p6 amino acid from cdr3 
mds_combined$cdr3a_p7 <- tolower (substr(mds_combined$cdr3a,7,7))

mds_combined$cdr3a_p6 <- tolower (substr(mds_combined$cdr3a,6,6))

mds_combined$cdr3b_p7 <- tolower (substr(mds_combined$cdr3b,7,7))

mds_combined$cdr3b_p6 <- tolower (substr(mds_combined$cdr3b,6,6))

# create index matrix and find hydrophobic or not

i <- cbind(match(mds_combined$cdr3a_p7, rownames(hydrochart)), match(mds_combined$cdr3a_p6, colnames(hydrochart)))
ib<- cbind(match(mds_combined$cdr3b_p7, rownames(hydrochart)), match(mds_combined$cdr3b_p6, colnames(hydrochart)))
mds_combined$cdr3a_hydro <- hydrochart[i]
mds_combined$cdr3b_hydro <- hydrochart[ib]
as.factor(mds_combined$cdr3a_hydro)
as.factor(mds_combined$cdr3b_hydro)


# Calculate the cdr3 aa length
mds_combined$cdr3a_length <- str_length(mds_combined$cdr3a)
mds_combined$cdr3b_length <- str_length(mds_combined$cdr3b)

##Cys index
##Find middle AA WIP, for now try N/2 as middle, later try N/2+1

f1 <- function(str1){
  N <- nchar(str1)
  if(!N%%2){
    res <- substr(str1, N/2-2, N/2+2) 
  }
  else{
    
    N1 <- median(sequence(N))
    res <- substr(str1, N1-2, N1+2)
  }
  res
}
mds_combined$cdr3amid <- 0
mds_combined$cdr3bmid <- 0

for (i in 1:nrow(mds_combined)){
  mds_combined$cdr3amid[i]<-f1(mds_combined$cdr3a[i])
}
for (i in 1:nrow(mds_combined)){
  mds_combined$cdr3bmid[i]<-f1(mds_combined$cdr3b[i])
}

#Find how many contain free cys first
table(mds_combined$genotype, str_detect(mds_combined$cdr3amid,"C"))


# get va vb pairing with this and find the most common ones
mds_combined$va_vb= str_c(mds_combined$va_gene_short, mds_combined$vb_gene_short, sep="_")
as.factor(mds_combined$va_vb)
va_vbcount = as.data.frame(table(mds_combined$va_vb))
sum(va_vbcount$Freq)
ggplot(va_vbcount[va_vbcount$Freq>10, ], aes(x=Var1, y=Freq)) + geom_col()

##try K means clustering??
mds_k <- mds[,1:2]
kmclusters <- kmeans(mds_k, 4)
k <- as.data.frame(kmclusters$cluster)
k <- mutate(k, clone_id= rownames(k))
mds_combined <- left_join(mds_combined, k, by="clone_id")

k_table <-table(clust=mds_combined$`kmclusters$cluster`, geno=mds_combined$genotype)

fviz_cluster(kmclusters, geom="point", data=mds_k)

write.csv(mds_combined, "mds_combined_NP311CD4.csv", row.names=FALSE)
write.csv(k_table, "mds_combined_PAd10_ktable.csv", row.names=FALSE)

####summary tables####
subject_pairing <- table(mouse=mds_combined_k$subject, geno=mds_combined_k$va_vb)
write.csv(subject_pairing, "mds_combined_NP311CD4_subjectpairtable.csv", row.names=TRUE)

subject_va <- table(mouse=mds_combined_k$subject, geno=mds_combined_k$va_gene_short)
write.csv(subject_va, "mds_combined_NP311CD4_subjectvatable.csv", row.names=TRUE)

subject_vb <- table(mouse=mds_combined_k$subject, geno=mds_combined_k$vb_gene_short)
write.csv(subject_vb, "mds_combined_NP311CD4_subjectvbtable.csv", row.names=TRUE)

####Save output####
write.csv(mds_combined, "mds_combined_JZLCK28-1.csv", row.names=FALSE)

####Distribution plots####
##First need to recreate the nbrdist files within repertoires
##All would use AB_wtd_nbrdist10
mds_combined$self_nbr<- ifelse(
  mds_combined$genotype == "WT",
  mds_combined$NP311_AB_wtd_nbrdist10,
  mds_combined$NP311_LCK_AB_wtd_nbrdist10
)

mds_combined$self_nbr<- ifelse(
  mds_combined$genotype == "WT",
  mds_combined$NP_AB_wtd_nbrdist10,
  mds_combined$NP_LCKfree_AB_wtd_nbrdist10
)

#for JZLCK12, needs debugging

mds_combined$self_nbr<- ifelse(
  mds_combined$genotype == "WT",
  mds_combined$NP_AB_wtd_nbrdist10,
  mds_combined$NP_LCKfree_AB_wtd_nbrdist10
)



ggplot(mds_combined, aes(x=self_nbr, color=subject))+
  geom_density(alpha = 0.2) + theme_bw()+ xlim(0,500)+
  #scale_color_manual(values=c("black", "red"))+
  facet_grid(cols= vars(epitope))
#summary table

Nbrtable <-tapply(mds_combined$self_nbr, mds_combined$subject, summary)

mds_combined$other_nbr<- ifelse(
  mds_combined$genotype == "WT",
  mds_combined$NP311_LCK_AB_wtd_nbrdist10,
  mds_combined$NP311_AB_wtd_nbrdist10
)

#For JZLCK28
mds_combined$self<-str_c(mds_combined$epitope,"_AB_wtd_nbrdist10")
df1<- select(mds_combined, clone_id, self)
df1<- bind_cols(df1, select(mds_combined, matches("_AB_wtd_nbrdist10")))
df1<- pivot_longer(df1, matches("_AB_wtd_nbrdist10"), names_to = "name", 
                   values_to = "self_nbr")
df1<-filter(df1, df1$self==df1$name)
mds_combined<- bind_cols(mds_combined, select(df1, self_nbr))


filter(mds_combined, str_detect(mds_combined$state, "CD4"))%>%
  ggplot(aes(x=self_nbr, color=factor(genotype,c("WT", "LCK")), linetype=state))+
  geom_density(alpha = 0.2) + theme_bw()+ xlim(0,500)+labs(color = "Genotype")+
  scale_color_manual(values=c("black", "red"))

mds_combined%>%
  group_by(subject, epitope)%>%
  summarise(median=median(self_nbr))%>%
  write.csv("JZLCK28 PA median nbr.csv", row.names = TRUE)

#hydrophobicity
filter(mds_combined, str_detect(mds_combined$state, "iCD4"))%>%
  ggplot(aes(x= factor(cdr3b_hydro), y=self_nbr, fill=factor(genotype,c("WT", "LCK"))))+
  geom_dotplot(binaxis = "y", stackdir = "center", position=position_dodge(0.8)) +
  theme_bw() +labs(fill = "Genotype")+
  scale_fill_manual(values=c("black", "red"))

filter(mds_combined, str_detect(mds_combined$state, "iCD4"))%>%
  ggplot(aes(genotype, fill= factor(cdr3a_hydro)))+
  geom_bar(position = "fill")+
  theme_bw() #+labs(fill = "Genotype")+
  scale_fill_manual(values=c("black", "red"))

  ggplot(mds_combined, aes(x= factor(cdr3b_hydro), y=cdr3b_length, fill=factor(genotype,c("WT", "LCKfree"))))+
    geom_dotplot(binaxis = "y", stackdir = "center", position=position_dodge(0.8)) +
    theme_bw() +labs(fill = "Genotype")+
    scale_fill_manual(values=c("black", "red"))
 
  #check correlation? 
  ggplot(mds_combined, aes(x= self_nbr, y=log(cdr3a_protseq_prob), color=factor(genotype,c("WT", "LCKfree"))))+
    geom_point() +
    theme_bw() +labs(color = "Genotype")+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    scale_color_manual(values=c("black", "red"))
  
  
  
  ggplot(mds_combined, aes(genotype, fill= factor(cdr3a_hydro)))+
    geom_bar(position = "fill")+
    theme_bw() #+labs(fill = "Genotype")+
  scale_fill_manual(values=c("black", "red"))



# generation probability
filter(mds_combined, str_detect(mds_combined$state, "CD4"))%>%
  ggplot(aes(x=-log(a_protseq_prob), color=factor(genotype,c("WT", "LCK")), linetype=state))+
  geom_density(alpha = 0.2) + theme_bw()+ #xlim(0,500)+
  labs(color = "Genotype")+
  scale_color_manual(values=c("black", "red"))+
  facet_wrap(~state, 4)




# for other stats like cdr3 length

##output LCK seqs that has >120
mds1 <- filter(mds_combined, genotype == "LCKfree", other_nbr>120)
write.csv(mds1, "filteredTCRs.csv")

#### Clone size distribution plot
mds_combined<- arrange(mds_combined, subject, clone_size)
mds_combined<-group_by(mds_combined, subject)%>%mutate(clone_rank=(1:n())/n())
mds_combined<-group_by(mds_combined, subject)%>%mutate(seq_rank=cumsum(clone_size)/sum(clone_size))
#info bomb
ggplot(data = mds_combined, aes(x=clone_rank, y=seq_rank, group=subject)) + geom_line(aes(color = genotype))+
  geom_point(aes(fill=vb_gene_short, size=clone_size), shape=21)+ theme_bw()+
  labs(colour = "Genotype")+
  scale_color_manual(name= "Genotype", values=c("WT"="lightgrey", "LCKfree"="#F8766D"))
#simple
ggplot(data = mds_combined, aes(x=clone_rank, y=seq_rank, group=subject, color=genotype)) + geom_line()+
  geom_point()+ theme_bw()+
  labs(colour = "Genotype")+
  scale_color_manual(name= "Genotype", values=c("WT"="lightgrey", "LCKfree"="#F8766D"))

mds_combined%>% group_by(subject)%>% summarise (count=n())
mds_combined<- filter(mds_combined, subject!="mouse_subject0006"&subject!="mouse_subject0029")


      
##### kPCA plots: Use 4x6 inch export as pdf for single PCA, 4x10 (or 750 x 300 .svg)for double, 5x7 for 2by2 #####
# Available modifications: + scale_x_reverse()+ ylim(-220,110). Use factor() to convert number to categories. str_length (). Use scale_color_manual(name= "", values=c(""="red")) to set color 
### Dot plots
#separated
ggplot(data = mds_combined, aes(x=V1, y=V2, colour = log(a_protseq_prob), size=clone_size)) + 
  labs(colour = "CDR3a probability")+scale_y_reverse()+
  geom_point()+ theme_bw()+scale_color_viridis()+
  facet_grid(~genotype)

ggplot(data = mds_combined, aes(x=V1, y=V2, colour = cdr3b_length, size=clone_size)) + 
  labs(colour = "CDR3b length")+scale_y_reverse()+
  geom_point()+ theme_bw()+scale_color_viridis()+
  facet_grid(~genotype)


ggplot(data = mds_combined, aes(x=V1, y=V2, colour = va_gene_short, size=clone_size)) + 
  labs()+theme_bw()+ 
  geom_point()+ 
  facet_wrap(~epitope, 2)

ggplot(data = mds_combined, aes(x=V1, y=V2, colour = factor(`kmclusters$cluster`), size=clone_size)) + 
  labs(colour = "K-means cluster")+theme_bw()+ 
  geom_point()+
  facet_wrap(~epitope, 2)

ggplot(data = mds_combined, aes(x=mds_combined$'1', y=mds_combined$'2', colour = vb_gene_short, size=clone_size)) + 
  labs()+
  geom_point()+ theme_bw()+
  facet_wrap(~epitope, 2)


ggplot(data = mds_combined, aes(x=V1, y=V2, colour = self_nbr, size=clone_size)) + 
  labs(colour = "NNdist10")+theme_bw()+ 
  geom_point()+scale_colour_viridis()+
  facet_grid(cols= vars(genotype))

ggplot(data = filter(mds_combined, epitope.x=="PA"), aes(x=V1, y=V2, colour = log10(X.488..695.40.A))) + 
  labs(colour = "Log(CD5 MFI)")+ scale_x_reverse()+ scale_y_reverse()+
  geom_point()+ theme_bw()+scale_colour_viridis()+
  facet_grid(cols= vars(genotype))


#overlay
ggplot(data = mds_combined, aes(x=V1, y=V2, colour = factor(genotype), size=clone_size)) + 
  scale_color_manual(name= "Genotype", values=c("WT"="dimgrey", "LCK"="#F8766D")) +theme_bw()+ 
  scale_shape_manual(values = c(21:24))+
  geom_point(shape =1, stroke=1)+ facet_wrap(~state, 1)

ggplot(data = mds_combined, aes(x=mds_combined$'1', y=mds_combined$'2', colour = factor(genotype), size=clone_size)) + 
  scale_color_manual(name= "Genotype", values=c("WT"="dimgrey", "LCK"="#F8766D")) +theme_bw()+ 
  scale_shape_manual(values = c(21:24))+
  geom_point(shape =1, stroke=1)+ facet_wrap(~state, 1)


filter(mds_combined, str_detect(mds_combined$state, "iCD4"))%>%ggplot(aes(x=V1, y=V2, colour = factor(genotype), size=clone_size)) + 
  labs(colour = "Genotype")+theme_bw()+ 
  scale_color_manual(name= "Genotype", values=c("WT"="dimgrey", "LCK"="#F8766D")) +
  geom_point(shape =1, stroke=1) + facet_wrap(~state, 1)

filter(mds_combined, mds_combined$state!="iCD8")%>%ggplot(aes(x=V1, y=V2, colour = state, size=clone_size)) + 
  labs(colour = "Genotype")+theme_bw()+ 
  scale_color_manual(values=c("grey", "#F8766D","blue")) +
  geom_point(shape =1, stroke=1) + facet_wrap(~genotype, 1)

##Check correlation?
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

ggplot(data = mds_combined, aes(x= self_nbr, y=log(cdr3b_protseq_prob),  colour = cdr3b_length, size=clone_size)) + 
  labs(colour = "CDR3b length")+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_text(x = 25, y = 300, label = lm_eqn(mds_combined), parse = TRUE)+
  geom_point()+ theme_bw()+scale_color_viridis()+
  facet_grid(~genotype)

##add annotations 

mds_combined$delabel <- NA
mds_combined$delabel[mds_combined$clone_size > 6] <- mds_combined$va_vb[mds_combined$clone_size > 6]

ggplot(data=mds_combined, aes(x=V1, y=V2, colour = factor(genotype), size=clone_size))+
  geom_point(shape =1, stroke=1)+ labs(colour = "Genotype")+theme_bw()+ 
  geom_text_repel(aes(label=delabel, size=6))+
  scale_color_manual(name= "Genotype", values=c("WT"="#00BFC480", "LCKfree"="#F8766D")) +
  theme(legend.position="none")

##ggplotly

p <- ggplot(data = mds_combined, aes(x=V1, y=V2, colour = self_nbr, shape = epitope,
                                       size=clone_size, mouse=subject, id= clone_id,
                                       CDR3a = cdr3a, CDR3b = cdr3b, 
                                     VA=va_gene_short, VB=vb_gene_short, 
                                     nbrdist=self_nbr, ahydro=cdr3a_hydro, bhydro=cdr3b_hydro)) + 
  labs()+theme_bw()+ scale_color_viridis()+
  geom_point() 

ggplotly(p)

## trying to make chord diagram. First create df of stats "counts". 
## Then get another df of the unique tracks and generate colours. Then subset df for each chord diagram


##optional: shorten va vb
mds_combined$va_gene_short<- str_match(mds_combined$va_gene_short, "[A-Z]*[0-9]*")
mds_combined$vb_gene_short<- str_match(mds_combined$vb_gene_short, "[A-Z]*[0-9]*")

###circos plots###
counts <- ddply(mds_combined, .(mds_combined$vb_gene,mds_combined$va_gene, mds_combined$genotype), nrow)
colnames(counts)<- c("vb", "va", "genotype","clones")

colours <- data.frame(rbind(as.matrix(table(counts$vb, counts$genotype)),as.matrix(table(counts$va, counts$genotype))))
colours$c <-rand_color(nrow(colours), transparency=0.8)
colours$c[rownames(colours)=="BV1"]="#EA782DE6"
colours$c[rownames(colours)=="BV3"]="#5386B3E6"

circos.clear()


##debug needed
counts_wt<- filter(counts,genotype=="WT")
colours_wt<- filter(colours,colours$WT!=0)

counts_lck<- filter(counts,genotype=="LCKfree")
colours_lck<- filter(colours,colours$LCKfree!=0)


chordDiagram(counts_wt, annotationTrack = "grid", grid.col=colours_wt$c, preAllocateTracks = 1)
chordDiagram(counts_lck, annotationTrack = "grid", grid.col=colours_lck$c, preAllocateTracks = 1)


circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1,sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.8)
  
}, bg.border = NA)


###archive

ggplot(data = mds_combined, aes(x=V1, y=V2, colour = log(cdr3a_protseq_prob), size=clone_size)) + 
  labs()+
  geom_point() + 
  facet_grid(cols= vars(genotype))

filter(mds_combined, str_detect(mds_combined$state, "iCD4"))%>%
  ggplot(aes(x=V1, y=V2, colour = vb_gene_short, size=clone_size)) + 
  labs()+
  geom_point() + 
  facet_grid(cols= vars(genotype))

ggplot(data = mds_combined, aes(x=vb_gene_short, y=X.488..695.40.A)) + 
  labs()+
  geom_jitter(height = 2, width = 2) +
  facet_grid(cols= vars(genotype))+scale_y_reverse()


#Plot a FACS plot, 49d PE vs 44 FITC
ggplot(data = mds_combined, 
  aes(x=X.561..780.60.A, y=X.488..530.30.A, color=vb_gene_short)) +
  scale_y_log10() +
  scale_x_log10()+
  geom_point() +
  facet_grid(cols= vars(genotype)) +
  labs(x="CD49d PE-Cy7", y="CD44 FITC")

  

#For coloring va_vb
ggplot(data = mds_combined, aes(x=V1, y=V2, colour =va_vb, size = clone_size)) + 
    geom_point() + theme_bw()+
    scale_color_manual(name= "va_vb", values=c("AV6-5_BV1"="#FF6464")) +
    facet_grid(cols= vars(genotype))

#coloring motifs that you found
ggplot(data = mds_combined, aes(x=V1, y=V2, size=clone_size)) + 
  geom_point(color="grey55") + 
  #geom_point(data=mds_combined[str_detect(mds_combined$cdr3a, "NNNN"), ], color="red")+
  #geom_point(data=mds_combined[str_detect(mds_combined$cdr3b, "DRG")&mds_combined$va_gene_short=="TRAV21", ], color="red")+
  geom_point(data=mds_combined[str_detect(mds_combined$cdr3a, "ANT"), ], color="blue")+
  facet_grid(cols= vars(genotype))

#For NP coloring va_vb
ggplot(data = mds_combined, aes(x=V1, y=V2, colour =va_vb, size = clone_size)) + 
  geom_point() + 
  scale_color_manual(name= "va_vb", values=c("TRAV16_TRBV13"="#FF6464", "TRAV7_TRBV2"="#F9E79F","TRAV16_TRBV5"="#3498DB","TRAV3_TRBV13"="#D2B4DE")) +
  facet_grid(cols= vars(genotype))

#coloring CDR3a length
ggplot(data = mds_combined, aes(x=V1, y=V2, colour = factor (cdr3a_length),size = clone_size)) + 
    geom_point()+ 
    scale_color_manual(name= "cdr3a_length", values=c("15"="red"))+
    facet_grid(cols= vars(genotype)) +scale_y_reverse()

#coloring hydrophobicity
ggplot(data = mds_combined, aes(x=V1, y=V2, colour = cdr3a_hydro)) + 
  geom_point()+ 
  scale_color_manual(name= "CDR3a_hydro", values=c("lipophilic"="red", "hydrophilic"="blue", "neutral"="grey60"))+
  facet_grid(cols= vars(genotype))#+scale_y_reverse()

#coloring protseqprob
ggplot(data = mds_combined, aes(x=V1, y=V2, colour = cdr3a_protseq_prob ,size = clone_size)) + 
  geom_point()+ 
  facet_grid(cols= vars(genotype))#+scale_y_reverse()


# Using plotly
p <- mds_combined %>% filter(genotype=="WT") %>%ggplot(aes(x=V1, y=V2, label = clone_id, CDR3a = cdr3a, CDR3b = cdr3b, va=va_gene_short, vb=vb_gene_short)) + geom_point() +scale_y_reverse()
ggplotly(p)

# WIP JZ Calculate the cumulative distribution of clone size, first create data frame s with clone size and genotype.

s <- data.frame("clone_size"=mds_combined$clone_size, "genotype"=mds_combined$genotype)


#####playing with stats#####
na.omit(mds_combined)
# Do some stats on cdr3 length
ggplot(mds_combined, aes(x= epitope, y= log(cdr3a_length))) +geom_boxplot()
anova1<-aov(cdr3a_length~ epitope, mds_combined)
plot(TukeyHSD(anova1), las=1)

# b and a protseq prob
ggplot(mds_combined, aes(x= epitope, y= log(b_protseq_prob))) +geom_boxplot()
mds_combined$logbprot <- log(mds_combined$b_protseq_prob)
mds_combined$logaprot <- log(mds_combined$a_protseq_prob)
mds_combined <-NaRV.omit(mds_combined)
mds_combined <-mds_combined[complete.cases(mds_combined),]

anova2<-aov(logbprot~ 1+epitope, mds_combined)
summary(anova2)
tuk2 <-TukeyHSD(anova2)
plot(TukeyHSD(anova2), las=1)

t.test(logbprot~ 1+epitope, filter(mds_combined, str_detect(mds_combined$state, "n")))

anova3<-aov(logaprot~ 1+epitope, mds_combined)
summary(anova3)
tuk3 <-TukeyHSD(anova3)
plot(TukeyHSD(anova3), las=1)

# nbrdist
anova4<- aov(self_nbr~ 1+epitope, mds_combined)
summary(anova4)
tuk4 <-TukeyHSD(anova4)
plot(TukeyHSD(anova4), las=1)

t.test(self_nbr~ 1+epitope, filter(mds_combined, str_detect(mds_combined$state, "n")))


# Hydro table
hydro_distro=as.data.frame (table(mds_combined$cdr3a_hydro, mds_combined$genotype))
ggplot(hydro_distro, aes(x=Var1, y=Freq)) +geom_dotplot() #+  facet_grid(cols= vars(Var2))




aprotlenfit = lm(na.omit(log(b_protseq_prob))~ na.omit(epitope), mds_combined)
cdr3alenfit1 = lm(cdr3a_length~ genotype + cdr3a_hydro, mds_combined)
summary(cdr3alenfit1)
confint(cdr3alenfit1)
anova(cdr3alenfit,cdr3alenfit1)

# exploring wt data...
# Is there correlation between cdr3a_length, cdr3a_hydro and cdr3a_protseq_prob
wt = filter(mds_combined, genotype=="wt")
lckfree = filter(mds_combined, genotype=="lckfree")

# linear models..
lenfit1= lm(cdr3a_length ~ cdr3a_hydro, wt)
confint(lenfit1)
lenfitlck1=lm(cdr3a_length ~ cdr3a_hydro, lckfree)
confint(lenfitlck1)
## so the results show in wt, hydrophilic a chains are longer, but not in lck?

lenfit2= lm(cdr3a_protseq_prob ~ cdr3a_length, data=wt)
confint(lenfit2)
plot(lenfit2)



ggplot(mds_combined, aes(x=cdr3a_hydro, y=cdr3a_length))+ 
  geom_boxplot()+
  facet_grid(cols= vars(genotype))
# I like the plot below
ggplot(mds_combined, aes(cdr3a_length))+ 
  geom_histogram()+
  facet_grid(rows= vars(cdr3a_hydro), cols= vars(genotype))

# correlations...
cor.test(wt$cdr3a_length, wt$cdr3a_hydro, method=c("pearson", "kendall", "spearman"))
plot(wt$cdr3a_length, wt$cdr3a_protseq_prob)




# get va vb pairing with this and find the most common ones
mds_combined$va_vb= str_c(mds_combined$va_gene_short, mds_combined$vb_gene_short, sep="_")
as.factor(mds_combined$va_vb)
va_vbcount = as.data.frame(table(mds_combined$va_vb))

ggplot(va_vbcount[va_vbcount$Freq>10, ], aes(x=Var1, y=Freq)) + geom_col()

## maybe get a more concise mds_combined?
mds_short = mds_combined[c("clone_id",
                           "genotype",
                           "cdr3a", 
                           "cdr3b",
                           "cdr3a_protseq_prob",
                           "cdr3b_protseq_prob",
                           "members", 
                           "clone_size",
                           "cdr3a_hydro",
                           "cdr3a_length", 
                           "va_gene_short",
                           "vb_gene_short")]
plot(mds_short)

mds_stats = mds_combined[c("genotype",
                          "cdr3a_protseq_prob",
                           "cdr3b_protseq_prob",
                           
                           "clone_size",
                           "cdr3a_hydro",
                           "cdr3a_length", 
                           "va_gene_short",
                           "vb_gene_short")]
plot(mds_stats)

#260224 plotting MFI on JZLCK24
joined <- read.csv("JZLCK24_joined.csv")
joined$genotype<- str_sub(joined$subject, start = 0, end=1)
ggplot(joined, aes(x=va_vb, y=X, color= genotype, shape=subject))+geom_boxplot()+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# #now i need to extract the "naive" or "immune" state from clone_id to identify TCRs from the naive and immune repertoires
# mds_combined$clone_id = str_replace(mds_combined$clone_id,pattern = "mouse_",replacement = "")
# corner(mds_combined)
# 
# mds_combined$state = unlist(lapply(str_split(mds_combined$clone_id,pattern = "_"),"[[",2))
# ggplot(data = mds_combined, aes(x=V1, y=V2, colour = state)) + 
#   geom_point()

# ## Use this code to filter out certain features from a column ie. naive/immune
# mds_combined$state[mds_combined$state == "immune.clone"]
# # OR
# mds_combined %>% filter(state=="naive.clone")

# #OR You can create a new object with filtered data if you prefer by using this code
# mds_combined_immune <- mds_combined %>% filter(state=="naive.clone")

## Use this code to colour based on state
# mds_combined$state = unlist(lapply(str_split(mds_combined$clone_id,pattern = "_"),"[[",2))
# ggplot(data = mds_combined_immune, aes(x=V1, y=V2, colour = state, size = clone_size)) + 
#   geom_point()
# 
# mds_combined$state = unlist(lapply(str_split(mds_combined$clone_id,pattern = "_"),"[[",2))
# ggplot(data = mds_combined, aes(x=V1, y=V2, colour = cdr3b_protseq_prob, size = clone_size)) + 
#   geom_point() + theme_linedraw() + scale_size(range = c(1,20)) 

#+ facet_wrap(vars(state), ncol = 3)
# 
# mds_combined$state = factor(mds_combined$state, levels = c("naive.clone", "immune.clone", "memory.clone"), labels = c("Naive", "Immune", "Memory"))

######### PLOTLY########
# Make the plot into an object, this code includes labeling of clone_id and subject
p <- ggplot(data = mds_combined, aes(x=V1, y=V2, colour = genotype, label = clone_id, CDR3a = cdr3a, CDR3b = cdr3b, VA=va_gene_short, VB=vb_gene_short)) +
  geom_point() + theme_linedraw() + scale_size(range = c(1,20)) + theme(text = element_text(size = 20)) + scale_color_brewer(palette="Set2")


p <- ggplot(data = mds_combined, aes(x=V1, y=V2, colour = cdr3b_protseq_prob, label = clone_id, size = clone_size, CDR3a = cdr3a, CDR3b = cdr3b)) +
  geom_point() + theme_linedraw() + scale_size(range = c(1,20)) + theme(text = element_text(size = 20)) + scale_fill_viridis_c()


p <- ggplot(data = mds_combined, aes(x=V1, y=V2, colour = state, label = clone_id, size = clone_size, CDR3a = cdr3a, CDR3b = cdr3b)) +
  geom_point() + theme_linedraw() + scale_size(range = c(1,20)) + scale_colour_viridis_d()
p

p <- ggplot(data = mds_combined, aes(x=V1, y=V2, colour = cdr3b_protseq_prob, label = clone_id, size = clone_size, CDR3a = cdr3a, CDR3b = cdr3b)) +
  geom_point() + theme_linedraw() + scale_size(range = c(1,10)) + theme(text = element_text(size = 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~state, ncol = 3)


p + scale_fill_gradient(low="blue", high="red")


p + scale_fill_viridis_c()



#I want hue
mds_combined$state = factor(mds_combined$state, levels = c("naive.clone", "immune.clone", "memory.clone"), labels = c("Naive", "Immune", "Memory"))
mds_combined$state = factor(mds_combined$state, levels = c("naive.clone", "immune.clone", "memory.clone"))

p + scale_colour_viridis_d()

 + scale_fill_gradient(low="blue", high="red")
  
# Use this for multiple plots
 + facet_wrap(vars(vb_gene_short), ncol = 3) 

## Use plotly on the plot object created above


ggplotly(p)

ggplot(p)


## igraph
tsne(f,colvec=c('gold'))
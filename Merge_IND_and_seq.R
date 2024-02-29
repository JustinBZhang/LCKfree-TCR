library(tidyverse)
library(dplyr)
library(readr)
library(plotly)
library(htmlwidgets)
library(CGPfunctions)
### using the TCR explore web tool, get and merge the IND files first. 
### to the IND csv file, add a column called "plate" and fill in p1, p2 etc. 

IND <- bind_rows(read.csv("JZLCK11_INDWT1-3.csv"), read.csv("JZLCK11_INDLCK1-3.csv"), read.csv("JZLCK12_IND.csv"))

IND <- read.csv("JZLCK21_IND.csv")
IND <- read.csv("JZLCK24 IND.csv")
IND <- read.csv("JZLCK28 IND.csv")


Loc_to_ID <- read_csv("Loc_to_ID.csv")
index_updated_IND <- merge(IND,Loc_to_ID,by=c("XLoc","YLoc"))

### create id as plate_well
index_updated_IND$id <- 
  str_c(index_updated_IND$plate, index_updated_IND$well, sep="_", collapse = NULL)

##for LCK24, create id as
index_updated_IND$id<-str_sub(index_updated_IND$name, 24, 27)
index_updated_IND$id<-str_trim(index_updated_IND$id)
index_updated_IND$id<-str_c(index_updated_IND$id, index_updated_IND$well, sep="_")
##for LCK28, create id as
index_updated_IND$id<-str_sub(index_updated_IND$name, 18, 19)
index_updated_IND$id<-str_trim(index_updated_IND$id)
index_updated_IND$id<-str_c(index_updated_IND$id, index_updated_IND$well, sep="_")

Seqs<- read_tsv("JZLCK24_221205_parsed_seqs_probs.tsv")
Seqs<- read_tsv("JZLCK28_naiveCD4_NP311_parsed_seqs_probs_mq20_clones_nbrdists.tsv")

joined <- left_join(Seqs, index_updated_IND, by= "id")
write.csv(joined, "JZLCK28_joined.csv", row.names=FALSE)
#shortcut
joined <- read.csv("JZLCK24_joined.csv")


## Use this code for shortening the Va gene alleles to gene only

joined$va_gene_short <- as.factor(str_sub(joined$va_gene, 1, -4))
joined$va_gene_short <- as.factor(str_replace(joined$va_gene_short, "TR", ""))

## Use this code for shortening the Vb gene alleles to gene only
joined$vb_gene_short <- as.factor(str_sub(joined$vb_gene, 1, -4))
joined$vb_gene_short <- as.factor(str_replace(joined$vb_gene_short, "TR", ""))


joined$va_vb <- str_c(joined$va_gene_short, "_",joined$vb_gene_short)
##detect V gene matches and compare stats
joined$match<-joined$va_vb=="AV6-7/DV9_BV3"
joined$ratio<-joined$X.561..582.15.A/joined$X.640..670.14.A
joined$genotype<- str_sub(joined$subject, start = 0, end=1)

ggplot(joined, aes(x=va_vb, y=ratio, color= genotype))+
  geom_point(shape =1, stroke=0.5)+
  stat_summary( geom = "point",
                fun = "mean",
                col = "black",
                shape="-",
                size=6
               )+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name= "Genotype", values=c("W"="dimgrey", "L"="#F8766D"))

ggplot(joined, aes(x=va_vb, y=ratio))+geom_boxplot()+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#filter(joined, epitope=="WT")%>%
ggplot(joined, aes(x=match, y=ratio))+
  geom_boxplot()+ geom_point()+
  theme_bw()+
  facet_grid(~subject)

Plot2WayANOVA(ratio ~ subject * match, joined, plottype = "line")

  filter(joined, epitope=="WT")%>%
    ggplot(aes(x=X.488..695.40.A, y=ratio, color=match))+
    geom_point()+ 
    theme_bw()#+
  facet_grid(~subject)
  
  
joined%>%
  group_by(subject, match)%>%
  summarise(median=median(ratio))

### determine NP or PA (note: NP-PE, PA-APC)
index_updated_IND <- 
  mutate(index_updated_IND, epitope = 
           ifelse(index_updated_IND$X.561..582.15.A>1000, "NP",
                  ifelse(index_updated_IND$X.640..670.14.A > 1500, "PA", "none")))

id_epitope <- select(index_updated_IND, c("id", "epitope"))
write.table(id_epitope, "id_epitope.csv", row.names=FALSE,sep="\t")


###now join the paired seq with epitope

fullseq= read.csv("JZLCK12_fullseq.csv")

fullseq= read.delim("JZLCK21_NP311CD8.tsv")

TCRdist <- left_join(fullseq, id_epitope, by= "id")
write.table(TCRdist, "joined.tsv", row.names=FALSE,sep="\t")
##join with the full sequence
TCRdist1 <- left_join(fullseq, index_updated_IND, by= "id")
write.table(TCRdist1, "JZLCK21_joined_IND.tsv", row.names=FALSE,sep="\t")

##or join with the TCRdist parsed seq

parsed <- read.delim("JZLCK12_Naive_NP_PA_subset_parsed_seqs_probs.tsv")
parsed_IND <-left_join(parsed, index_updated_IND, by= "id")
write.table(parsed_IND, "JZLCK12_joined_parsed_IND.tsv", row.names=FALSE,sep="\t")


parsed <- read.delim("JZLCK21_NP311CD8_parsed_seqs_probs.tsv")
parsed_IND <-left_join(parsed, index_updated_IND, by= "id")
write.table(parsed_IND, "JZLCK21_joined_parsed_IND.tsv", row.names=FALSE,sep="\t")

### note 15/6: done joining, now need to "gate" real TCRs

p <- ggplot(data = parsed_IND, aes(x=X.561..582.15.A, y=X.640..670.14.A, colour = vb_rep, va=va_rep, label = id, CDR3a = cdr3a, CDR3b = cdr3b, subject = subject, Tet=X.561..582.15.A, TCRb=X.640..670.14.A)) +
  geom_point() + theme_linedraw() + scale_size(range = c(1,20)) + theme(text = element_text(size = 20)) + scale_fill_viridis_c()+
  scale_x_log10()+scale_y_log10()+labs(x = "tetramer PE", y = "TCRb APC")


ggplotly(p)


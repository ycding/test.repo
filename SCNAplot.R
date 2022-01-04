####to make GISTIC peak plots using thresholded GISTIC data
peak <-read.csv("X:/LAB/lab/Latina_breast_cancer/GISTIC/latina146/gistic_peaks_s5m7q05v2.csv")
thresholded <- read.delim("X:/LAB/lab/Latina_breast_cancer/GISTIC/latina146/seg5mar7q05broad/s5m7q05borad.all_thresholded.by_genes.txt")
gif <-read.csv("X:/LAB/lab/OvarianCancer/geneAssociation/GTF_withEntrezID.csv", head=T)
gif$Chrom <-as.numeric(gsub(gif$chrom, pattern = "chr", replacement = ""))
thres.gene <-merge(thresholded, gif, by.x="Locus.ID", by.y="ENTREZID")
amp.peak<-peak$Descriptor[which(peak$SCNA=="gain")]
#loss.peak<-peak$Descriptor[which(peak$SCNA=="loss")]

mytheme <- theme(strip.background =element_rect(fill="white"),
                 strip.text.x = element_text(size = 14, colour = "brown"),
                 plot.title=element_text(face="bold.italic",size="14", color="brown"),
                 axis.title=element_text(face="bold.italic",size=16, color="brown"),
                 axis.text=element_text(face="bold", size=12,color="darkblue"),
                 panel.background=element_rect(fill="white",color="darkblue"),
                 panel.grid.major.y=element_line(color="grey",linetype=1),
                 panel.grid.minor.y=element_line(color="grey",linetype=2),
                 panel.grid.minor.x=element_blank(),legend.position="top", 
                 legend.title = element_text(size=14, face="bold"),
                 legend.text = element_text(size=14, face="bold"))

#amplification plots
for (i in 1:length(amp.peak)) {
  i=28
  peak_name<- amp.peak[i]
  gene_name<-peak$candidate_gene[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
  #gene_name <-"MIR4728"
  TSS<-thres.gene$txStart[which(thres.gene$gene==gene_name)]
  reg.start <-peak$region_start[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
  reg.end <-peak$region_end[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
  chr<-peak$chr[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
  chr_nam<-paste("Chr",chr)
  Region<-subset(thres.gene, Chrom==chr & txStart > reg.start -35000000 & txEnd < reg.end +5000000)
  Region <-Region[!duplicated(Region$Gene.Symbol),]
  row.names(Region) <- Region$Gene.Symbol
  Region <-Region[order(Region$txStart),]
  Region$sum2 <- sapply(1:nrow(Region), function (i) length(grep(+2, Region[i, 4:149])))
  Region$sum_m1 <- sapply(1:nrow(Region), function (i) length(grep(-1, Region[i, 4:149])))
  Region$sum_m2 <- sapply(1:nrow(Region), function (i) length(grep(-2, Region[i, 4:149])))
  region <-Region[,c(155,156,164:166)]# generate input file with a few relevant variable for ggplot
  region$sum2 <- region$sum2 -region$sum_m2
  gene_plot <- ggplot(region) + 
    geom_area(aes(x=`txStart`, y=(`sum2`/1.46)), fill="red", stat="identity")  +
    #geom_area(aes(x=`txStart`, y=-((`sum_m1`+`sum_m2`)/1.46)), fill="blue", stat="identity") +  # not good to plot CN loss using this method
    geom_vline(xintercept=TSS) +
    ylim(0,max(region$sum2/1.46 +2)) + 
    xlab(chr_nam) + 
    scale_x_continuous(breaks = seq(from = min(region$txStart), to = max(region$txEnd), by = 2500000))+
    ylab("% CN Amplifications in 146 samples") +
    geom_text(mapping = aes(x = TSS, y = max(sum2/1.46 + 0.5), label = gene_name, hjust = -.5, vjust = -.5)) + mytheme
  png(paste0("X:/LAB/lab/Latina_breast_cancer/GISTIC/latina146/plot/SCNA2_5Mb_test", peak_name, gene_name, ".png"),width = 10*200,height = 7*200,res = 300,pointsize = 8 )       
  plot(gene_plot)
  dev.off()
  write.csv(region, paste0("X:/LAB/lab/Latina_breast_cancer/GISTIC/latina146/plot/SCNA2_gene_5Mb", peak_name,gene_name, ".csv"))
}

i=24
peak_name<- amp.peak[i]
gene_name<-peak$candidate_gene[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
gene_name1 <-"MIR4728"
TSS1<-thres.gene$txStart[which(thres.gene$gene==gene_name1)]
gene_name2 <-"ERBB2"
TSS2<-thres.gene$txStart[which(thres.gene$gene==gene_name2)]
gene_name3 <-"KIAA0100"
TSS3<-thres.gene$txStart[which(thres.gene$gene==gene_name3)]

gene_name4 <-"ZNF652"
TSS4<-thres.gene$txStart[which(thres.gene$gene==gene_name4)]
gene_name5 <-"PTRH2"
TSS5<-thres.gene$txStart[which(thres.gene$gene==gene_name5)]
gene_name6 <-"DDX5"
TSS6<-thres.gene$txStart[which(thres.gene$gene==gene_name6)]
gene_name7 <-"UBE2O"
TSS7<-thres.gene$txStart[which(thres.gene$gene==gene_name7)]

reg.start <-peak$region_start[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
reg.end <-peak$region_end[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
chr<-peak$chr[which(peak$Descriptor==peak_name & peak$SCNA=="gain")]
chr_nam<-paste("Chr",chr)
Region<-subset(thres.gene, Chrom==chr & txStart > reg.start - 15000000 & txEnd < reg.end + 45000000)
Region <-Region[!duplicated(Region$Gene.Symbol),]
row.names(Region) <- Region$Gene.Symbol
Region <-Region[order(Region$txStart),]
Region$sum2 <- sapply(1:nrow(Region), function (i) length(grep(+2, Region[i, 4:149])))
Region$sum_m1 <- sapply(1:nrow(Region), function (i) length(grep(-1, Region[i, 4:149])))
Region$sum_m2 <- sapply(1:nrow(Region), function (i) length(grep(-2, Region[i, 4:149])))
region <-Region[,c(155,156,164:166)]# generate input file with a few relevant variable for ggplot
region$sum2 <- region$sum2 -region$sum_m2
gene_plot <- ggplot(region) + 
  geom_area(aes(x=`txStart`, y=(`sum2`/1.46)), fill="red", stat="identity")  +
  #geom_area(aes(x=`txStart`, y=-((`sum_m1`+`sum_m2`)/1.46)), fill="blue", stat="identity") +  # not good to plot CN loss using this method
  geom_vline(xintercept=TSS1) +
  geom_vline(xintercept=TSS2, color="blue") +
  geom_vline(xintercept=TSS3) +
  geom_vline(xintercept=TSS4) +
  geom_vline(xintercept=TSS5) +
  geom_vline(xintercept=TSS6) +
  geom_vline(xintercept=TSS7) +
  ylim(0,max(region$sum2/1.46 +2)) + 
  xlab(chr_nam) + 
  scale_x_continuous(breaks = seq(from = min(region$txStart), to = max(region$txEnd), by = 10000000))+
  ylab("% CN Amplifications in 146 samples") +
  geom_text(mapping = aes(x = TSS1, y = 12.6, label = gene_name1, hjust = -.05, vjust = -.5),fontface="italic") + 
  geom_text(mapping = aes(x = TSS2-10000, y = 12.6,  label = gene_name2, hjust = +1.1, vjust = -.5),color="blue",fontface="italic") + 
  geom_text(mapping = aes(x = TSS3, y = 7.4, label = gene_name3, hjust = -.05, vjust = -.5),fontface="italic") +
  geom_text(mapping = aes(x = TSS4, y = 11, label = gene_name4, hjust = -.05, vjust = -.5),fontface="italic") + 
  geom_text(mapping = aes(x = TSS5, y = 12.6, label = gene_name5, hjust = +0.4, vjust = -.5),fontface="italic") + 
  geom_text(mapping = aes(x = TSS6, y = 11, label = gene_name6, hjust = -.05, vjust = -.5),fontface="italic") + 
  geom_text(mapping = aes(x = TSS7, y = 5, label = gene_name7, hjust = -.05, vjust = -.5),fontface="italic") + 
  mytheme




png(paste0("X:/LAB/lab/Latina_breast_cancer/GISTIC/latina146/plot/SCNA2_Chr17_6regions", ".png"),width = 10*200,height = 7*200,res = 300,pointsize = 8 )       
plot(gene_plot)
dev.off()
write.csv(region, paste0("X:/LAB/lab/Latina_breast_cancer/GISTIC/latina146/plot/SCNA2_gene_Chr17_6regions", ".csv"))
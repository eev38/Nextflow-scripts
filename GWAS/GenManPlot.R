#!/usr/bin/env Rscript

library(dplyr)
library(qqman)

# Set up variable to control command line arguments
args <- commandArgs(trailingOnly= TRUE)


GWASResults_t <- read.table(args[1], header = TRUE) %>%  
    mutate(rs = paste0("rs",chr, ps))

#FDR 

p.adj <- p.adjust(GWASResults_t$p_wald, method = "fdr")


GWASResults <- GWASResults_t %>%
    mutate(q_val = p.adj)

FDRThr <- GWASResults %>%
    filter(q_val <= 0.05)




### Extract top 100 hits from raw p-values

TopResults <- GWASResults %>%
    select(chr, rs, p_wald, q_val) %>%
         slice_min(p_wald, n=100)  %>% # smallest 100 p-values
         arrange(p_wald) #sort from lowest to highest p-value

write.csv(TopResults, "TopResults.csv", row.names = FALSE)


#optional add FDR threshold line to manhattan plot 
#for genomewideline max p-value within subset where qval <0.05 
pval <- -log10(max(FDRThr$p_wald))


#generate manhattan plot 

pdf("ManhattanPlot.pdf")
manhattan(GWASResults, chr="chr", bp="ps",  p="p_wald", snp ="rs", col=c("plum4", "powderblue"), genomewideline= pval, suggestiveline=FALSE) 
dev.off()



#generate qqplot

pdf("qqPlot.pdf")
qqplot <- qq(GWASResults$p_wald)
dev.off()


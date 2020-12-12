#########################################################
#written by Peihua Zhao
#mail: peihua.zhao@kuleuven.be or peihuazhao@outlook.com
#########################################################


LongRangeInteraction <- function(filename, vp.chr, vp.pos, window_size = 400, min_reads = 2, masked_region_length = 5000000){ 
          masked_region_start = vp.pos - masked_region_length
          masked_region_end = vp.pos + masked_region_length
          mydata <- read.table(filename , header = T)
          chromosomes=unique(mydata[,1])
          interaction.out=c()
          for (each_chr in chromosomes){ 
                  chr_mydata = mydata[mydata[,1]==each_chr,]
                  if (each_chr == vp.chr){
                              chr_mydata = chr_mydata[chr_mydata[,2] <  masked_region_start | chr_mydata[,2] > masked_region_end, ]
                  }
                  background_proportion = sum(chr_mydata[,3]>=min_reads) / dim(chr_mydata)[1]
                  pvalue.vector <- c()
                  window_start_pos.vector <- c()
                  window_end_pos.vector <- c()
                  observed_fragments.vector <- c()
                  for (i in seq(1,dim(chr_mydata)[1]-window_size,window_size)){
                             window_start=i
                             window_end=i + window_size -1
                             window_start_pos = chr_mydata[window_start , 2]
                             window_end_pos = chr_mydata[window_end , 2]
                             if  ((each_chr == vp.chr) & (window_start_pos <= masked_region_start) & (window_end_pos >= masked_region_end)){
                                        next
                             }
                             window_mydata = chr_mydata[window_start :  window_end, ]
                             observed_fragments= sum(window_mydata[,3]>=min_reads)
                             pvalue<-(binom.test(observed_fragments, window_size, p = background_proportion, alternative = c("greater"),))$p.value
                             pvalue.vector <- c(pvalue.vector , pvalue)
                             window_start_pos.vector <- c(window_start_pos.vector , window_start_pos)
                             window_end_pos.vector <- c(window_end_pos.vector , window_end_pos)
                             observed_fragments.vector <- c(observed_fragments.vector , observed_fragments)
                  }
                  each_chr_output=data.frame(Chr= each_chr, Interaction_start=window_start_pos.vector , Interaction_end=window_end_pos.vector , Number_Observed_Fragments =observed_fragments.vector  , Background_Proportion=background_proportion , P_Value=pvalue.vector) 
                  interaction.out=rbind(interaction.out , each_chr_output)
          }
          FDR=p.adjust(interaction.out$P_Value , method = "fdr")
          interaction.out=cbind(interaction.out , FDR)
          out.file <- paste(filename,".long-range.interaction.txt",sep="")
          write.table(interaction.out, out.file, quote=F, row=F, sep="\t")
}

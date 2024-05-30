rm(list =ls())
gc()

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


## merge table-----------------------------------
all_file = list.files(path = "./")
match_str = "_clean_star_ens_mm9_sort_featurecount_gene_countpair.txt"
all_file_target = all_file[grepl(pattern = match_str,all_file,fixed = T)]



for(i in 1:length(all_file_target)){
  file_name = all_file_target[i]
  sample_name = gsub(pattern = match_str,replacement = "",file_name)
  if(i ==1){
    sub_df =  read.csv(file = file_name,sep = "\t",skip = 1)
    colnames(sub_df)[7] = sample_name
  }else{
    sub_df1 = read.csv(file = file_name,sep = "\t",skip = 1)
    colnames(sub_df1)[7] =  sample_name
    sub_df[,sample_name] = sub_df1[,sample_name]
  }
 
}

saveRDS(sub_df,file = "star_ens_mm9_sort_featurecount_gene_countpair_rnaseq_count.rds")

calculate_fpkm <- function(count) {
  ## 每kb
   len = count[,1]/1000
   count_filt = count[,-1]
   
   out_list = list()
   for(j in 1:ncol(count_filt)){
     count_filt1 = count_filt[,j]
     ## 计算每个样本library size
     total = as.numeric(sum(count_filt1))
     len_total = len * total
     # 计算FPKM
     nor = (count_filt1*1e6)/len_total
     out_list[[j]] = nor
   }
   out_df =as.data.frame( do.call(cbind,out_list))
   colnames(out_df) = paste0(colnames(count_filt),"_FPKM")
   return(out_df)
}  

mat = calculate_fpkm(sub_df[,c(6,7,8)])

combined_df = cbind(sub_df,mat)

combined_df$simple_id = sub("\\..*", "", combined_df$Geneid)

saveRDS(combined_df,file = "star_ens_mm9_sort_featurecount_gene_countpair_rnaseq_count_fpkm.rds")

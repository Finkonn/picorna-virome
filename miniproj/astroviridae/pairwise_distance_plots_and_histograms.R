setwd("C:/MALICHAVA/Virome/picorna-virome/miniproj/astroviridae")

library('ape')
library('ggplot2')
library('stringr')
library('ggpubr')
library('spaa')


#Correspondence between nt and aa distances

# col_sp - position of species in sequence name
# OQ325671/pigeon/2025/Sp1 - in this example col_sp=4

plot_aa_nt_hists_sp = function(aln_object, pos_start, pos_end, col_sp=0, pairwise=TRUE){
  
  dna_slice = aln_object[1:length(aln_object[,1]), seq(from = pos_start, to = pos_end, by=1)]
  dna_slice_char=as.character.DNAbin(dna_slice)
  dna_slice_char[dna_slice_char=='-'] <- NA
  
  
  aa_slice = trans(dna_slice)
  aa_slice_char=as.character.AAbin(aa_slice)
  aa_slice_char[aa_slice_char=='X'] <- NA
  
  dist_nn = dist.gene(dna_slice_char, method = "percentage", pairwise.deletion = pairwise)
  dist_aa = dist.gene(aa_slice_char, method = "percentage", pairwise.deletion = pairwise)
  
  dist_nn= dist2list(dist_nn)
  dist_aa= dist2list(dist_aa)
  
  if (col_sp != 0){
    dist_nn = check_param(dist_nn,col_sp)
    colnames(dist_nn)[length(colnames(dist_nn))] = "same_species"
    d = merge(dist_nn, dist_aa, by=c("row", "col"))
    
    #scatterplot with marginal histograms
    p <- ggplot(d,aes(value.x,value.y, color=same_species)) + geom_point(alpha=0.1) + theme_bw() +
      xlab("nucleotide distance")+ylab("amino acid distance")
  }else{
    d = merge(dist_nn, dist_aa, by=c("row", "col"))
    p <- ggplot(d,aes(value.x,value.y)) + geom_point(alpha=0.1) + theme_bw() +
      xlab("nucleotide distance")+ylab("amino acid distance")
  }
  
  
  
  
  
  return(list(d, p))
  
}

# Plots distributions of pairwise nucleotide and amino acid distances

plot_hists = function(df_dist1, title1){
  
  if ("same_species" %in% colnames(df_dist1)){
    hist1_nt = ggplot(df_dist1, aes(x=value.x, fill = same_species)) +
      geom_histogram(alpha=0.5, bins=50) +
      theme_bw() + xlab("nucleotide distance") + ggtitle(title1) +
      theme(legend.position="none")+
      theme(legend.justification=c(0.1,0.9), legend.position=c(0.1,0.9))
    
    
    hist1_aa = ggplot(df_dist1, aes(x=value.y, fill = same_species)) +
      geom_histogram(alpha=0.5, bins=50) +
      theme_bw() + xlab("amino acid distance") + ggtitle(" ")  + theme(legend.position="none")
  }else{
    
    hist1_nt = ggplot(df_dist1, aes(x=value.x)) +
      geom_histogram(alpha=0.5, bins=50) +
      theme_bw() + xlab("nucleotide distance") + ggtitle(title1) +
      theme(legend.position="none")+
      theme(legend.justification=c(0.1,0.9), legend.position=c(0.1,0.9))
    
    
    hist1_aa = ggplot(df_dist1, aes(x=value.y)) +
      geom_histogram(alpha=0.5, bins=50) +
      theme_bw() + xlab("amino acid distance") + ggtitle(" ")  + theme(legend.position="none")
    
  }
  
  
  g1 = ggarrange(hist1_nt, hist1_aa,
                 ncol = 2, nrow = 1)
  
  
  
  return(g1)
}

# Checks whether values in specific position in sequence titles are equal or not
# Note - separator /

check_param = function(d, colnum){
  
  spl_df1 = str_split(d$row, '/', simplify = TRUE)
  spl_df2 = str_split(d$col, '/', simplify = TRUE)
  spl_df1_param = spl_df1[,colnum]
  spl_df2_param = spl_df2[,colnum]
  
  any_NA = spl_df1_param == "NA" | spl_df2_param == "NA"
  same_param = spl_df1_param == spl_df2_param
  
  vector_param = vector("character", length(spl_df1_param))
  vector_param[same_param == TRUE] = "yes"
  vector_param[same_param == FALSE] = "no"
  vector_param[any_NA == TRUE] = "unknown"
  d = cbind(d, vector_param)
  return(d)
}



# Example alignment

aln = read.dna("Astroviridae_pcr_1B_noamb_mammalia_aln_gr_only_ictv_headers.fna", format="fasta")

# species is not indicated in sequences names, so col_sp=0
ph_sp  = plot_aa_nt_hists_sp(aln, 1, length(aln[1,]), 5)
df_dist= ph_sp[[1]]

# Пары с nt p-distance > 0.2
#pairs_nt_gt02 <- df_dist[df_dist$value.x > 0.2, c("row", "col", "value.x", "value.y")]
pairs_nt_gt02 <- df_dist[, c("row", "col", "value.x", "value.y")]
pairs_nt_gt02 <- pairs_nt_gt02[order(-pairs_nt_gt02$value.x), ]
print(pairs_nt_gt02)

write.csv(pairs_nt_gt02, "pairs_pcr_orf1b_only_ICTV_headers.csv", row.names = FALSE)

# plot with correspondence beetween nt and aa distances
ph_sp[[2]]


hists = plot_hists(df_dist, "ORF1B_pcr")
hists

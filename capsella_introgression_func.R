# Functions for Capsella population definitions
# M. Wilson Brown
# July 10, 2024

##### Setup In/Out----------
# another option for pink color is #C5839A; I was using "pink" before (FFC0CB), #E2A2B3, or just 'mistyrose'
#anc.cols <- c(bursa_pastoris = "#a9ba9d", rubella = "#F1D1D9", heterozygous = "#96938E")
anc.cols <- c(bursa_pastoris = "#6D9636", rubella = "#FFCCFD", heterozygous = "#96938E")
# a warm gray or green grey need to be quite dark for it to be appropriate for all forms of colorblindness

# paths
ahmm_vit_path = "~/Documents/PhD/Research/capsella_introgression/ahmm_output/"
temp_plotdir <- "~/Documents/PhD/Research/capsella_introgression/alt_plots_temp/"
bed_files <- "~/Documents/PhD/Research/capsella_introgression/bed_files/"

##### Data Wrangling-------
bed_by_sample <- function(df, anc, population, outdir){
  # select specified ancestry paths
  tracts <- df[which(df$k3pop_sm %in% population & df$ancestry == anc),]
  # select only useful columns
  temp <- tracts[,c("vcf_sample_name","chrom","start","end", "ancestry")]
  
  # break up data frame by sample
  by_sample <- split(temp, temp$vcf_sample_name)
  
  # drop the sample id column from all
  by_sample <- lapply(by_sample, function(x){
    x[,-"vcf_sample_name"]
  })
  
  # simplify ancestry name if NY and NJ provided
  #ifelse(length(population > 1), population <- population[1])
  
  # write each to file
  for(i in 1:length(by_sample)){
    write.table(by_sample[[i]],
                paste0(outdir, #output directory
                       paste(names(by_sample)[i], #IDs
                             population,
                             anc, 
                             sep = "_"),
                       ".bed"), #file suffix
                col.names = F, quote = F, row.names = F, sep = "\t")
  }
  # return list
  return(by_sample)
}

###### Plotting Functions--------
# Viterbi decoded ancestry columns plot for a single scaffold
viterbi_columns_plot <- function(df, scaffold_num, population){
  pl <- ggplot() + geom_segment(data = subset(df, (chrom == paste0("jlSCF_", scaffold_num) & k3pop_sm %in% population)),
                          aes(color= ancestry, x=sample_name, xend=sample_name, y=start, yend=end),
                          linewidth = 8) +
    scale_color_manual(values=anc.cols) +
    ggtitle(paste(population, ": Scaffold", scaffold_num)) +
    ylab("position") + xlab("Sample Name") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))
  return(pl)
}

viterbi_columns_plot_allChr <- function(df, population){
  pl <- ggplot() + geom_segment(data = subset(df, k3pop_sm %in% population),
                                aes(color= ancestry, x=sample_name, xend=sample_name, y=start, yend=end),
                                linewidth = 8) +
    facet_grid(~chrom, scales = "free_y") +
    scale_color_manual(values=anc.cols) +
    ggtitle(paste(population," introgression patterns")) +
    ylab("position") + xlab("Sample Name") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))
  return(pl)
}

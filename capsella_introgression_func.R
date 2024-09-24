# Functions for Capsella introgression and local ancestry inference 
# M. Wilson Brown
# July 24, 2024

# scratch notes
# light.brown <- c("#f3cea9","#330000")
# greens <- c("#E5FFC5","#002C00")

##### Setup In/Out----------
# Colors for plotting
anc.cols <- c(bursa_pastoris = "#6D9636", rubella = "#FFCCFD", heterozygous = "#96938E")
# a warm gray or green grey need to be quite dark for it to be appropriate for all forms of colorblindness

species.cols <- c('Capsella\ grandiflora' = "#ffdb58",
               'Capsella\ rubella' = "#e99a9bff",
               'Capsella\ orientalis' = "#5d8aa8ff",
               'Capsella\ bursa-pastoris' = "#96938E",
               'E_Asia' = "#B2F2FD",
               'MENA' = "#996330",
               'N_Europe' = "#7FC55F",
               'NYC' = "#585380")

# label mapping for figures
pop.labels <- c('E_Asia' = "E. Asia",
                'MENA' = "Mediterranean",
                'N_Europe' = "N. Eurasia",
                'NYC' = "NYC Metro",
                'Capsella\ grandiflora' = "C. graniflora",
                'Capsella\ rubella' = "C. rubella",
                'Capsella\ orientalis' = "C. orientalis",
                'Capsella\ bursa-pastoris' = "C. bursa-pastoris")

# paths
ahmm_vit_path <- "~/Documents/PhD/Research/capsella_introgression/ahmm_output/"
temp_plotdir <- "~/Documents/PhD/Research/capsella_introgression/alt_plots_temp/"
plot_dir <- "~/Documents/PhD/Research/capsella_introgression/plots/"
bed_files <- "~/Documents/PhD/Research/capsella_introgression/bed_files/"

## pretty chromosome names for facets
pretty_chrom <- paste0("Chr. ", 9:16)
names(pretty_chrom) <- paste0("jlSCF_", 9:16)

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

viterbi_columns_plot_allChr <- function(df, population, facet.col = NULL, facet.txt.col = NULL){
  if(!is.null(facet.col)){
    pl <- ggplot() + geom_segment(data = subset(df, k3pop_sm %in% population),
                                aes(color= ancestry, x=sample_name, xend=sample_name, y=start, yend=end),
                                linewidth = 8) +
    facet_wrap(~chrom, scales = "free",
               labeller = labeller(chrom = pretty_chrom),
               nrow = 2) +
    scale_color_manual(values=anc.cols) +
    ggtitle(paste0(population," introgression patterns")) +
    ylab("position") + xlab("Sample Name") +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          strip.background =element_rect(fill= facet.col),
          strip.text = element_text(colour = facet.txt.col))
  } else {
    pl <- ggplot() + geom_segment(data = subset(df, k3pop_sm %in% population),
                                  aes(color= ancestry, x=sample_name, xend=sample_name, y=start, yend=end),
                                  linewidth = 8) +
      facet_wrap(~chrom, scales = "free",
                 labeller = labeller(chrom = pretty_chrom),
                 nrow = 2) +
      scale_color_manual(values=anc.cols) +
      ggtitle(paste0(population," introgression patterns")) +
      ylab("position") + xlab("Sample Name") +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))
  }
  return(pl)
}

# PIXY functions
pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read.delim(pixy_files[i], sep = "\t")
      df <- df %>%
        pivot_longer(cols = -c(pop, window_pos_1, window_pos_2, chromosome),
               names_to = "statistic", values_to = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read.delim(pixy_files[i], sep = "\t")
      df <- df %>%
        pivot_longer(cols = -c(pop1, pop2, window_pos_1, window_pos_2, chromosome),
                     names_to = "statistic", values_to = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}

# pixy Labeller
# custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)

# Functions for Capsella population definitions
# M. Wilson Brown
# July 10, 2024

# Viterbi decoded ancestry columns plot for a single scaffold
viterbi_columns_plot <- function(df, scaffold_num, k4vsNY_population){
  pl <- ggplot() + geom_segment(data = subset(df, (chrom == paste0("SCF_", scaffold_num) & k4vsNY_pop == k4vsNY_population)),
                          aes(color= ancestry, x=vcf_sample_name, xend=vcf_sample_name, y=start, yend=end),
                          linewidth = 8) +
    scale_color_manual(values=anc.cols) +
    ggtitle(paste(k4vsNY_population, ": Scaffold", scaffold_num)) +
    ylab("position") + xlab("Sample Name") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))
  return(pl)
}
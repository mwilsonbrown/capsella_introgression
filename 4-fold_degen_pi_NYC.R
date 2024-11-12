# Pi at 4-fold degenerate sites
# M. Wilson Brown
# November 12, 2024


# After running Pixy on all chromosomes separately, I concatenated all the data and then only selected rows for which there was a comparison made
# Libraries
library(ggplot2)

# Load in the data (~4 million lines)
pi <- read.delim("~/Documents/PhD/Research/capsella_introgression/pixy/nyc_allSites_CBP_NoNA.txt",
                 col.names = c("pop","chromosome","window_pos_1","window_pos_2","avg_pi","no_sites","count_diffs","count_comparisons","count_missing"))
# add an indicator for the subgenome each chromosome belongs to
pi$subgenome <- ifelse(pi$chromosome %in% c(paste0("jlSCF_", 1:8)), "CO", "CR")
                       
pi_poly <- ggplot() + geom_freqpoly(data = pi, aes(x=avg_pi, color = pop)) +
  theme_classic() +
  facet_grid(pop~subgenome)
ggsave(filename = "~/Documents/PhD/Research/capsella_introgression/plots/pi_poly.png",
       pi_poly,
       width = 5, height = 6, units = "in")

pi_pop <- ggplot() + geom_violin(data = pi, aes(x=pop, y=avg_pi, fill = subgenome),
                                 position = "dodge") + theme_classic() + scale_fill_manual(values = c("CO" = "lightblue",
                                                                                                      "CR" = "pink"))


# Caluclate average PI at 4 fold degenerate sites
# sum of the numerator and the sum of the denominator and divide
# Pi is the sum of the differences divided by the number of comparisons so

#for NYC
nyc_pi <- pi[which(pi$pop == "NYC"),]
# average Pi in NYC
sum(nyc_pi$count_diffs)/sum(nyc_pi$count_comparisons)

nyc_co <- nyc_pi[which(nyc_pi$chromosome %in% c(paste0("jlSCF_", 1:8))),]
nyc_cr <- nyc_pi[which(nyc_pi$chromosome %in% c(paste0("jlSCF_", 9:16))),]

sum(nyc_co$count_diffs)/sum(nyc_co$count_comparisons)
sum(nyc_cr$count_diffs)/sum(nyc_cr$count_comparisons)


#for N. Eurasia
neu_pi <- pi[which(pi$pop == "N_Europe"),]
# average Pi in NYC
sum(neu_pi$count_diffs)/sum(neu_pi$count_comparisons)

# Pi-4 is lower in NYC compared to the Northern Eurasian Population indicating a smaller effective population size

as_pi <- pi[which(pi$pop == "E_Asia"),]
# average Pi in NYC
sum(as_pi$count_diffs)/sum(as_pi$count_comparisons)

med_pi <- pi[which(pi$pop == "MENA"),]
# average Pi in NYC
sum(med_pi$count_diffs)/sum(med_pi$count_comparisons)



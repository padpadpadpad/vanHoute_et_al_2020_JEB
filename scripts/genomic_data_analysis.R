#-------------------------------------------------------------------------------------------------#
# analysis of genomic data from the pooled sequencing of compost and compost-water microcosms  ####
#-------------------------------------------------------------------------------------------------#

# load in packages ####

# clean workspace
rm(list = ls())

library(stringr)
library(vegan)
library(broom)
library(patchwork)
library(ggforce)
library(ggvegan) # remotes::install_github("gavinsimpson/ggvegan")
library(tidyverse)

# function for calculating distance from 00
dist_from_00 <- function(x, y){
  return(sqrt((0 - x)^2+(0-y)^2))
}

# load in data ####
d <- read.csv('data/snps_indels_data.csv', stringsAsFactors = FALSE) %>%
  select(., Chr, refseq, type, change_type, pos, ref, alt, annotation_impact, everything()) %>%
  mutate_at(., vars(starts_with('C_')), function(x)round(x, 1))

d_manhattan <- gather_(d, 'population', 'prop', colnames(d)[grepl('C_', colnames(d))]) %>%
  separate(., population, c('Something', 'bact_phage', 'mixed', 'rep'), sep = '_', remove = FALSE) %>%
  select(., -Something) %>%
  mutate(., mixed = plyr::mapvalues(.$mixed, from = c('S', 'M'), to = c('Compost', 'Compost-water')))

# manhattan plot of changes
manhattan <- ggplot() +
  geom_point(aes(start, prop), col = 'black', alpha = 0.5, filter(d_manhattan, change_type == 'none') %>% sample_n(10000)) +
  geom_point(aes(start, prop, col = mixed, shape = change_type), filter(d_manhattan, change_type != 'none' & prop >= 0.1), size = 5) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.8)) +
  ylab('proportion') +
  xlab('position in genome') +
  ylim(c(0,1))

# stack data frame ####

# keep only new genetic variants
d2 <- d[apply(d[grepl('C_', colnames(d))] > 0, 1, any),]

d_stack <- gather(d2, 'population', 'prop', starts_with('C_')) %>%
  separate(., population, c('Something', 'bact_phage', 'mixed', 'rep'), sep = '_', remove = FALSE) %>%
  select(., -Something) %>%
  mutate(., mixed = plyr::mapvalues(.$mixed, from = c('S', 'M'), to = c('Compost', 'Compost-water')))

# create summary stats for number of SNPs and distance from ancestor
d_SNP_sum <- d_stack %>%
  filter(., prop > 0) %>%
  group_by(., mixed, rep) %>%
  summarise(., dist_ancest = sum(prop),
            num_SNPs = n(),
            .groups = 'drop')

d_SNP_sum_mean <- d_SNP_sum %>%
  group_by(mixed) %>%
  summarise(across(c(dist_ancest, num_SNPs), .fns = list(mean = mean, sd = sd)), .groups = 'drop')

# plot 1 - Number of SNPs and distance from ancestor 
plot_dist <- ggplot(d_SNP_sum) +
  geom_linerange(aes(mixed, ymin = dist_ancest_mean - dist_ancest_sd, ymax = dist_ancest_mean + dist_ancest_sd), d_SNP_sum_mean) +
  geom_point(aes(mixed, dist_ancest_mean, shape = mixed), fill = 'black', size = 7, d_SNP_sum_mean) +
  geom_point(aes(mixed, dist_ancest, shape = mixed), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  ylab('genotypic distance from reference') +
  ggtitle('(a) genotypic distance') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        axis.title.x = element_blank()) +
  ylim(c(0, 0.8)) +
  scale_shape_manual(values = c(21, 24))

# plot 1b ####
plot_snpnum <- ggplot(d_SNP_sum) +
  geom_linerange(aes(mixed, ymin = num_SNPs_mean - num_SNPs_sd, ymax = num_SNPs_mean + num_SNPs_sd), d_SNP_sum_mean) +
  geom_point(aes(mixed, num_SNPs_mean, shape = mixed), fill = 'black', size = 7, d_SNP_sum_mean) +
  geom_point(aes(mixed, num_SNPs, shape = mixed), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  ylab('number of SNPs / indels') +
  ggtitle('(b) number of SNPs / indels') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        axis.title.x = element_blank()) +
  ylim(c(0,5)) +
  scale_shape_manual(values = c(21, 24))

# wilcoxon test for non-normal data
mod_dist <- wilcox.test(dist_ancest ~ mixed, d_SNP_sum, exact = FALSE)
mod_numsnps <- wilcox.test(num_SNPs ~ mixed, d_SNP_sum, exact = FALSE)

by(d_SNP_sum, d_SNP_sum$mixed, summary)

# calculate in how many replicates each variant occurs in
d_stack <- group_by(d_stack, name) %>%
  mutate(., pres = ifelse(prop > 0, 1, 0),
         tot_pres = sum(pres)) %>%
  ungroup()

d_stack_sum <- group_by(d_stack, name, prop, mixed, tot_pres) %>%
  tally() %>%
  ungroup()

p_genes <- ggplot(d_stack_sum, aes(forcats::fct_reorder(name, tot_pres, .desc = TRUE), prop, shape = mixed, group = mixed, size = n)) +
  geom_point(position = position_dodge(preserve = "total", width = 0.8), fill = 'white') +
  theme_bw(base_size = 12) +
  xlab('gene') +
  ylab('proportion') +
  scale_size_continuous(range = c(1,4)) +
  scale_shape_manual(values = c(21, 24), guide = FALSE) +
  ggtitle('(d) presence and proportion of SNPs/indels') +
  theme(legend.direction = 'horizontal',
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.7, 0.8)) +
  guides(size = guide_legend(title = 'Number of populations:', title.position = "top", override.aes = list(shape = 21, fill = 'white')))

p_genes

#------------#
# do nmds ####
#------------#

# use subset for only samples with mutations

#subset data for just the columns with samples in 
d_clust <- select(d2, colnames(d2)[grepl('C_', colnames(d2))])
row.names(d_clust) <- d2$name

# transpose the dataframe
d_clust <- t(d_clust)

# get variables for clustering
d_vars <- data.frame(label = row.names(d_clust)) %>%
  separate(., label, c('Something', 'bact_phage', 'mixed', 'rep'), sep = '_', remove = FALSE) %>%
  mutate(mixed = ifelse(mixed == 'S', 'Compost', 'Compost-water')) %>%
  select(., -c(Something, bact_phage))

# do nmds
nmds <- metaMDS(d_clust, distance = 'euclidean')

stressplot(nmds)
plot(nmds)

# get data from nmds
d_nmds <- fortify(nmds) %>%
  janitor::clean_names()

# wrangle sites
d_sample <- filter(d_nmds, score== 'sites') %>%
  merge(., d_vars, by = 'label')

# calculate centroid
d_sample <- group_by(d_sample, mixed) %>%
  mutate(centroid_nmds1 = mean(nmds1),
         centroid_nmds2 = mean(nmds2)) %>%
  ungroup()

# wrangle species
d_gene <- filter(d_nmds, score == 'species') %>%
  mutate(dist = dist_from_00(nmds1, nmds2))

# fancy biplot ####
p_nmds <- ggplot() +
  geom_text(aes(nmds1, nmds2, label = label, hjust = 0.5*(1 - sign(nmds1)), vjust = 0.5*(1-sign(nmds2))), col = 'red4', d_gene) +
  geom_segment(aes(x = 0, y = 0, yend = nmds2, xend = nmds1, group = score, alpha = dist), d_gene, arrow = arrow(length = unit(0.01, "npc")), col = 'red4') +
  geom_segment(aes(x = centroid_nmds1, y = centroid_nmds2, yend = nmds2, xend = nmds1, group = rep), d_sample) + 
  geom_point(aes(nmds1, nmds2, shape = mixed), d_sample , size = 3, fill = 'white') +
  geom_point(aes(centroid_nmds1, centroid_nmds2, shape = mixed), size = 5, data = select(d_sample, centroid_nmds1, centroid_nmds2, mixed) %>% distinct(), fill = 'black') +
  scale_x_continuous(expand = c(.15, .15)) +
  scale_y_continuous(expand = c(.15, .15)) +
  scale_alpha(guide = FALSE) +
  ggtitle('(e) NMDS of euclidean distances') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  scale_shape_manual(values = c(21,24))

# multidimensional scaling and permutational anova

# distance matrices
euclid_matrix <- dist(d_clust)

# make some character vectors factors
d_vars <- mutate_at(d_vars, c('mixed'), as.factor)
row.names(d_vars) <- d_vars$row

# adonis
mod_euclid <- adonis(euclid_matrix ~ mixed, data = d_vars, permutations = 9999)

# Looking at variance across groups
# can only do one factor so lets use id for every combination
mod_dispers <- betadisper(euclid_matrix, d_vars$mixed)

# plot of model
plot(mod_dispers)
boxplot(mod_dispers)

# anova
anova(mod_dispers)

# Permutation test for F
pmod <- permutest(mod_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod_dispers)

# alpha diversity - heterozygosity ####

# calculate Hardy Weinberg equilibrium for each SNP
alpha_div <- mutate(d_stack, p = prop,
                    q = 1-p,
                    div = 1 - p^2 - q^2,
                    nuc_div = q*p) %>%
  group_by(., mixed, rep) %>%
  summarise(., diversity = sum(div),
            nuc_div = sum(nuc_div)/n()) %>%
  data.frame()

alpha_div_mean <- group_by(alpha_div, mixed) %>%
  summarise(., across(diversity, list(mean = mean, sd = sd))) %>%
  ungroup()

# plot alpha_div ####
plot_alpha_div <- ggplot(alpha_div) +
  geom_linerange(aes(mixed, ymin = diversity_mean - diversity_sd, ymax = diversity_mean + diversity_sd), alpha_div_mean) +
  geom_point(aes(mixed, diversity_mean, shape = mixed), fill = 'black', size = 7, alpha_div_mean) +
  geom_point(aes(mixed, diversity, shape = mixed), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ggtitle('(c) alpha diversity') +
  xlab('') +
  ylab('alpha diversity') +
  theme(legend.position = 'none') +
  scale_shape_manual(values = c(21,24)) +
  ylim(c(0,1.1))

mod_div <- wilcox.test(diversity ~ mixed, alpha_div, exact = FALSE)
by(alpha_div, alpha_div$mixed, summary)

# plot
plot_dist + plot_snpnum + plot_alpha_div + plot_layout(ncol = 1) 

p_final <- ((plot_dist /plot_snpnum / plot_alpha_div) | (p_genes/p_nmds)) + plot_layout(widths = c(0.4, 0.6))
p_final
ggsave('figures/figure4.pdf', p_final, width = 8, height = 9)
ggsave('figures/figure4.png', p_final, width = 8, height = 9)
ggsave('figures/figure4.eps', p_final, width = 8, height = 9, dpi = 800, device=cairo_ps)

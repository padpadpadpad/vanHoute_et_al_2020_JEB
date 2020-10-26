#-----------------------------------------------------#
# analysis of resource-use data from biolog plates ####
#-----------------------------------------------------#

# recreates Figure 2

# load packages
library(tidyverse)
library(patchwork)
library(widyr)
library(corrr)
library(MicrobioUoE) # remotes::install_github('padpadpadpad/MicrobioUoE')

# plot path
fig_path <- 'figures'

# load in raw data ####
d_raw <- read.csv('data/biolog_plate_data.csv', na.strings = '')

# make NAs the last logical value
d_raw <- fill(d_raw, c(TREATMENT, REPLICATE), .direction = 'down') %>%
  rename(., treat = TREATMENT, rep = REPLICATE) %>%
  unite(., id, c(treat, rep, clone), sep = '_', remove = FALSE)

# turn into long format
Carb_cols <- colnames(d_raw)[grepl('X', colnames(d_raw))]
d_stack <- gather_(d_raw, 'C_source', 'OD', Carb_cols) %>%
  mutate(., C_source = as.numeric(gsub('X', '', C_source)))

# add column for mean OD per substrate
d_stack <- group_by(d_stack, C_source) %>%
  mutate(., mean = mean(OD),
         treat = tolower(treat)) %>%
  ungroup()

# filter out poor growth substrates
d_stack_poor <- filter(d_stack, C_source %in% c(21, 37, 63, 76))

ggplot(d_stack_poor, aes(OD)) +
  geom_histogram(fill = 'white', col = 'black') +
  facet_wrap(~C_source) +
  theme_bw()

# 76, 21, 63 & 37
d_stack_filt <- filter(d_stack, ! C_source %in% c(76, 21, 63, 37)) %>%
  mutate(., treat = ifelse(treat == 'static', 'Compost', 'Compost-water'))

# remove any other substrates where a clones had growth < 0.1 OD, no growth
d_stack_filt <- group_by(d_stack_filt, C_source) %>%
  filter(., min(OD) >= 0.1) %>%
  ungroup()

# carbon sources dropped
# 21 37 63 76
unique(d_stack_filt$C_source) %>% length()
# 92 substrates remain

# plot biolog plate data
p1 <- ggplot(d_stack_filt, aes(as.numeric(forcats::fct_reorder(factor(C_source), mean, .desc = TRUE)), OD)) +
  facet_wrap(~treat, labeller = labeller(treat = label_facets)) +
  geom_line(aes(group = interaction(rep, clone, treat)), alpha = 0.1) +
  stat_summary(aes(group = interaction(treat, rep)), fun = mean, geom = 'line') +
  xlab('Carbon source (ranked)') +
  ylab(expression(OD[660])) +
  theme_bw(base_size = 13*1.2) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

#----------------------------------#
# Calculate phenotypic variance ####
#----------------------------------#

# Phenotypic diversity, VP , is calculated as the average Euclidean distance across all pairs of clones in a population (Hall & Colegrave, 2006). For two genotypes from the same population, the Euclidean distance is the square root of the sum of the squared differences in OD between the two genotypes across all substrates. Biologically, this measures the differences in the metabolic profiles of two genotypes.


# Take all the distance matrices of each pairwise combination of clones in each evolved
# average Euclidean distance
pop_dists_df <- d_stack_filt %>%
  group_by(., treat, rep) %>%
  pairwise_dist(clone, C_source, OD, upper = FALSE) %>%
  rename(clone_i = item1, clone_j = item2)

# calculate average phenotypic diversity per population
V_P <- group_by(pop_dists_df, treat, rep) %>%
  summarise(., V_P = mean(distance)) %>%
  data.frame()
V_P_mean <- group_by(V_P, treat) %>%
  summarise(., V_P_mean = mean(V_P),
            sd = sd(V_P)) %>%
  data.frame()

#---------------------------------------#
# calculate V_G - genotypic variance ####
#---------------------------------------#

# Genotypic variance VG is calculated as the average variance of clone performance across all substrates (Venail et al., 2008).
V_G <- group_by(d_stack_filt, treat, C_source, rep) %>%
  summarise(V_G = var(OD)) %>%
  data.frame()
V_G_pop <- group_by(V_G, treat, rep) %>%
  summarise(V_G = mean(V_G)) %>%
  data.frame()

V_G_pop_mean <- group_by(V_G_pop, treat) %>%
  summarise(V_G_mean = mean(V_G),
            sd = sd(V_G)) %>%
  ungroup()

#-------------------------------------------#
# calculate V_E - environmental variance ####
#-------------------------------------------#

# Environmental variance, VE is calculated as the variance in the average clone performance across all substrates.

# calculate mean performance on each substrate per population
# calculate variance of mean performance across substrates
V_E <- group_by(d_stack_filt, treat, C_source, rep) %>%
  summarise(V_E = mean(OD)) %>%
  data.frame()
V_E_pop <- group_by(V_E, treat, rep) %>%
  summarise(V_E = var(V_E)) %>%
  data.frame()
V_E_pop_mean <- group_by(V_E_pop, treat) %>%
  summarise(V_E_mean = mean(V_E),
            sd = sd(V_E)) %>%
  ungroup()

# plot phenotypic, genotypic and environmental variance across treatments

# plot phenotypic variance
V_P_plot <- ggplot(V_P) +
  geom_point(aes(treat, V_P_mean, shape = treat), size = 7, fill = 'black', V_P_mean) +
  geom_linerange(aes(x = treat, ymin = V_P_mean - sd, ymax = V_P_mean + sd), V_P_mean) +
  geom_point(aes(treat, V_P, shape = treat), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  ylab('phenotypic variance') +
  xlab('') +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 12, color = 'black')) +
  ggtitle(expression((a)~Phenotypic~variance~(V[P]))) +
  scale_shape_manual(values = c(21, 24))

# plot V_G
V_G_plot <- ggplot(V_G_pop) +
  geom_point(aes(treat, V_G_mean, shape = treat), size = 7, fill = 'black', V_G_pop_mean) +
  geom_linerange(aes(x = treat, ymin = V_G_mean - sd, ymax = V_G_mean + sd), V_G_pop_mean) +
  geom_point(aes(treat, V_G, shape = treat), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  ylab('genotypic variance') +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        title = element_text(size = 12),
        axis.title.x = element_blank()) +
  ggtitle(expression((c)~Genotypic~variance~(V[G]))) +
  scale_shape_manual(values = c(21, 24))

# plot V_E
V_E_plot <- ggplot(V_E_pop) +
  geom_point(aes(treat, V_E_mean, shape = treat), size = 7, fill = 'black', V_E_pop_mean) +
  geom_linerange(aes(x = treat, ymin = V_E_mean - sd, ymax = V_E_mean + sd), V_E_pop_mean) +
  geom_point(aes(treat, V_E, shape = treat), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  ylab('environmental variance') +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        title = element_text(size = 12),
        axis.title.x = element_blank()) +
  ggtitle(expression((d)~Environmental~variance~(V[E]))) +
  scale_shape_manual(values = c(21, 24))

#----------------------------------------------------#
# Calculate G x E interaction for each population ####
#----------------------------------------------------#

# see Barrett et al. 2005 Am Nat and Venail et al. 2008 Nature

# 1. calculate responsiveness - indicates differences in the environmental variances and thus measures diversity of resource exploitation strategies (specialists and generalists)
# sum (sd_j - sd_i)^2/(2*n_genotypes(n_genotypes - 1))

# create dataframe for standard deviation per clone across environments
d_sd <- group_by(d_stack_filt, treat, rep, clone) %>%
  summarise(., sd_E = sd(OD)) %>%
  data.frame()

# create 2 copies of this for merging later
sd_j_clone <- dplyr::rename(d_sd, clone_j = clone, sd_j = sd_E) 
sd_i_clone <- rename(d_sd, clone_i = clone, sd_i = sd_E)

# create every pairwise combination of 1:n (clones/genotypes) for each population
d_R <- group_by(d_sd, treat, rep) %>%
  do(data.frame(expand.grid(clone_j = .$clone, clone_i = .$clone))) %>%
  ungroup() %>%
  filter(., clone_j > clone_i) %>%
  merge(., sd_j_clone, by = c('clone_j', 'treat', 'rep')) %>%
  merge(., sd_i_clone, by = c('clone_i', 'treat', 'rep'))

# calculate R for each pairwise combination
d_R <- group_by(d_R, treat, rep) %>%
  mutate(., R_comb = (sd_j - sd_i)^2/(2*n())*(n()-1)) %>%
  ungroup()

# calculate responsiveness for each population
# sum of all the pairwise combinations
d_R_pop <- group_by(d_R, treat, rep) %>%
  summarise(., R_pop = sum(R_comb)) %>%
  data.frame()

d_R_pop_mean <- group_by(d_R_pop, treat) %>%
  summarise(., R_pop_mean = mean(R_pop),
            sd = sd(R_pop)) %>%
  ungroup()

# plot responsiveness
r_plot <- ggplot(d_R_pop) +
  geom_point(aes(treat, R_pop_mean, shape = treat), size = 7, fill = 'black', d_R_pop_mean) +
  geom_linerange(aes(x = treat, ymin = R_pop_mean - sd, ymax = R_pop_mean + sd), d_R_pop_mean) +
  geom_point(aes(treat, R_pop, shape = treat), fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  ylab('responsiveness') +
  xlab('culture conditions') +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        title = element_text(size = 12)) +
  ggtitle(expression((e)~Responsiveness)) +
  scale_shape_manual(values = c(21, 24))

# 2. calculate inconsistency: Inconsistency indicates non-correlations between genotypes over different environments and indicated resource specialisation.

# calculate correlations between each pairs of clones within a population
d_pearson <- group_by(d_stack_filt, treat, rep) %>%
  nest() %>%
  mutate(., cor_col = purrr::map(data, pairwise_cor, item = clone, feature = C_source, value = OD, upper = FALSE)) %>%
  unnest(cor_col) %>%
  rename(clone_i = item1, clone_j = item2) 

# merge dataframe to responsiveness dataframe
d_inconsist <- merge(d_pearson, sd_i_clone, by = c('treat', 'rep', 'clone_i'), all.x = TRUE) %>%
  merge(., sd_j_clone, by = c('treat', 'rep', 'clone_j'), all.x = TRUE) %>%
  group_by(., treat, rep) %>%
  mutate(., i = (sd_j*sd_i*(1-correlation))/(n()*(n()-1))) %>%
  summarise(., I_pop = sum(i),
            pear_pop = mean(correlation)) %>%
  data.frame()

d_inconsist_mean <- group_by(d_inconsist, treat) %>%
  summarise(I_pop_mean = mean(I_pop),
            sd = sd(I_pop)) %>%
  ungroup()

# plot inconsistency
I_plot <- ggplot(d_inconsist) +
  geom_point(aes(treat, I_pop_mean, shape = treat), size = 7, fill = 'black', d_inconsist_mean) +
  geom_linerange(aes(x = treat, ymin = I_pop_mean - sd, ymax = I_pop_mean + sd), d_inconsist_mean) +
  geom_point(aes(treat, I_pop, shape = treat), fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  ylab('inconsistency') +
  xlab('culture conditions') +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none',
        title = element_text(size = 12)) +
  ggtitle(expression((f)~Inconsistency)) +
  scale_shape_manual(values = c(21, 24))

# stats tests ####

# 1. phenotypic diversity
mod_vp <- lm(V_P ~ treat, V_P)
mod_vp2 <- lm(V_P ~ 1, V_P)
anova(mod_vp, mod_vp2)
plot(mod_vp)

# 2. genotypic diversity
mod_vg <- lm(V_G ~ treat, V_G_pop)
mod_vg2 <- lm(V_G ~ 1, V_G_pop)
anova(mod_vg, mod_vg2)
plot(mod_vg)

# 3. environmental diversity
mod_ve <- lm(V_E ~ treat, V_E_pop)
mod_ve2 <- lm(V_E ~ 1, V_E_pop)
anova(mod_ve, mod_ve2)

# 4. responsiveness
mod_r <- lm(R_pop ~ treat, d_R_pop)
mod_r2 <- lm(R_pop ~ 1, d_R_pop)
anova(mod_r, mod_r2)

# 5. inconsistency
mod_i <- lm(I_pop ~ treat, d_inconsist)
mod_i2 <- lm(I_pop ~ 1, d_inconsist)
anova(mod_i, mod_i2)

# plot for paper
fig2 <- p1 + 
{V_G_plot + V_E_plot} +
{r_plot + I_plot} +
  plot_layout(ncol = 1,
              heights = c(0.4, 0.3, 0.3))

ggsave(file.path('figures', 'figure2.pdf'), fig2, height = 12, width = 9)
ggsave(file.path('figures', 'figure2.png'), fig2, height = 12, width = 9)
ggsave(file.path('figures', 'figure2.eps'), fig2, height = 12, width = 9, dpi = 800, device=cairo_ps)

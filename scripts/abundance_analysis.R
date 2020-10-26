#-------------------------------------------------------------------------#
# Analysis of abundance data from compost and compost-water microcosms ####
#-------------------------------------------------------------------------#

# recreates Figure 1

# load packages ####
library(tidyverse)
library(patchwork)
library(emmeans)

# load in abundance data
d <- read.csv('data/abundance_analysis.csv') %>%
  janitor::clean_names() %>%
  mutate(., treatment = tolower(treatment),
         treatment = ifelse(treatment == 'mixed', 'compost-water', 'compost'),
         # calculate total abundance by multiplying by volume of soil used
         cfu_total = ifelse(treatment == 'compost-water', cfu_g_soil * 3, cfu_g_soil * 25))

# calculate mean abundance - total and per g of soil
d_mean <- group_by(d, treatment) %>%
  summarise(cfu_mean = mean(cfu_g_soil),
            sd = sd(cfu_g_soil),
            cfu_total_mean = mean(cfu_total),
            cfu_sd = sd(cfu_total),
            .groups = 'drop')

# Make Figure 1b: abundance per gram of soil
p2 <- ggplot(d) +
  geom_point(aes(treatment, cfu_mean, shape = treatment), size = 7, fill = 'black', d_mean) +
  geom_linerange(aes(x = treatment, ymin = cfu_mean - sd, ymax = cfu_mean + sd), d_mean) +
  geom_point(aes(treatment, cfu_g_soil, shape = treatment), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^5, 10^8),
                minor_breaks = NULL) +
  theme_bw(base_size = 16) +
  xlab('Evolution environment') +
  ylab(expression(productivity~(cfu~g^-1~soil))) +
  scale_shape_manual(values = c(21, 24)) +
  theme(legend.position = 'none') +
  labs(title = '(b)')

# linear model
model <- lm(log10(cfu_g_soil) ~ treatment, d)
model2 <- lm(log10(cfu_g_soil) ~ 1, d)
anova(model, model2)

# calculate confidence intervals for each treatment
emmeans(model, pairwise~treatment)

# plot of total abundance
p1 <- ggplot(d) +
  geom_point(aes(treatment, cfu_total_mean, shape = treatment), size = 7, fill = 'black', d_mean) +
  geom_linerange(aes(x = treatment, ymin = cfu_total_mean - cfu_sd, ymax = cfu_total_mean + cfu_sd), d_mean) +
  geom_point(aes(treatment, cfu_total, shape = treatment), fill = 'white', size = 3, position = position_jitter(width = 0.1)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^5, 10^8),
                minor_breaks = NULL) +
  theme_bw(base_size = 16) +
  xlab('Evolution environment') +
  ylab(expression(total~bacterial~density~(cfu))) +
  scale_shape_manual(values = c(21, 24)) +
  theme(legend.position = 'none') +
  labs(title = '(a)')

# linear model
model <- lm(log10(cfu_total) ~ treatment, d)
model2 <- lm(log10(cfu_total) ~ 1, d)
anova(model, model2)

emmeans(model, pairwise~treatment)

p <- p1 + p2

ggsave('figures/figure1.pdf', p, height = 5, width = 10)
ggsave('figures/figure1.png', p, height = 5, width = 10)
ggsave('figures/figure1.eps', p, height = 5, width = 10, dpi = 800)

# calculate an estimate of the number of generations
days_total = 48
days_transfer = 6

# work out mean density at end
mean(d$cfu_total)

# calculate number of generations in total
(log(mean(d$cfu_total)) - log(mean(d$cfu_total)/3)) / log(2) * (days_total/days_transfer)


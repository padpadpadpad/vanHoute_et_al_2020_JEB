#--------------------------------------------------------------------------------#
# Analysis of relative fitness data from compost and compost-water microcosms ####
#--------------------------------------------------------------------------------#

# recreates Figure 3 and Table 1

# load packages ####
library(tidyverse)
library(lme4)
library(gt)
library(htmltools)
library(emmeans)

# function for facet labels
label_facets <- function(string){
  len <- length(string)
  string = paste('(', letters[1:len], ') ', string, sep = '')
  return(string)
}

# load in fitness assay data ####
d_fit <- read.csv('data/relative_fitness_analysis.csv') %>%
  janitor::clean_names() %>%
  mutate(., environment = tolower(environment)) %>%
  separate(., treatment, c('long_term_mixed', 'rep'), sep = 1) %>%
  mutate(., long_term_mixed = ifelse(long_term_mixed == 'S', 'compost evolved', 'compost-water evolved'),
         environment = ifelse(environment == 'static', 'compost', 'compost-water'))

# calculate relative fitness
d_fit <- mutate(d_fit,
                rel_fit = (fraction_white_1*(1 - fraction_white))/(fraction_white*(1 - fraction_white_1)))

# plot relative fitness
ggplot(d_fit) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_line(aes(environment, rel_fit, group = rep), size = 0.6) +
  geom_point(aes(environment, rel_fit, shape = long_term_mixed), fill = 'white', size = 4) +
  facet_wrap(~ long_term_mixed, labeller = labeller(long_term_mixed = label_facets)) +
  xlab ('competition environment') +
  ylab('relative fitness') +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = 'none') +
  scale_shape_manual(values = c(21, 24)) +
  ylim(c(0.5, 3.1))

ggsave('figures/figure3.pdf', last_plot(), height = 4, width = 8)
ggsave('figures/figure3.png', last_plot(), height = 4, width = 8)
ggsave('figures/figure3.eps', last_plot(), height = 4, width = 8, dpi = 800)

d_fit <- mutate(d_fit, pop = group_indices(d_fit, long_term_mixed, rep))

# lets do an analysis on these
d_fit_mod <- lmer(rel_fit ~ environment*long_term_mixed + (1|pop), d_fit)
d_fit_mod2 <- lmer(rel_fit ~ environment+long_term_mixed + (1|pop), d_fit)
anova(d_fit_mod, d_fit_mod2)

# look at pairwise differences
em_mod <- emmeans::emmeans(d_fit_mod, ~ long_term_mixed*environment)

contr_mat <- coef(pairs(em_mod))[, c('c.2', 'c.3', 'c.5')]

emmeans::emmeans(d_fit_mod, pairwise ~ long_term_mixed*environment)
emmeans::emmeans(d_fit_mod, pairwise ~ long_term_mixed)

emmeans(d_fit_mod, ~ long_term_mixed*environment, contr = contr_mat, adjust = 'Holm')

# try and make output table

d_table <- emmeans(d_fit_mod, ~ long_term_mixed*environment, contr = contr_mat, adjust = 'Hochberg')$contrasts %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(., p.value = round(p.value, 4)) %>%
  mutate_at(., c('estimate', 'SE', 'df', 't.ratio'), function(x){round(x, 2)}) %>%
  rename(`p value` = p.value, `t-ratio` = t.ratio, `d.f.` = df) %>%
  mutate(., contrast = case_when(contrast == 'c.2' ~ 'compost,compost - compost,compost-water',
                                 contrast == 'c.3' ~ 'compost,compost - compost-water,compost-water',
                                 contrast == 'c.5' ~ 'compost-water,compost - compost-water,compost-water'))

d_table

table <- gt(d_table) %>%
  cols_align('center') %>%
  tab_source_note(
    source_note = "P value adjustment: Holm-Bonferroni method for 3 tests."
  ) %>%
  tab_source_note(
    source_note = "Degrees-of-freedom method: kenward-roger."
  ) %>%
  tab_style(
    style = cell_text(
      weight = "bold"),
    locations = cells_body(
      rows = `p value` < 0.05)
  ) %>%
  tab_spanner(label = 'contrast', columns = 'contrast') %>%
  cols_label(contrast = 'evolution,competition') %>%
  gt:::as.tags.gt_tbl()

# change font
table[[1]]$children[[1]] <- gsub(
  "font-family: [[:print:]]*\n",
  "font-famuly: 'Times New Roman';\n",
  table[[1]]$children[[1]]
)


html_print(table, 'table_1.png')

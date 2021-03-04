## Generating options and choosing between them rely on distinct forms of value representation
## Adam Morris, Jonathan Phillips, Karen Huang, Fiery Cushman
## Analysis code for Studies 1 & 2
# Note: This code must be run with R version 3 (it's been specifically tested with R 3.6.3).
# It will break with R version >=4.

# setup -------------------------------------------------------------------

# load packages with groundhog (http://groundhogr.com/)
# start a new R session before doing this!
# if you get this (or any other) error:
# "groundhog says: 11 of the 21 packages needed by 'dplyr_0.8.4' are currently loaded, but not with the version that is needed."
# then run "rm(list=ls())", restart your R session, and try again.
#  if you still get the error, then do all the following steps in order:
# switch to R version 4, restart your R session, run "library(groundhog); groundhog.library('dplyr', '2020-06-01')", restart your R session, run "library(groundhog); groundhog.library('dplyr', '2020-03-01')", switch back to R version 3, restart your R session, and try running this script again.
# (I think this is a bug in groundhogr, and I have no idea why this fixes it, but that's what worked for me.)
library(groundhog)
pkgs = c('dplyr', 'tidyr', 'ggplot2', 'ggExtra', 'lme4', 'lmerTest', 'mlogit', 'stringdist', 'rje')
groundhog.library(pkgs, '2020-03-01')

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# update ggplot theme
theme_update(strip.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             plot.background = element_blank(),
             axis.text=element_text(size=30, colour = "black"),
             axis.title=element_text(size=18, face = "bold"),
             axis.title.x = element_text(vjust = 0),
             legend.title = element_text(size = 24, face = "bold"),
             legend.text = element_text(size = 20),
             plot.title = element_text(size = 26, face = "bold", vjust = 1),
             panel.margin = unit(1.0, "lines"), 
             plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
             axis.line = element_line(colour = "black", size = 2),
             axis.ticks = element_line(color = 'black', size = 3),
             axis.ticks.length = unit(.25, 'cm')
)

# make some handy functions
se = function(x) {return(sd(x, na.rm = T) / sqrt(sum(!is.na(x))))}
se.p = function(x) {return(sqrt(mean(x, na.rm = T) * (1 - mean(x, na.rm = T)) / sum(!is.na(x))))}
ci = function(m,s) {return(c(m - 1.96*s, m + 1.96*s))}
dodge <- position_dodge(width=0.9)
betterLine = function(data, formula, color = '#105db0') {
  lg = summary(lm(formula, data))$coefficients
  return(c(geom_abline(intercept = lg[1,1], slope = lg[2,1], color = color, size = 1.25),
           geom_abline(intercept = lg[1,1] - lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8),
           geom_abline(intercept = lg[1,1] + lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8)))
}

# This function is to compute standardized beta coefficients from a mixed effects model.
# It is copied verbatim from the package 'sjstats', version 0.17.2.
# The reason it's copied here (instead of loading the package) is that, when I did the analysis,
# I was using this latest version of sjstats, but with an old version of R (3.6.3).
# Groundhog won't let me use the latest version with an old version of R, and the older version was giving different answers,
# so I decided the simplest solution was just to copy the relevant function into my script.
# Citation:
# Lüdecke D (2019). _sjstats: Statistical Functions for Regression Models (Version 0.17.5)_. doi: 10.5281/zenodo.1284472 (URL: https://doi.org/10.5281/zenodo.1284472), <URL: https://CRAN.R-project.org/package=sjstats>.
std_beta <- function(fit, ci.lvl = .95, ...) {
  # compute ci, two-ways
  ci <- 1 - ((1 - ci.lvl) / 2)
  
  # code from Ben Bolker, see
  # http://stackoverflow.com/a/26206119/2094622
  sdy <- stats::sd(lme4::getME(fit, "y"))
  sdx <- apply(lme4::getME(fit, "X"), 2, sd)
  sc <- lme4::fixef(fit) * sdx / sdy
  se.fixef <- stats::coef(summary(fit))[, "Std. Error"]
  se <- se.fixef * sdx / sdy
  
  data_frame(
    term = names(lme4::fixef(fit)),
    std.estimate = sc,
    std.error = se,
    conf.low = sc - stats::qnorm(ci) * se,
    conf.high = sc + stats::qnorm(ci) * se
  ) %>%
    dplyr::slice(-1)
  
}

# choose which experiment -------------------------------------------------

expt = as.numeric(readline(prompt="Which study do you want to analyze? (Enter 1 or 2): "))

# import data -------------------------------------------------------------
expt.data = paste0('study', expt, '_data.csv')
expt.output = paste0('study', expt, '_output.rdata')

# import data from qualtrics readout
df = read.csv(expt.data, stringsAsFactors = F) %>% filter(DistributionChannel != 'preview', Status == 0) %>%
  dplyr::select(ResponseId, choice_30, choicetime_30_Page.Submit,
                cs_1, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_13,
                cs_time_Page.Submit,
                val_q_1, val_q_2, val_q_3, val_q_4, val_q_5, val_q_6, val_q_7, val_q_8, val_q_9,
                val_time_Page.Submit,
                val_other_q_1, val_other_q_2, val_other_q_3, val_other_q_4, val_other_q_5, val_other_q_6, val_other_q_7, val_other_q_8, val_other_q_9,
                Q159_Page.Submit)

# make a nicer data frame
# throughout this script, 'val' refers to general/context-free value and 'val.spec' refers to context-specific value
df2 = data.frame(subject = character(), food = character(), chosen = logical(), val = numeric(), val.spec = numeric(), choice = character())

cs.names = c('cs_1', 'cs_3', 'cs_4', 'cs_5', 'cs_6', 'cs_7', 'cs_8', 'cs_13')
val.names = c('val_q_1', 'val_q_2', 'val_q_3', 'val_q_4', 'val_q_5', 'val_q_6', 'val_q_7', 'val_q_8', 'val_q_9')
val.spec.names = c('val_other_q_1', 'val_other_q_2', 'val_other_q_3', 'val_other_q_4', 'val_other_q_5', 'val_other_q_6', 'val_other_q_7', 'val_other_q_8', 'val_other_q_9')
for (i in 1:nrow(df)) {
  subj = df$ResponseId[i]
  choice = as.character(df$choice_30[i])
  
  all.foods = as.character(df[i, cs.names])
  choice.ind = 1
  
  val = as.numeric(as.character(df[i, val.names[1]]))
  val.spec = as.numeric(as.character(df[i, val.spec.names[1]]))
  df2 = rbind(df2, data.frame(subject = subj, food = choice, chosen = 1, val = val, val.spec = val.spec, choice = choice))
  
  for (j in 1:length(cs.names)) {
    food = as.character(df[i, cs.names[j]])
    if (nchar(food) > 0) {
      val = as.numeric(as.character(df[i, val.names[j+1]]))
      val.spec = as.numeric(as.character(df[i, val.spec.names[j+1]]))
      chosen = 0
      df2 = rbind(df2, data.frame(subject = subj, food = food, chosen = chosen, val = val, val.spec = val.spec, choice = choice))
    }
  }
}

# exclude people
df.subj = df2 %>% group_by(subject) %>%
  summarize(cs.size = n(), nas = any(is.na(val) | is.na(val.spec)),
            var.choice = sd(c(val, val.spec)),
            val = mean(val), val.spec = mean(val.spec))

df2.filt = df2 %>% filter(!(subject %in% df.subj$subject[df.subj$nas == T | df.subj$var.choice == 0]))
df2.subj.filt = df.subj %>% filter(!nas & var.choice != 0)


# graphs and analysis -----------------------------------------------------

# check out value distributions of generated options
pct.gen.val.high = mean(df2.filt$val > 4)
pct.gen.val.high
ci(pct.gen.val.high, se.p(df2.filt$val > 4))
pct.gen.val.low = mean(df2.filt$val < 4)
pct.gen.val.low
ci(pct.gen.val.low, se.p(df2.filt$val < 4))
pct.gen.valspec.high = mean(df2.filt$val.spec > 4)
pct.gen.valspec.high
ci(pct.gen.valspec.high, se.p(df2.filt$val.spec > 4))
pct.gen.valspec.low = mean(df2.filt$val.spec < 4)
pct.gen.valspec.low
ci(pct.gen.valspec.low, se.p(df2.filt$val.spec < 4))

# check out value distributions of chosen options
pct.chosen.val.high = mean(df2.filt$val[df2.filt$chosen == 1] > 4)
pct.chosen.val.high
ci(pct.chosen.val.high, se.p(df2.filt$val[df2.filt$chosen == 1] > 4))
pct.chosen.val.low = mean(df2.filt$val[df2.filt$chosen == 1] < 4)
pct.chosen.val.low
ci(pct.chosen.val.low, se.p(df2.filt$val[df2.filt$chosen == 1] < 4))
pct.chosen.valspec.high = mean(df2.filt$val.spec[df2.filt$chosen == 1] > 4)
pct.chosen.valspec.high
ci(pct.chosen.valspec.high, se.p(df2.filt$val.spec[df2.filt$chosen == 1] > 4))
pct.chosen.valspec.low = mean(df2.filt$val.spec[df2.filt$chosen == 1] < 4)
pct.chosen.valspec.low
ci(pct.chosen.valspec.low, se.p(df2.filt$val.spec[df2.filt$chosen == 1] < 4))

# test how correlated general & specific value are
cor.test(df2.filt$val, df2.filt$val.spec)

# graph raw data (for supplement)
p = ggplot(data = df2.filt %>% mutate(chosen = ifelse(chosen == T, T, NA)),
           aes(x = val, y = val.spec, color = chosen, group = chosen)) +
  geom_jitter(size = 2) +
  geom_vline(xintercept = 4) +
  geom_hline(yintercept = 4) +
  labs(x = '', y = '') +
  guides(color = F)
ggMarginal(p, type='histogram', xparams = list(bins=8), yparams = list(bins = 8))

# test whether general value > midpoint of scale (4), controlling for specific value
m1 = lmer(val ~ val.spec + (val.spec | subject), data = df2.filt %>% mutate(val = val - 4, val.spec = val.spec - 4))
summary(m1)

# test whether specific value > midpoint of scale (4), controlling for general value
m2 = lmer(val.spec ~ val + (val | subject), data = df2.filt %>% mutate(val = val - 4, val.spec = val.spec - 4))
summary(m2)

## add info to df for implied probability analysis (reported in main text, described in supplement)

# add ranks
subjs = unique(df2.filt$subject)
for (i in 1:length(subjs)) {
  subj = subjs[i]
  df2.rows = df2.filt$subject == subj
  vals = -df2.filt$val[df2.rows]
  vals.spec = -df2.filt$val.spec[df2.rows]
  df2.filt$val.rank[df2.rows] = rank(vals, ties.method = 'min')
  df2.filt$val.spec.rank[df2.rows] = rank(vals.spec, ties.method = 'min')
}
df2.filt = df2.filt %>%
  mutate(val.rank = ifelse(val.rank > 5, 5, val.rank), val.spec.rank = ifelse(val.spec.rank > 5, 5, val.spec.rank),
         val.rank = max(val.rank) - val.rank + 1, val.spec.rank = max(val.spec.rank) - val.spec.rank + 1)

# assign every other midpoint to high or low (for dichotomization)
c1 = 0
c2 = 0
for (i in 1:nrow(df2.filt)) {
  if (df2.filt$val[i] == 4) {
    df2.filt$val.four[i] = c1
    c1 = ifelse(c1 == 0, 1, 0)
  } else {
    df2.filt$val.four[i] = -1
  }
  
  if (df2.filt$val.spec[i] == 4) {
    df2.filt$val.spec.four[i] = c2
    c2 = ifelse(c2 == 0, 1, 0)
  } else {
    df2.filt$val.spec.four[i] = -1
  }
}

df2.filt = df2.filt %>%
  mutate(val.adj = ifelse(val == 4, ifelse(val.four == 0, 3.9, 4.1), val),
         val.spec.adj = ifelse(val == 4, ifelse(val.spec.four == 0, 3.9, 4.1), val.spec),
         val.fac = factor(val.adj > 4, c(F,T), c('Low', 'High')),
         val.spec.fac = factor(val.spec.adj > 4, c(F,T), c('Low', 'High')))


## get inferred probabilities of inclusion in CS

# for graph, do marginal probability of each
graph.cs = df2.filt %>% select(subject, val.fac, val.spec.fac) %>%
  rename(val = val.fac, val.spec = val.spec.fac) %>%
  gather("type", "rank", -subject) %>%
  group_by(type, subject) %>%
  summarize(low = mean(rank == 'Low'), high = mean(rank == 'High')) %>%
  gather("rank", "prob", -c(subject,type)) %>%
  group_by(type, rank) %>%
  summarize(prob.m = mean(prob), prob.se = se(prob)) %>%
  group_by() %>%
  mutate(type = factor(type, c('val', 'val.spec')), rank = factor(rank, c('low', 'high')))
 
ggplot(graph.cs, aes(x = rank, y = prob.m, group = type, color = type)) +
  geom_point(size = 4) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = prob.m - prob.se, ymax = prob.m + prob.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(limits = c(.15,.85), breaks = c(.2, .8)) +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#d11a02', '#105db0'))

# for stats, do it more precisely: controlling for each other
df.model.cs = df2.filt %>% group_by(subject) %>%
  summarize(ll = mean(val.fac == 'Low' & val.spec.fac == 'Low'),
            lh = mean(val.fac == 'Low' & val.spec.fac == 'High'),
            hl = mean(val.fac == 'High' & val.spec.fac == 'Low'),
            hh = mean(val.fac == 'High' & val.spec.fac == 'High')) %>%
  gather("type", "prob", -subject) %>%
  mutate(val.high = type %in% c('hl', 'hh'), val.spec.high = type %in% c('lh', 'hh'))

m.cs = lmer(prob ~ val.high + val.spec.high + (1 | subject), data = df.model.cs)
summary(m.cs)
std_beta(m.cs)

## analyze choice from consideration set

# graph it
df.chosen.mf = df2.filt %>% group_by(val.rank, subject) %>%
  summarize(chosen = mean(chosen)) %>%
  group_by(val.rank) %>%
  summarize(chosen.m = mean(chosen), chosen.se = se(chosen))
df.chosen.mb = df2.filt %>% group_by(val.spec.rank, subject) %>%
  summarize(chosen = mean(chosen)) %>%
  group_by(val.spec.rank) %>%
  summarize(chosen.m = mean(chosen), chosen.se = se(chosen))
df.chosen = rbind(df.chosen.mf, df.chosen.mb %>% mutate(val.rank = val.spec.rank) %>% select(-val.spec.rank)) %>%
  mutate(type = rep(c(0,1), each = nrow(df.chosen.mf)), type.fac = factor(type, c(0,1), c('MF', 'MB')))
ggplot(df.chosen, aes(x = val.rank, y = chosen.m, color = type.fac, group = type.fac)) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2) +
  xlab('') +
  ylab('') +
  scale_color_manual(values = c('#d11a02', '#105db0')) +
  scale_y_continuous(breaks = c(0, .6), limits = c(0, .67)) +
  scale_x_continuous(labels = NULL, breaks = c(1,5)) +
  theme(legend.position = 'none')

# run multinomial logit
df2.logit = df2.filt %>% mutate(subject = as.numeric(subject))
for (i in 1:nrow(df2.logit)) {
  s = df2.logit$subject[i]
  if (i == 1 || s != last.s) {
    fc = 1
  }
  df2.logit$food.id[i] = fc
  fc = fc + 1
  last.s = s
}
df2.logit2 = mlogit.data(df2.logit, choice = "chosen", shape = "long", id.var = "subject", alt.var = "food.id", chid.var = 'subject')

m.selection = mlogit(chosen ~ val.rank + val.spec.rank | -1, df2.logit2, panel = T,
                     rpar = c(val.rank = "n", val.spec.rank = "n"), halton = NA, R = 1000, tol = .001)
m.selection.summary = summary(m.selection)
m.selection.summary

# For the mlogit model, there was no previous function to get the standardized coefficients, so I just did it myself.
# I'm not dividing by the sd of the dependent variable because it doesn't make sense for a multinomial logistic regression.
m.selection.summary$CoefTable[1,1] * sd(df2.logit2$val.rank)
m.selection.summary$CoefTable[1,2] * sd(df2.logit2$val.rank)
m.selection.summary$CoefTable[2,1] * sd(df2.logit2$val.spec.rank)
m.selection.summary$CoefTable[2,2] * sd(df2.logit2$val.spec.rank)

# get Bayes factor for null (used in Study 2)
m.selection.null = mlogit(chosen ~ val.spec.rank | -1, df2.logit2, panel = T,
                     rpar = c(val.spec.rank = "n"), halton = NA, R = 1000, tol = .001)
summary(m.selection.null)

ll1 = logLik(m.selection)
BIC1 = attr(ll1, 'df') * log(length(m.selection$fitted.values)) - 2 * as.numeric(ll1)

ll0 = logLik(m.selection.null)
BIC0 = attr(ll0, 'df') * log(length(m.selection.null$fitted.values)) - 2 * as.numeric(ll0)

BFnull.s1choice = exp((BIC1 - BIC0) / 2)
BFnull.s1choice
1 / BFnull.s1choice

# save --------------------------------------------------------------------

save.image(expt.output)

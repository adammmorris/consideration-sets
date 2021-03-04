## Generating options and choosing between them rely on distinct forms of value representation
## Adam Morris, Jonathan Phillips, Karen Huang, Fiery Cushman
## Analysis code for Study 8
# Note: This code must be run with R version 3 (it's been specifically tested with R 3.6.3).
# It will break with R version >=4.

# setup -------------------------------------------------------------------

# load packages with groundhog (http://groundhogr.com/)
# if you get this error:
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

se = function(x) {return(sd(x, na.rm = T) / sqrt(length(x)))}
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

# import data -------------------------------------------------------------

df = read.csv('study8_data.csv', stringsAsFactors = F) %>% filter(DistributionChannel != 'preview', Status == 0) %>%
  dplyr::select(ResponseId, choice_30, choicetime_30_Page.Submit, Q164, Q165_Page.Submit,
                cs_1, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_13,
                cs_time_Page.Submit,
                val_q_1, val_q_2, val_q_3, val_q_4, val_q_5, val_q_6, val_q_7, val_q_8, val_q_9, val_q_10, val_q_11,
                val_time_Page.Submit)


df2 = data.frame(subject = character(), food = character(), chosen = logical(), val = numeric(), choice = character(), cond = numeric())

choice.names = c('choice_30', 'Q164')
cond.names = c('choicetime_30_Page.Submit', 'Q165_Page.Submit')
cs.names = c('cs_1', 'cs_3', 'cs_4', 'cs_5', 'cs_6', 'cs_7', 'cs_8', 'cs_13')
val.names = c('val_q_1', 'val_q_3', 'val_q_4', 'val_q_5', 'val_q_6', 'val_q_7', 'val_q_8', 'val_q_9', 'val_q_10', 'val_q_11')
for (i in 1:nrow(df)) {
  subj = df$ResponseId[i]
  
  # figure out condition
  choice = NULL
  cond = NULL
  for (j in 1:length(cond.names)) {
    cond = as.character(df[i, cond.names[j]])
    if (nchar(cond) > 0) {
      choice = as.character(df[i, choice.names[j]])
      cond = j
      break
    }
  }
  
  all.foods = as.character(df[i, cs.names])
  
  # get num cs
  cs.size = sum(nchar(all.foods) > 0)
  
  val = as.numeric(as.character(df[i, val.names[cond]]))
  df2 = rbind(df2, data.frame(subject = subj, food = choice, chosen = 1, val = val, choice = choice, cond = cond))
  
  for (j in 1:length(cs.names)) {
    food = as.character(df[i, cs.names[j]])
    if (nchar(food) > 0) {
      val = as.numeric(as.character(df[i, val.names[j+2]]))
      chosen = 0
      df2 = rbind(df2, data.frame(subject = subj, food = food, chosen = chosen, val = val, choice = choice, cond = cond))
    }
  }
}

df2 = df2 %>% mutate(cond = factor(cond, 1:2, c('Same', 'Opposite')))


# do analysis -------------------------------------------------------------------

# exclude subjects
df.subj = df2 %>% group_by(subject) %>%
  summarize(var.choice = sd(val))

exclude = (df.subj %>% filter(var.choice == 0 | is.na(var.choice)))$subject

df2.filt = df2 %>% filter(!(subject %in% exclude)) %>%
  mutate(chosen = factor(chosen, c(0,1), c('Not chosen', 'Chosen')),
         val.opp = ifelse(cond == 'Same', val < 4, val > 4))

# get ranks
subjs = unique(df2.filt$subject)
for (i in 1:length(subjs)) {
  subj = subjs[i]
  df2.rows = df2.filt$subject == subj
  if (any(df2.filt$cond[df2.rows] == 'Same')) {
    vals = -df2.filt$val[df2.rows]
  } else {
    vals = df2.filt$val[df2.rows]
  }
  df2.filt$val.rank[df2.rows] = rank(vals, ties.method = 'min')
}
df2.filt = df2.filt %>% mutate(val.rank = ifelse(val.rank > 5, 5, val.rank),
                               val.rank = max(val.rank) - val.rank + 1)

df2.subj.filt = df.subj %>% filter(!(subject %in% exclude))

## analysis of "intrusions"
# graph
df2.test = df2.filt %>%
  group_by(cond, subject) %>%
  summarize(val = mean(val, na.rm = T), val.opp = mean(val.opp, na.rm = T),
            val.rank = mean(val.rank, na.rm = T)) %>%
  group_by(cond) %>%
  summarize(val.m = mean(val), val.se = se(val),
            val.opp.m = mean(val.opp), val.opp.se = se(val.opp),
            val.rank.m = mean(val.rank), val.rank.se = se(val.rank))

ggplot(df2.test, aes(x = cond, y = val.opp.m, fill = cond)) +
  geom_bar(stat = 'identity', position = dodge) +
  geom_errorbar(aes(ymin = val.opp.m - val.opp.se, ymax = val.opp.m + val.opp.se), width = .2, position = dodge, color = 'black') +
  ylab('') +
  xlab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = c(0,.2,.4)) +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('#8E44AD', '#58D68D'))

# stats
m = glmer(val.opp ~ cond + (1 | subject), family = 'binomial', data = df2.filt)
summary(m)
std_beta(m)

# graph for supplement
df.supp.graph = df2.filt %>%
  mutate(val.corrected = ifelse(cond == 'Same', val, max(val) - val + 1)) %>%
  group_by(cond, val.corrected) %>%
  summarize(n = n())

whichCond = 'Opposite'
ggplot(df.supp.graph %>% filter(cond == whichCond), aes(x = val.corrected, y = n)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 7), labels = c('', '')) +
  scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0,350)) +
  betterLine(df.supp.graph %>% filter(cond == whichCond), n ~ val.corrected)

## analysis of choice (for checking that people aren't more confused in opposite condition)
df2.chosen = df2.filt %>%
  group_by(cond, val.rank, subject) %>%
  summarize(chosen = mean(chosen == 'Chosen')) %>%
  group_by(cond, val.rank) %>%
  summarize(chosen.m = mean(chosen), chosen.se = se(chosen))
ggplot(df2.chosen, aes(x = val.rank, y = chosen.m)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2, position = dodge) +
  facet_wrap(~cond) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 5), labels = c('', '')) +
  scale_y_continuous(breaks = c(0, .6), labels = c(0, .6)) +
  betterLine(df2.chosen, chosen.m ~ val.rank)

ggplot(df2.chosen %>% filter(cond == 'Opposite'), aes(x = val.rank, y = chosen.m)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2, position = dodge) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 5), labels = c('', '')) +
  scale_y_continuous(breaks = c(0, .6), labels = c(0, .6), limits = c(0, .63)) +
  betterLine(df2.chosen %>% filter(cond == 'Opposite'), chosen.m ~ val.rank)

df2.filt %>% filter(chosen == 'Chosen') %>%
  group_by(cond, subject) %>%
  summarize(chosen = mean(val.rank == 5)) %>%
  group_by(cond) %>%
  summarize(chosen.m = mean(chosen), chosen.se = sqrt(chosen.m * (1-chosen.m) / n()))


# save --------------------------------------------------------------------

save.image('study8_output.RData')

## Generating options and choosing between them rely on distinct forms of value representation
## Adam Morris, Jonathan Phillips, Karen Huang, Fiery Cushman
## Analysis code for Studies 4-6: aggregated
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

getIndex = function(x, list) {
  y = numeric(length(x))
  for (j in 1:length(x)) {
    if (any(list %in% x[j])) {
      y[j] = which(list %in% x[j])
    } else {
      y[j] = NA
    }
  }
  return(y)
}

as.string.vector = function(x) {
  temp = strsplit(substr(x,2,nchar(x)-1), split=",")[[1]]
  return(substr(temp, 2, nchar(temp) - 1))
}

as.string.vector.noquotes = function(x) {
  temp = strsplit(substr(x,2,nchar(x)-1), split=",")[[1]]
  return(temp)
}

as.numeric.vector = function(x) {
  return(as.numeric(strsplit(substr(x,2,nchar(x)-1), split=",")[[1]]))
}

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

# load data ----------------------------------------------------------------

load('study4_output.rdata')

df1 = df.words.filt %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1_value_rank, s2_value_rank, s1value.fac, s2value.fac, in.cs, subject, chosen, word_ind)

load('study5_output.rdata')

df2 = df.words.filt %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1_value_rank, s2_value_rank, s1value.fac, s2value.fac, in.cs, subject, chosen, word_ind)

load('study6_output.rdata')

df3 = df.words.filt %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1_value_rank, s2_value_rank, s1value.fac, s2value.fac, in.cs, subject, chosen, word_ind)

df3.all = df.words.all %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1_value_rank, s2_value_rank, s1value.fac, s2value.fac, in.cs, subject, chosen, word_ind)

df = rbind(df1,df2,df3) %>% mutate(s1value.num = as.numeric(s1value.fac), s2value.num = as.numeric(s2value.fac))

df.selection = rbind(df1,df2,df3.all) %>% mutate(s1value.num = as.numeric(s1value.fac), s2value.num = as.numeric(s2value.fac))


# make df for multinomial logit
df.logit = df.selection %>% filter(in.cs == T) %>% select(subject, s1_value_rank, s2_value_rank, chosen, word_ind) %>%
  mutate(chosen = as.logical(chosen), intercept = 1) %>%
  rowwise() %>%
  mutate(subject.id = which(as.character(subject) == unique(df.selection$subject))) %>%
  group_by(subject.id) %>%
  mutate(option.id = row_number()) %>%
  mlogit.data(choice = "chosen", shape = "long", id.var = "subject.id", alt.var = "option.id", chid.var = "subject.id")


# analysis ----------------------------------------------------------------

## effect of stage 1 value

# on generation
df.s1 = df %>% group_by(s1value.num, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(s1value.num) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(in.cs))

# version for main text
ggplot(df.s1, aes(x = s1value.num, y = in.cs.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(.35,.45), limits = c(.32,.46)) + 
  scale_x_continuous(breaks = c(1,6), labels = c('', ''))
# version for supplement
ggplot(df.s1, aes(x = s1value.num, y = in.cs.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(.3,.45), limits = c(.3,.46)) + 
  scale_x_continuous(breaks = c(1,6), labels = c('1', '6')) +
  betterLine(df.s1, in.cs.m ~ s1value.num)

m.s1.cs = glmer(in.cs ~ s1value.num + (s1value.num | subject) + (s1value.num | word_ind), data = df, family = 'binomial')
summary(m.s1.cs)
std_beta(m.s1.cs)

# is there a quadratic component to the generation effect?
m.checkmark = glmer(in.cs~s1value.num+I(s1value.num^2)+(s1value.num+I(s1value.num^2)|subject)+(s1value.num+I(s1value.num^2)|word_ind),
                     data = df %>% mutate(s1value.num = s1value.num - mean(s1value.num)),
                     family='binomial')
summary(m.checkmark)
std_beta(m.checkmark)

# on choice (for main text)
df.s1.rank = df %>% group_by(s1_value_rank, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(s1_value_rank) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(in.cs))

ggplot(df.s1.rank, aes(x = s1_value_rank, y = chosen.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0,.8), limits = c(0,.82)) + 
  scale_x_continuous(breaks = c(1,5), labels = c('', ''))


## effect of stage 2 value
df.s2 = df %>% group_by(s2_value, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(s2_value) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(in.cs))

# on generation
ggplot(df.s2, aes(x = s2_value, y = in.cs.m)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .02) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(2, 25), labels = c(2, 25))+
  betterLine(df.s2, in.cs.m ~ s2_value)

m.s2.cs = glmer(in.cs ~ s2value.num + (s2value.num | subject) + (s2value.num | word_ind), data = df, family = 'binomial')
summary(m.s2.cs)
std_beta(m.s2.cs)

# on choice (for main text)
df.s2.rank = df %>% group_by(s2_value_rank, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(s2_value_rank) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(in.cs))

ggplot(df.s2.rank, aes(x = s2_value_rank, y = chosen.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0,.8), limits = c(0,.82)) + 
  scale_x_continuous(breaks = c(1,5), labels = c())

## choice analysis for supplement

# supplement graph
df.s1.rank.comb = df %>% group_by(s1_value_rank, s2_value_rank, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(s1_value_rank, s2_value_rank) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(in.cs))
ggplot(df.s1.rank.comb, aes(x = s2_value_rank, y = chosen.m, group = s1_value_rank, color = s1_value_rank)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .02) +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) + 
  scale_x_continuous(breaks = c(1,5), labels = c('1', '5')) +
  theme(legend.position = 'none')

# stats
m.choice = mlogit(chosen ~ s1_value_rank + s2_value_rank | -1 + word_ind, df.logit, panel = T,
                     rpar = c(s1_value_rank = "n", s2_value_rank = "n"), halton = NA, R = 1000, tol = .001)
m.choice.summary = summary(m.choice)
m.choice.summary

# For the mlogit model, there was no previous function to get the standardized coefficients, so I just did it myself.
# I'm not dividing by the sd of the dependent variable because it doesn't make sense for a multinomial logistic regression.
m.choice.summary$CoefTable[1,1] * sd(df.logit$s1_value_rank)
m.choice.summary$CoefTable[1,2] * sd(df.logit$s1_value_rank)
m.choice.summary$CoefTable[2,1] * sd(df.logit$s2_value_rank)
m.choice.summary$CoefTable[2,2] * sd(df.logit$s2_value_rank)


m.choice.null = mlogit(chosen ~ s2_value_rank | -1 + word_ind, df.logit, panel = T,
                  rpar = c(s2_value_rank = "n"), halton = NA, R = 1000, tol = .001)

ll1 = logLik(m.choice)
BIC1 = attr(ll1, 'df') * log(length(m.choice$fitted.values)) - 2 * as.numeric(ll1)

ll0 = logLik(m.choice.null)
BIC0 = attr(ll0, 'df') * log(length(m.choice.null$fitted.values)) - 2 * as.numeric(ll0)

BFnull = exp((BIC1 - BIC0) / 2)
BFnull

# save --------------------------------------------------------------------

save.image('study4-6_aggregated_output.rdata')

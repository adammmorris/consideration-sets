## Generating options and choosing between them rely on distinct forms of value representation
## Adam Morris, Jonathan Phillips, Karen Huang, Fiery Cushman
## Analysis code for Study 9
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

# import data -------------------------------------------------------------

expt = 9
cs.flipped = F
mid.value = 6
numTrials = 132
numQuestions = 2
numWords = 12;
minNAs = 1;

# Load data
df.demo = read.csv(paste0('study', expt, '_demo.csv'), stringsAsFactors = F) %>% arrange(subject) %>% mutate(total_time_real = total_time / 60000)
df.words.raw = read.csv(paste0('study', expt, '_words.csv'), stringsAsFactors = F) %>% arrange(subject, word_ind)
df.s1.raw = read.csv(paste0('study', expt, '_s1.csv'), stringsAsFactors = F) %>% arrange(subject)
df.s2.raw = read.csv(paste0('study', expt, '_s2.csv'), stringsAsFactors = F) %>% arrange(subject, question_order)

subjlist = df.demo$subject

# words
df.words = df.words.raw %>% filter(subject %in% subjlist) %>%
  mutate(repeated = word_ind == lead(word_ind)) %>% 
  filter(!repeated)

# s2
df.s2 = df.s2.raw %>% filter(subject %in% subjlist)
df.s2$choice = toupper(df.s2$choice)
df.s2$scratch = gsub("[.]", ",", toupper(as.character(df.s2$scratch)))
df.s2$all_values = as.character(df.s2$all_values)
df.s2$comp_check_pass = as.numeric(df.s2$comp_check_pass)

for (i in 1:nrow(df.s2)) {
  subj.name = df.s2$subject[i]
  wordlist = (df.words %>% filter(subject == subj.name))$word
  c = gsub("\n.*","",df.s2$choice[i])
  creal = wordlist[amatch(c, wordlist, maxDist = 2)]
  cind = getIndex(creal, wordlist)
  
  if (is.na(cind)) {
    # try scratch
    c = gsub("\n.*","",df.s2$scratch[i])
    creal = wordlist[amatch(c, wordlist, maxDist = 2)]
    cind = getIndex(creal, wordlist)
  }
  
  df.s2$choice_real[i] = creal
  df.s2$choice_real_ind[i] = cind
  
  df.s2$s2_value[i] = ifelse(is.na(cind), NA, df.words$s1_value[df.words$subject == subj.name & df.words$word == creal])
}

df.s2.subj = df.s2 %>% filter(subject %in% subjlist) %>%
  group_by(subject) %>%
  summarize(comp_check_pass = mean(comp_check_pass[question_order == 0]),
            numNAs = sum(is.na(choice_real)),
            numTrials = n(),
            cond = cond[1])
# s1
df.s1 = df.s1.raw %>% filter(subject %in% subjlist) %>% mutate(choice = as.numeric(choice))
df.s1$word_chosen = ifelse(df.s1$choice, df.s1$alt, df.s1$word)
df.s1$correct_choice = ifelse(df.s1$value == df.s1$value2, T, ifelse(df.s1$value > df.s1$value2, df.s1$choice == 0, df.s1$choice == 1))

df.s1.subj = df.s1 %>% group_by(subject) %>%
  summarize(pctCorrect_choice = mean(correct_choice, na.rm = T), numTrials = n())

# compute exclusion -------------------------------------------------------

include_rows = NULL
include_names = NULL

for (subj in 1:length(subjlist)) {
  subj.name = subjlist[subj]
  df.s1.subj.temp = df.s1.subj %>% filter(subject == subj.name)
  df.s1.temp = df.s1.subj %>% filter(subject == subj.name)
  df.s2.subj.temp = df.s2.subj %>% filter(subject == subj.name)
  df.demo.temp = df.demo %>% filter(subject == subj.name)
  df.words.temp = df.words %>% filter(subject == subj.name)
  
  exclude = df.demo.temp$write_down == 'Yes' || 
    df.demo.temp$use_s1 == 'Yes' || 
    df.s2.subj.temp$numNAs > minNAs ||
    df.s2.subj.temp$numTrials != numQuestions ||
    df.s1.subj.temp$numTrials != numTrials
  if (exclude) {
    include_rows[subj] = FALSE
  } else {
    include_rows[subj] = TRUE
    include_names = c(include_names, subj.name)
  }
}

# data prep --------------------------------------------------------------

df.s1.filt = df.s1 %>% filter(subject %in% include_names)

## words

df.words.filt = df.words %>% filter(subject %in% include_names)
for (i in 1:nrow(df.words.filt)) {
  
  subj = df.words.filt$subject[i]
  choice = (df.s2 %>% filter(subject == subj, question_order == 0))$choice_real
  df.s1.temp = df.s1 %>% filter(subject == subj)
  
  df.words.filt$freq[i] = sum(df.words.filt$word[i] == df.s1.temp$word_chosen)
  
  # get in.cs
  cs = (df.s2 %>% filter(subject == subj & question == 'choice-set'))$choice
  if (length(cs) > 0) {
    cs = as.string.vector.noquotes(gsub('\"', '', cs))
    for (j in 1:length(cs)) {
      cs.split = strsplit(cs[j], ":")[[1]]
      word = cs.split[1]
      val = cs.split[2]
      word_rows = df.words.filt$word == word & df.words.filt$subject == subj
      
      # in later versions, "1" means no and "0" means yes
      val = ifelse(cs.flipped, ifelse(val == "1", "0", ifelse(val == "0", "1", val)), val)
      
      if (length(choice) > 0) {
        df.words.filt$in.cs[word_rows] = ifelse(val == "1" | word == choice, T, ifelse(val == "0", F, NA))
      } else {
        df.words.filt$in.cs[word_rows] = NA
      }
    }
  } else {
    df.words.filt$in.cs[word_rows] = NA
  }
}

df.words.filt = df.words.filt %>% mutate(chosen = ifelse(in.cs, 0, NA),
                                         chosen_noNA = 0,
                                         high_s1value = factor(s1_value > mid.value, c(F,T), c('Low', 'High')),
                                         cond.fac = factor(cond, c(1,0), c('Best', 'Worst')),
                                         val.opp = ifelse(cond.fac == 'Best', s1_value < 6, s1_value > 6))

for (i in 1:nrow(df.words.filt)) {
  subj = df.words.filt$subject[i]
  df.temp = df.words.filt %>% filter(in.cs, subject == subj) %>% select(word_ind, s1_value, cond.fac) %>%
    mutate(s1_value = ifelse(cond.fac == 'Best', -s1_value, s1_value),
           val.rank = rank(s1_value, ties.method = 'min'))
  
  if (df.words.filt$in.cs[i]) {
    df.words.filt$s2_value_rank[i] = df.temp$val.rank[df.temp$word_ind == df.words.filt$word_ind[i]]
  } else {
    df.words.filt$s2_value_rank[i] = NA
  }
}


df.words.filt = df.words.filt %>% mutate(s2_value_rank = ifelse(s2_value_rank > 5, 5, s2_value_rank),
                                         s2_value_rank = max(s2_value_rank, na.rm = T) - s2_value_rank + 1)

## s2
df.s2.filt = df.s2 %>% filter(subject %in% include_names & question_order == 0)
for (i in 1:nrow(df.s2.filt)) {
  subj.name = df.s2.filt$subject[i]
  wordlist = (df.words.filt %>% filter(subject == subj.name))$word
  
  cind = df.s2.filt$choice_real_ind[i]
  creal = df.s2.filt$choice_real[i]
  
  word_rows = subj.name == df.words.filt$subject & creal == df.words.filt$word
  
  df.words.filt$chosen[word_rows] = 1
  df.words.filt$chosen_noNA[word_rows] = 1
}

df.demo.filt = df.demo %>% filter(subject %in% include_names)

# results ----------------------------------------------------------

# analysis of intrusions
df.intrusions = df.words.filt %>% filter(in.cs) %>%
  group_by(cond.fac, subject) %>%
  summarize(val.opp = mean(val.opp)) %>%
  group_by(cond.fac) %>%
  summarize(val.opp = mean(val.opp), val.opp.se = sqrt(val.opp * (1 - val.opp) / n()))

ggplot(df.intrusions, aes(x = cond.fac, y = val.opp, fill = cond.fac)) +
  geom_bar(stat = 'identity', position = dodge) +
  geom_errorbar(aes(ymin = val.opp - val.opp.se, ymax = val.opp + val.opp.se), width = .2, position = dodge, color = 'black') +
  theme(axis.title = element_text(size = 24)) +
  ylab('') +
  xlab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = c(0,.2,.4), limits = c(0,.4)) +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('#8E44AD', '#58D68D'))

m.intrusions = glmer(val.opp ~ cond.fac + (1 | subject), family = 'binomial', data = df.words.filt %>% filter(in.cs))
summary(m.intrusions)
std_beta(m.intrusions)

# first instance
m.intrusions1 = glmer(val.opp ~ cond.fac + (1 | subject), family = 'binomial', data = df.words.filt %>% filter(in.cs, id < 185741))
summary(m.intrusions1)
std_beta(m.intrusions1)
# second instance
m.intrusions2 = glmer(val.opp ~ cond.fac + (1 | subject), family = 'binomial', data = df.words.filt %>% filter(in.cs, id >= 185741))
summary(m.intrusions2)
std_beta(m.intrusions2)

# analysis of choice
df2.chosen = df.words.filt %>% filter(in.cs) %>%
  group_by(cond.fac, s2_value_rank, subject) %>%
  summarize(chosen = mean(chosen, na.rm = T)) %>%
  group_by(cond.fac, s2_value_rank) %>%
  summarize(chosen.m = mean(chosen), chosen.se = se(chosen))
ggplot(df2.chosen, aes(x = s2_value_rank, y = chosen.m)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2, position = dodge) +
  facet_wrap(~cond.fac) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 5), labels = c('', '')) +
  scale_y_continuous(breaks = c(0, .6), labels = c(0, .6)) +
  betterLine(df2.chosen, chosen.m ~ s2_value_rank)

ggplot(df2.chosen %>% filter(cond.fac == 'Best'), aes(x = s2_value_rank, y = chosen.m)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2, position = dodge) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 5), labels = c('', '')) +
  scale_y_continuous(breaks = c(0, .6), labels = c(0, .6), limits = c(0, .6)) +
  betterLine(df2.chosen %>% filter(cond.fac == 'Best'), chosen.m ~ s2_value_rank)
ggplot(df2.chosen %>% filter(cond.fac == 'Worst'), aes(x = s2_value_rank, y = chosen.m)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2, position = dodge) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 5), labels = c('', '')) +
  scale_y_continuous(breaks = c(0, .6), labels = c(0, .6), limits = c(0, .61)) +
  betterLine(df2.chosen %>% filter(cond.fac == 'Worst'), chosen.m ~ s2_value_rank)

df.words.filt %>% filter(chosen == 1) %>%
  group_by(cond.fac, subject) %>%
  summarize(s2_value_rank = mean(s2_value_rank == 5)) %>%
  group_by(cond.fac) %>%
  summarize(s2_value_rank = mean(s2_value_rank), s2_value_rank.se = sqrt(s2_value_rank * (1-s2_value_rank) / n()))

# graph for supplement
df.supp.graph = df.words.filt %>% filter(in.cs) %>%
  mutate(val.corrected = ifelse(cond.fac == 'Best', s1_value, max(s1_value) - s1_value + 1)) %>%
  group_by(cond.fac, val.corrected) %>%
  summarize(n = n())

whichCond = 'Worst'
ggplot(df.supp.graph %>% filter(cond.fac == whichCond), aes(x = val.corrected, y = n)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1, 12), labels = c('', '')) +
  scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0,160)) +
  betterLine(df.supp.graph %>% filter(cond.fac == whichCond), n ~ val.corrected)

# save --------------------------------------------------------------------

save.image(paste0('study', expt, '_output.rdata'))


# setup -------------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(tidyr)
require(lme4)
require(lmerTest)
require(mlogit)
require(stringdist)
require(rje)

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

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

# import data -------------------------------------------------------------

expt = '7'
cs.flipped = F 
mid.value = 6
numTrials = 132

numWords = 12
minNAs = 1
numQuestions = 2

# Load data
df.demo.raw = read.csv(paste0(expt, '/demo.csv'), stringsAsFactors = F) %>% arrange(subject) %>% mutate(total_time_real = total_time / 60000, stamp = as.Date(stamp))
df.words.raw = read.csv(paste0(expt, '/words.csv'), stringsAsFactors = F) %>% arrange(subject, word_ind)
df.s1.raw = read.csv(paste0(expt, '/s1.csv'), stringsAsFactors = F) %>% arrange(subject)
df.s2.raw = read.csv(paste0(expt, '/s2.csv'), stringsAsFactors = F) %>% arrange(subject, question_order)

subjlist = df.demo.raw$subject

df.demo = df.demo.raw %>% filter(subject %in% subjlist)

# words
df.words = df.words.raw %>% filter(subject %in% subjlist) %>%
  mutate(repeated = word_ind == lead(word_ind),
         s1_value = ifelse(cond == 1, s1_value, s1_value + 13)) %>%
  filter(!repeated)

# stage 2 prep
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
  
  df.s2$s2_value[i] = ifelse(is.na(cind), NA, df.words$s2_value[df.words$subject == subj.name & df.words$word == creal])
}

df.s2.subj = df.s2 %>% filter(subject %in% subjlist) %>%
  group_by(subject) %>%
  summarize(comp_check_pass = mean(comp_check_pass[question_order == 0]),
            numNAs = sum(is.na(choice_real)),
            numTrials = n(),
            cond = cond[1])

# stage 1 prep
df.s1 = df.s1.raw %>% filter(subject %in% subjlist) %>% mutate(choice = as.numeric(choice))
df.s1$word_chosen = ifelse(df.s1$choice, df.s1$alt, df.s1$word)
df.s1$correct_choice = ifelse(df.s1$value == df.s1$value2, T, ifelse(df.s1$value > df.s1$value2, df.s1$choice == 0, df.s1$choice == 1))

df.s1.subj = df.s1 %>% group_by(subject, cond) %>%
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
  
  exclude = df.demo.temp$write_down == 'Yes' || df.s2.subj.temp$comp_check_pass < 1 ||
    df.s2.subj.temp$numNAs > minNAs ||
    df.s2.subj.temp$numTrials != numQuestions ||
    subj.name == 'A1ZZ1RAYRQXPBG' || # someone who glitched
    df.s1.subj.temp$numTrials != numTrials ||
    df.s1.temp$pctCorrect_choice < .7
  if (exclude) {
    include_rows[subj] = FALSE
  } else {
    include_rows[subj] = TRUE
    include_names = c(include_names, subj.name)
  }
}

# data prep --------------------------------------------------------------

df.s1.filt = df.s1 %>% filter(subject %in% include_names)
df.s1.subj.filt = df.s1.subj %>% filter(subject %in% include_names)

## prep words df

df.words.filt = df.words %>% filter(subject %in% include_names)
for (i in 1:nrow(df.words.filt)) {
  
  subj = df.words.filt$subject[i]
  choice = (df.s2 %>% filter(subject == subj, question_order == 0))$choice_real
  df.s1.temp = df.s1 %>% filter(subject == subj)
  
  df.words.filt$freq[i] = sum(df.words.filt$word[i] == df.s1.temp$word_chosen)
  
  # get in.cs
  cs = (df.s2 %>% filter(subject == subj & question == 'choice-set'))$choice
  order = (df.s2 %>% filter(subject == subj & question == 'choice-set'))$scratch
  if (length(cs) > 0) {
    cs = as.string.vector.noquotes(gsub('\"', '', cs))
    order = as.string.vector.noquotes(gsub('\"', '', order))
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
      
      order.split = strsplit(order[j], ":")[[1]]
      order.word = order.split[1]
      order.val = order.split[2]
      word_rows = df.words.filt$word == order.word & df.words.filt$subject == subj
      
      df.words.filt$order[word_rows] = as.numeric(order.val)
    }
  } else {
    df.words.filt$in.cs[word_rows] = NA
    df.words.filt$order[word_rows] = NA
  }
}

df.words.filt = df.words.filt %>% mutate(chosen = ifelse(in.cs, 0, NA), chosen_noNA = 0,
                                         high_s1value = factor(s1_value > mid.value, c(F,T), c('Low', 'High')),
                                         high_s2value = factor(s2_value > 16, c(F,T), c('Low', 'High')))

# get ranks
for (i in 1:nrow(df.words.filt)) {
  if (!is.na(df.words.filt$in.cs[i]) && df.words.filt$in.cs[i] == T) {
    subj = df.words.filt$subject[i]
    df.words.temp = df.words.filt %>% filter(subject == subj, in.cs)
    
    s1.vals = df.words.temp$s1_value
    s2.vals = df.words.temp$s2_value
    
    s1.vals.rank = rank(-s1.vals, ties.method = 'min')
    s2.vals.rank = rank(-s2.vals, ties.method = 'min')
    
    cur.s1.rank = s1.vals.rank[s1.vals == df.words.filt$s1_value[i]]
    cur.s2.rank = s2.vals.rank[s2.vals == df.words.filt$s2_value[i]]
    df.words.filt$s1_value_rank[i] = ifelse(cur.s1.rank >= 5, 5, cur.s1.rank)
    df.words.filt$s2_value_rank[i] = ifelse(cur.s2.rank >= 5, 5, cur.s2.rank)
  } else {
    df.words.filt$s1_value_rank[i] = NA
    df.words.filt$s2_value_rank[i] = NA
  }
}

df.words.filt = df.words.filt %>% mutate(s1_value_rank = max(s1_value_rank, na.rm = T) - s1_value_rank + 1,
                                         s2_value_rank = max(s2_value_rank, na.rm = T) - s2_value_rank + 1)
df.words.filt = df.words.filt %>% mutate(cond.fac = factor(cond, c(1,0), c('Gains', 'Losses')))
df.words.filt = df.words.filt %>% mutate(s1_value_fac = cut(s1_value, 6, 1:6))

## stage 2 data frame prep
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

# create data frame for multinomial logit (for testing what influences choice out of consideration set)

df.logit = df.words.filt %>% filter(in.cs == T) %>% select(subject, s1_value, s2_value, s1_value_rank, s2_value_rank, chosen, word_ind) %>%
  mutate(chosen = as.logical(chosen), intercept = 1) %>%
  rowwise() %>%
  mutate(subject.id = which(as.character(subject) == unique(df.words.filt$subject))) %>%
  group_by(subject.id) %>%
  mutate(option.id = row_number()) %>%
  mlogit.data(choice = "chosen", shape = "long", id.var = "subject.id", alt.var = "option.id", chid.var = "subject.id")
   
# results ----------------------------------------------------------

## graph predictions of different hypotheses

# predictions of "high value" hypothesis
graph.s1.binary.hv = graph.s1.binary
graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Losses' & graph.s1.binary$high_s1value == 'Low'] = graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'Low']
graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Losses' & graph.s1.binary$high_s1value == 'High'] = graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'High']

ggplot(graph.s1.binary.hv, aes(x = high_s1value, y = in.cs.m, group = 1, linetype = cond.fac)) +
  geom_point(size = 5) +
  stat_summary(fun.y=sum, geom="line", size = 1.5) +
  #geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = cs.s1.binary.breaks, limits = cs.s1.binary.lims) + 
  theme(legend.position = 'none') +
  facet_wrap(~cond.fac) +
  theme(strip.text = element_blank())

# predictions of "absolute value" hypothesis
graph.s1.binary.av = graph.s1.binary
graph.s1.binary.av$in.cs.m[graph.s1.binary$cond.fac == 'Losses' & graph.s1.binary$high_s1value == 'High'] = graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'Low']
graph.s1.binary.av$in.cs.m[graph.s1.binary$cond.fac == 'Losses' & graph.s1.binary$high_s1value == 'Low'] = graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'High']

ggplot(graph.s1.binary.av, aes(x = high_s1value, y = in.cs.m, group = 1, linetype = cond.fac)) +
  geom_point(size = 5) +
  stat_summary(fun.y=sum, geom="line", size = 1.5) +
  #geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = cs.s1.binary.breaks, limits = cs.s1.binary.lims) + 
  theme(legend.position = 'none') +
  facet_wrap(~cond.fac) +
  theme(strip.text = element_blank())

# predictions of mixture
graph.s1.binary.mix = graph.s1.binary
graph.s1.binary.mix$in.cs.m[graph.s1.binary$cond.fac == 'Losses' & graph.s1.binary$high_s1value == 'High'] = (graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'Low'] + graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'High']) / 2
graph.s1.binary.mix$in.cs.m[graph.s1.binary$cond.fac == 'Losses' & graph.s1.binary$high_s1value == 'Low'] = (graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'Low'] + graph.s1.binary.hv$in.cs.m[graph.s1.binary$cond.fac == 'Gains' & graph.s1.binary$high_s1value == 'High']) / 2

ggplot(graph.s1.binary.mix, aes(x = high_s1value, y = in.cs.m, group = 1, linetype = cond.fac)) +
  geom_point(size = 5) +
  stat_summary(fun.y=sum, geom="line", size = 1.5) +
  #geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = cs.s1.binary.breaks, limits = cs.s1.binary.lims) + 
  theme(legend.position = 'none') +
  facet_wrap(~cond.fac) +
  theme(strip.text = element_blank())


## effects of stage 1 value -> generation

cs.s1.binary.breaks = c(.35, .4)
cs.s1.binary.lims = c(.34, .41)

# effect on generation -- dichotomized
graph.s1.binary = df.words.filt %>% group_by(cond.fac, high_s1value, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(cond.fac, high_s1value) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(chosen))

ggplot(graph.s1.binary, aes(x = high_s1value, y = in.cs.m, group = 1, linetype = cond.fac)) +
  geom_point(size = 5) +
  stat_summary(fun.y=sum, geom="line", size = 1.5) +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2, linetype = 1) +
  xlab('') + ylab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = cs.s1.binary.breaks, limits = cs.s1.binary.lims) + 
  theme(legend.position = 'none') +
  facet_wrap(~cond.fac) +
  theme(strip.text = element_blank())

# effect on generation -- undichotomized
graph.s1.full = df.words.filt %>% group_by(cond.fac, s1_value) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), in.cs.se = sqrt(in.cs * (1-in.cs) / n()),
            chosen = mean(chosen, na.rm = T), chosen.se = sqrt(chosen * (1-chosen) / n()),
            order.m = mean(order, na.rm = T), order.se = se(order))

ggplot(graph.s1.full, aes(x = s1_value, y = in.cs)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs - in.cs.se, ymax = in.cs+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(min(graph.s1.full$s1_value), max(graph.s1.full$s1_value))) +
  scale_y_continuous(breaks = c(.3,.45), limits = c(.3,.45)) +
  facet_wrap(~cond.fac) +
  betterLine(graph.s1.full, in.cs ~ s1_value)

ggplot(graph.s1.full %>% filter(cond.fac == 'Gains'), aes(x = s1_value, y = in.cs)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs - in.cs.se, ymax = in.cs+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(min(graph.s1.full$s1_value), max(graph.s1.full$s1_value))) +
  scale_y_continuous(breaks = c(.3,.45), limits = c(.3,.45)) +
  betterLine(graph.s1.full, in.cs ~ s1_value)

ggplot(graph.s1.full %>% filter(cond.fac == 'Losses'), aes(x = s1_value, y = in.cs)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs - in.cs.se, ymax = in.cs+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(min(graph.s1.full$s1_value), max(graph.s1.full$s1_value))) +
  scale_y_continuous(breaks = c(.3,.45), limits = c(.3,.45)) +
  betterLine(graph.s1.full, in.cs ~ s1_value)

# stats for generation
m.cs.s1 = glmer(in.cs~s1_value*cond+(s1_value|subject) + (s1_value*cond|word_ind),
               data = df.words.filt %>% mutate(cond = cond - .5, s1_value = s1_value - mean(s1_value)),
               family='binomial')
summary(m.cs.s1)

m.cs.s1.gains = glmer(in.cs~s1_value+(s1_value|subject) + (s1_value|word_ind),
                data = df.words.filt %>% filter(cond.fac == 'Gains'),
                family='binomial')
summary(m.cs.s1.gains)

m.cs.s1.losses = glmer(in.cs~s1_value+(s1_value|subject) + (s1_value|word_ind),
                      data = df.words.filt %>% filter(cond.fac == 'Losses'),
                      family='binomial')
summary(m.cs.s1.losses)

## effects of stage 2 value -> generation

# effect on generation
graph.s2 = df.words.filt %>% group_by(cond.fac, s2_value, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(cond.fac, s2_value) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(chosen))

ggplot(graph.s2, aes(x = s2_value, y = in.cs.m)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  labs(x = '', y = '')+
  scale_y_continuous(breaks = c(.2, .5), limits = c(.18, .57)) +
  scale_x_continuous(breaks = c(2, 25), limits = c(2,25)) +
  betterLine(graph.s2, in.cs.m ~ s2_value) +
  facet_wrap(~cond.fac)

ggplot(graph.s2 %>% filter(cond.fac == 'Gains'), aes(x = s2_value, y = in.cs.m)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  labs(x = '', y = '')+
  scale_y_continuous(breaks = c(.2, .5), limits = c(.18, .57)) +
  scale_x_continuous(breaks = c(2, 25), limits = c(2,25)) +
  betterLine(graph.s2, in.cs.m ~ s2_value)

ggplot(graph.s2 %>% filter(cond.fac == 'Losses'), aes(x = s2_value, y = in.cs.m)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  labs(x = '', y = '')+
  scale_y_continuous(breaks = c(.2, .5), limits = c(.18, .57)) +
  scale_x_continuous(breaks = c(2, 25), limits = c(2,25)) +
  betterLine(graph.s2, in.cs.m ~ s2_value)

# stats for generation
m.cs.s2 = glmer(in.cs~s2_value*cond+(s2_value|subject) + (cond|word_ind),
                data = df.words.filt %>% mutate(cond = cond - .5, s2_value = s2_value - mean(s2_value)),
                family='binomial')
summary(m.cs.s2)

## effects on choice

# combined
graph.choice = df.words.filt %>% group_by(cond.fac, s1_value_rank, s2_value_rank) %>% filter(in.cs == T) %>%
  summarize(chosen = mean(chosen, na.rm = T), chosen.se = sqrt(chosen * (1-chosen) / n()))
ggplot(graph.choice, aes(x = s2_value_rank, y = chosen, group = s1_value_rank, color = s1_value_rank)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = chosen - chosen.se, ymax = chosen+chosen.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1,5), labels = c(1,5)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  facet_wrap(~cond.fac)

ggplot(graph.choice %>% filter(cond.fac == 'Gains'), aes(x = s2_value_rank, y = chosen, group = s1_value_rank, color = s1_value_rank)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = chosen - chosen.se, ymax = chosen+chosen.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1,5), labels = c(1,5)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  theme(legend.position = 'none')

ggplot(graph.choice %>% filter(cond.fac == 'Losses'), aes(x = s2_value_rank, y = chosen, group = s1_value_rank, color = s1_value_rank)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = chosen - chosen.se, ymax = chosen+chosen.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(1,5), labels = c(1,5)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  theme(legend.position = 'none')

# get percentage of people choosing the best option from their consideration set in each condition
graph.choice %>% filter(s2_value_rank == 5) %>% group_by(cond.fac) %>% summarize(chosen.m = mean(chosen))

# stats
m.choice = mlogit(chosen ~ s1_value_rank + s2_value_rank | -1, df.logit, panel = T,
                  rpar = c(s1_value_rank = "n", s2_value_rank = "n"), halton = NA, R = 1000, tol = .001)
summary(m.choice)

m.choice.null = mlogit(chosen ~ s2_value_rank | -1, df.logit, panel = T,
                       rpar = c(s2_value_rank = "n"), halton = NA, R = 1000, tol = .001)

ll1 = logLik(m.choice)
BIC1 = attr(ll1, 'df') * log(length(m.choice$fitted.values)) - 2 * as.numeric(ll1)

ll0 = logLik(m.choice.null)
BIC0 = attr(ll0, 'df') * log(length(m.choice.null$fitted.values)) - 2 * as.numeric(ll0)

BFnull.s1choice = exp((BIC1 - BIC0) / 2)
BFnull.s1choice
# save --------------------------------------------------------------------

save.image(paste0(expt, '/analysis.RData'))

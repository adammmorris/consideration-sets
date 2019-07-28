# setup -------------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(lme4)
require(lmerTest)
require(mlogit)
require(lattice)
require(stringdist)
require(ggstatsplot)
require(plotly)
require(rsm)
require(rje)
require(BayesFactor)
require(ggthemr)

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

# choose which experiment -------------------------------------------------

expt = readline(prompt="Which months experiment do you want to analyze? (Enter 1, 2, or 3): ")

# import data -------------------------------------------------------------

if (expt != '3') {
  cs.flipped = F 
  mid.value = 6
  numTrials = 132
} else {
  cs.flipped = T
  mid.value = 8
  numTrials = 120
}

numWords = 12
minNAs = 1
numQuestions = 2

# Load data
df.demo = read.csv(paste0(expt, '/demo.csv'), stringsAsFactors = F) %>% arrange(subject) %>% mutate(total_time_real = total_time / 60000)
df.words.raw = read.csv(paste0(expt, '/words.csv'), stringsAsFactors = F) %>% arrange(subject, word_ind)
df.s1.raw = read.csv(paste0(expt, '/s1.csv'), stringsAsFactors = F) %>% arrange(subject)
df.s2.raw = read.csv(paste0(expt, '/s2.csv'), stringsAsFactors = F) %>% arrange(subject, question_order)

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
  
  df.s2$s2_value[i] = ifelse(is.na(cind), NA, df.words$s2_value[df.words$subject == subj.name & df.words$word == creal])
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
  
  exclude = df.demo.temp$write_down == 'Yes' || df.s2.subj.temp$comp_check_pass < 1 ||
    df.s2.subj.temp$numNAs > minNAs ||
    df.s2.subj.temp$numTrials != numQuestions ||
    df.s1.subj.temp$numTrials != numTrials ||
    df.s1.temp$pctCorrect_choice < .7 ||
    subj.name == 'bb22886bb143c2600650c8a4c2069c8e' # someone who glitched
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

df.words.filt = df.words.filt %>% mutate(chosen = ifelse(in.cs, 0, NA), chosen_noNA = 0,
                                         high_s1value = factor(s1_value > mid.value, c(F,T), c('Low', 'High')))

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

# create data frame for multinomial logit (for testing what influences choice out of consideration set)

df.logit = df.words.filt %>% filter(in.cs == T) %>% select(subject, s1_value, s2_value, chosen) %>%
  mutate(chosen = as.logical(chosen), intercept = 1) %>%
  rowwise() %>%
  mutate(subject.id = which(as.character(subject) == unique(df.words.filt$subject))) %>%
  group_by(subject.id) %>%
  mutate(option.id = row_number()) %>%
  mlogit.data(choice = "chosen", shape = "long", id.var = "subject.id", alt.var = "option.id", chid.var = "subject.id")

# results ----------------------------------------------------------

if (expt == '3') {
  df.words.all = df.words.filt
  df.words.freq = df.words.filt %>% filter(cond > 0)
  df.words.filt = df.words.filt %>% filter(cond == 0)
}

## effects of stage 1 value

# effect on generation -- dichotomized
graph.s1.binary = df.words.filt %>% group_by(high_s1value, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(high_s1value) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(chosen))

if (expt != '3') {
  cs.s1.binary.breaks = c(.35, .4)
  cs.s1.binary.lims = c(.34, .41)
} else {
  cs.s1.binary.breaks = c(.35, .45)
  cs.s1.binary.lims = c(.35, .47)
}
ggplot(graph.s1.binary, aes(x = high_s1value, y = in.cs.m)) +
  geom_point(size = 5, color = 'black') +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2, color = 'black') +
  xlab('') + ylab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = cs.s1.binary.breaks, limits = cs.s1.binary.lims) + 
  geom_smooth(method='lm')

# effect on generation -- undichotomized
graph.s1.full = df.words.filt %>% group_by(s1_value) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), in.cs.se = sqrt(in.cs * (1-in.cs) / n()),
            chosen = mean(chosen, na.rm = T), chosen.se = sqrt(chosen * (1-chosen) / n()))

ggplot(graph.s1.full, aes(x = s1_value, y = in.cs)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs - in.cs.se, ymax = in.cs+in.cs.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(min(graph.s1.full$s1_value), max(graph.s1.full$s1_value))) +
  scale_y_continuous(breaks = c(.3,.45), limits = c(.29,.52)) +
  betterLine(graph.s1.full, in.cs ~ s1_value)

# stats for generation
m.cs.s1 = glmer(in.cs~s1_value+(s1_value||subject),
               data = df.words.filt,
               family='binomial')
summary(m.cs.s1)

# effect on choice
ggplot(graph.s1.full, aes(x = s1_value, y = chosen)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = chosen - chosen.se, ymax = chosen+chosen.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_continuous(breaks = c(min(graph.s1.full$s1_value), max(graph.s1.full$s1_value))) +
  scale_y_continuous(breaks = c(.15, .3), limits = c(.13, .3)) +
  betterLine(graph.s1.full, chosen ~ s1_value)

# stats for choice

m.choice.s1 = mlogit(chosen ~ s1_value | -1, df.logit, panel = T,
                     rpar = c(s1_value = "n"), correlation = F, halton = NA, R = 1000, tol = .001)
summary(m.choice.s1)

ll1 = logLik(m.choice.s1)
BIC1 = attr(ll1, 'df') * log(length(m.choice.s1$fitted.values)) - 2 * as.numeric(ll1)

m.choice.s1.null = mlogit(chosen ~ intercept | -1, df.logit, correlation = F, halton = NA, R = 1000, tol = .001)
summary(m.choice.s1.null)

ll0 = logLik(m.choice.s1.null)
BIC0 = attr(ll0, 'df') * log(length(m.choice.s1.null$fitted.values)) - 2 * as.numeric(ll0)

BFnull = exp((BIC1 - BIC0) / 2)
BFnull

## effects of stage 2 value

# effect on generation
graph.s2 = df.words.filt %>% group_by(s2_value, subject) %>%
  summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T)) %>%
  group_by(s2_value) %>%
  summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(chosen))

ggplot(graph.s2, aes(x = s2_value, y = in.cs.m)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
  labs(x = '', y = '')+
  scale_y_continuous(breaks = c(.2, .5), limits = c(.18, .57)) +
  scale_x_continuous(breaks = c(2, 25), limits = c(2,25)) +
  betterLine(graph.s2, in.cs.m ~ s2_value)

# stats for generation
m.cs.s2 = glmer(in.cs~s2_value+(s2_value||subject),
                data = df.words.filt,
                family='binomial')
summary(m.cs.s2)

# effect on choice
ggplot(graph.s2, aes(x = s2_value, y = chosen.m)) +
  geom_point(size = 5) + geom_line() +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .2) +
  labs(x = '', y = '')+
  scale_y_continuous(breaks = c(0, .8), limits = c(0, .87)) +
  scale_x_continuous(breaks = c(2, 25), limits = c(2,25)) +
  betterLine(graph.s2, chosen.m ~ s2_value)

# stats for choice
m.choice.s2 = mlogit(chosen ~ s2_value | -1, df.logit, panel = T,
                     rpar = c(s2_value = "n"), correlation = F, halton = NA, R = 1000, tol = .001)
summary(m.choice.s2)

## analyze stuff unique to the different experiments

if (expt != '3') {
  ## in the first & second months experiments, we analyzed for "checkmark" shape
  # how does effect change if you drop the most extreme values?
  m.checkmark1 = glmer(in.cs~s1_value+(s1_value||subject),
                       data = df.words.filt %>% filter(!(s1_value %in% c(1,12))),
                       family='binomial')
  summary(m.checkmark1)
  
  # is there a quadratic component to the effect?
  m.checkmark2 = glmer(in.cs~s1_value+I(s1_value^2)+(s1_value+I(s1_value^2)||subject),
                       data = df.words.filt %>% mutate(s1_value = s1_value - mean(s1_value)),
                       family='binomial')
  summary(m.checkmark2) 
} else {
  ## in the third months experiment, we analyzed for an effect of choice frequency on option generation
  df.words.freq = df.words.freq %>% mutate(cond_fac = factor(cond > 3, levels = c(F, T), labels = c('Low', 'High')))
  graph.freq = df.words.freq %>%
    group_by(cond_fac, subject) %>%
    summarize(in.cs = mean(in.cs, na.rm = T), chosen = mean(chosen, na.rm = T), freq = mean(freq, na.rm = T)) %>%
    group_by(cond_fac) %>%
    summarize(in.cs.m = mean(in.cs, na.rm = T), in.cs.se = se(in.cs),
              chosen.m = mean(chosen, na.rm = T), chosen.se = se(chosen),
              freq.m = mean(freq), freq.se = se(freq))
  ggplot(graph.freq, aes(x = cond_fac, y = in.cs.m)) +
    geom_point(size = 5) + geom_line() +
    geom_smooth(method='lm')+
    geom_errorbar(aes(ymin = in.cs.m - in.cs.se, ymax = in.cs.m+in.cs.se), width = .2) +
    xlab('Choice frequency') + ylab('Prob. of\nword coming to mind') +
    scale_y_continuous(breaks = c(.3,.4), limits = c(.3,.41)) +
    theme(axis.title = element_text(size = 24))
  
  m.freq = glmer(in.cs~cond_fac+(cond_fac||subject),
                 data = df.words.freq,
                 family='binomial')
  summary(m.freq)
  
  m.freq.null = glmer(in.cs~1+(cond_fac||subject),
                 data = df.words.freq,
                 family='binomial')
  
  # bayes factor in favor of null
  BFnull.freq = exp((BIC(m.freq) - BIC(m.freq.null)) / 2)
  BFnull.freq
}

# save --------------------------------------------------------------------

save.image(paste0(expt, '/analysis.RData'))

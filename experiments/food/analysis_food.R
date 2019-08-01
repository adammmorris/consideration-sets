# setup -------------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(ggExtra)
require(tidyr)
require(lme4)
require(lmerTest)
require(mlogit)
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

se = function(x) {return(sd(x, na.rm = T) / sqrt(sum(!is.na(x))))}
dodge <- position_dodge(width=0.9)

betterLine = function(data, formula, color = '#105db0') {
  lg = summary(lm(formula, data))$coefficients
  return(c(geom_abline(intercept = lg[1,1], slope = lg[2,1], color = color, size = 1.25),
           geom_abline(intercept = lg[1,1] - lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8),
           geom_abline(intercept = lg[1,1] + lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8)))
}

# choose which experiment -------------------------------------------------

expt = readline(prompt="Which dinner experiment do you want to analyze? (Enter 1 or 2): ")

# import data -------------------------------------------------------------
df = read.csv(paste0(expt, '/data.csv'), stringsAsFactors = F) %>% filter(DistributionChannel != 'preview', Status == 0) %>%
  dplyr::select(ResponseId, choice_30, choicetime_30_Page.Submit,
                cs_1, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_13,
                cs_time_Page.Submit,
                val_q_1, val_q_2, val_q_3, val_q_4, val_q_5, val_q_6, val_q_7, val_q_8, val_q_9,
                val_time_Page.Submit,
                val_other_q_1, val_other_q_2, val_other_q_3, val_other_q_4, val_other_q_5, val_other_q_6, val_other_q_7, val_other_q_8, val_other_q_9,
                Q159_Page.Submit)


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

df.subj = df2 %>% group_by(subject) %>%
  summarize(cs.size = n(), nas = any(is.na(val) | is.na(val.spec)),
            var.choice = sd(c(val, val.spec)),
            val = mean(val), val.spec = mean(val.spec))

df2.filt = df2 %>% filter(!(subject %in% df.subj$subject[df.subj$nas == T | df.subj$var.choice == 0]))
df2.subj.filt = df.subj %>% filter(!nas & var.choice != 0)


# graphs and analysis -----------------------------------------------------

# test how correlated general & specific value are
cor.test(df2.filt$val, df2.filt$val.spec)

## raw data
p = ggplot(data = df2.filt %>% mutate(chosen = ifelse(chosen == T, T, NA)), aes(x = val, y = val.spec, color = chosen, group = chosen)) +
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

## add info to df

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

# for graph, do marginal probability
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

# for stats, we're going to do it controlling for each other
df.model.cs = df2.filt %>% group_by(subject) %>%
  summarize(ll = mean(val.fac == 'Low' & val.spec.fac == 'Low'),
            lh = mean(val.fac == 'Low' & val.spec.fac == 'High'),
            hl = mean(val.fac == 'High' & val.spec.fac == 'Low'),
            hh = mean(val.fac == 'High' & val.spec.fac == 'High')) %>%
  gather("type", "prob", -subject) %>%
  mutate(val.high = type %in% c('hl', 'hh'), val.spec.high = type %in% c('lh', 'hh'))

m.cs = lmer(prob ~ val.high + val.spec.high + (1 | subject), data = df.model.cs)
summary(m.cs)

## analyze choice from CS

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
summary(m.selection)

# save --------------------------------------------------------------------


save.image(paste0(expt, '/analysis.RData'))

# setup -------------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(ggExtra)
require(lme4)
require(lmerTest)
require(mlogit)
require(lattice)
require(stringdist)
require(ggstatsplot)
require(plotly)
require(rsm)
require(rje)
require(tidyr)

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

df2 = df2 %>% mutate(val.diff = val - val.spec)

df.subj = df2 %>% group_by(subject) %>%
  summarize(cs.size = n(), nas = any(is.na(val) | is.na(val.spec)),
            var.choice = sd(c(val, val.spec)), val.diff = mean(val.diff),
            val = mean(val), val.spec = mean(val.spec))

df2.filt = df2 %>% filter(!(subject %in% df.subj$subject[df.subj$nas == T | df.subj$var.choice == 0]))
df2.subj.filt = df.subj %>% filter(!nas & var.choice != 0)


# graphs and analysis -----------------------------------------------------

cor.test(df2.filt$val, df2.filt$val.spec)

## raw data
p = ggplot(data = df2.filt %>% mutate(chosen = ifelse(chosen == T, T, NA)), aes(x = val, y = val.spec, color = chosen, group = chosen)) +
  geom_jitter(size = 2) +
  geom_vline(xintercept = 4) +
  geom_hline(yintercept = 4) +
  labs(x = '', y = '') +
  guides(color = F)
ggMarginal(p, type='histogram', xparams = list(bins=8), yparams = list(bins = 8))

m1 = lmer(val ~ val.spec + (val.spec | subject), data = df2.filt %>% mutate(val = val - 4, val.spec = val.spec - 4))
summary(m1)

m2 = lmer(val.spec ~ val + (val | subject), data = df2.filt %>% mutate(val = val - 4, val.spec = val.spec - 4))
summary(m2)

## transformation
nvals = length(unique(df2.filt$val))
probs.mf = numeric(nvals)
probs.mb = numeric(nvals)
probs.mf.chosen = numeric(nvals)
probs.mb.chosen = numeric(nvals)
df.figure = data.frame()

cs.size = mean((df2.filt %>% group_by(subject) %>% summarize(cs.size = n()))$cs.size)

for (i in 1:nvals) {
  probs.mf[i] = mean(df2.filt$val == i)
  probs.mb[i] = mean(df2.filt$val.spec == i)
  
  df2.chosen = df2.filt %>% filter(chosen == 1)
  probs.mf.chosen[i] = mean(df2.chosen$val == i)
  probs.mb.chosen[i] = mean(df2.chosen$val.spec == i)
  
  df.figure = rbind(df.figure,
                    data.frame(rank = i,
                               chosen = 0,
                               type = 0,
                               prob = probs.mf[i]))
  df.figure = rbind(df.figure,
                    data.frame(rank = i,
                               chosen = 0,
                               type = 1,
                               prob = probs.mb[i]))
  df.figure = rbind(df.figure,
                    data.frame(rank = i,
                               chosen = 1,
                               type = 0,
                               prob = probs.mf.chosen[i] / probs.mf[i] / cs.size))
  df.figure = rbind(df.figure,
                    data.frame(rank = i,
                               chosen = 1,
                               type = 1,
                               prob = probs.mb.chosen[i] / probs.mb[i] / cs.size))
}

df.figure = df.figure %>%
  mutate(chosen.fac = factor(chosen, c(0,1), c('Considered', 'Chosen')),
         type.fac = factor(type, c(0,1), c('MF', 'MB'))) %>%
  select(-c(chosen, type))
  
ggplot(df.figure %>% spread(chosen.fac, prob), aes(x = rank, y = Considered, color = type.fac, fill = type.fac, group = type.fac)) +
  geom_area(position = 'identity', alpha = .7, color = NA) +
  geom_smooth(method='lm', se = F, size = 1.5) +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = c('#e06958', '#71a4d9')) +
  scale_color_manual(values = c('#d11a02', '#105db0')) +
  scale_y_continuous(breaks = c(0, .3), limits = c(0, .3)) +
  scale_x_continuous(labels = c(1, 7), breaks = c(1,7)) +
  theme(legend.position = 'none')

ggplot(df.figure %>% spread(chosen.fac, prob), aes(x = rank, y = Chosen, color = type.fac, fill = type.fac, group = type.fac)) +
  geom_area(position = 'identity', alpha = .7, color = NA) +
  geom_smooth(method='lm', se = F, size = 1.5) +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = c('#e06958', '#71a4d9')) +
  scale_color_manual(values = c('#d11a02', '#105db0')) +
  scale_y_continuous(breaks = c(0, .5), limits = c(0, .52)) +
  scale_x_continuous(labels = c(1, 7), breaks = c(1,7)) +
  theme(legend.position = 'none')

df.chosen.mf = df2.filt %>% group_by(val, subject) %>%
  summarize(chosen = mean(chosen)) %>%
  group_by(val) %>%
  summarize(chosen.m = mean(chosen), chosen.se = se(chosen))
df.chosen.mb = df2.filt %>% group_by(val.spec, subject) %>%
  summarize(chosen = mean(chosen)) %>%
  group_by(val.spec) %>%
  summarize(chosen.m = mean(chosen), chosen.se = se(chosen))
df.chosen = rbind(df.chosen.mf, df.chosen.mb %>% mutate(val = val.spec) %>% select(-val.spec)) %>%
  mutate(type = rep(c(0,1), each = 7), type.fac = factor(type, c(0,1), c('MF', 'MB')))
ggplot(df.chosen, aes(x = val, y = chosen.m, color = type.fac, fill = type.fac, group = type.fac)) +
  #geom_area(position = 'identity', alpha = .7, color = NA) +
  #geom_smooth(method='lm', se = F, size = 1.5) +
  geom_line() +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = c('#e06958', '#71a4d9')) +
  scale_color_manual(values = c('#d11a02', '#105db0')) +
  #scale_y_continuous(breaks = c(0, .5), limits = c(0, .52)) +
  #scale_x_continuous(labels = c(1, 7), breaks = c(1,7)) +
  theme(legend.position = 'none')

## dichotomized, by-subject transformation

subjs = length(unique(df2.filt$subject))
df.figure.subj = data.frame()
for (subj in 1:subjs) {
  subj.name = df2.filt$subject[subj]
  df2.temp = df2.filt %>% filter(subject == subj.name)
  df2.temp.chosen = df2.temp %>% filter(chosen == 1)
  cs.size = nrow(df2.temp)
  
  df2.temp$val[df2.temp$val == 4] = 4 + sample(c(-.5, .5), 1)
  df2.temp$val.spec[df2.temp$val.spec == 4] = 4 + sample(c(-.5, .5), 1)
  
  cutoff = 4
  p = mean(df2.temp$val <= cutoff)
  df.figure.subj = rbind(df.figure.subj,
                         data.frame(subject = subj,
                                    rank = 0,
                                    type = 0,
                                    considered = p,
                                    chosen = ifelse(p == 0, NA, mean(df2.temp.chosen$val <= cutoff) / p / cs.size)))
  p = mean(df2.temp$val > cutoff)
  df.figure.subj = rbind(df.figure.subj,
                         data.frame(subject = subj,
                                    rank = 1,
                                    type = 0,
                                    considered = p,
                                    chosen = ifelse(p == 0, NA, mean(df2.temp.chosen$val > cutoff) / p / cs.size)))
  p = mean(df2.temp$val.spec <= cutoff)
  df.figure.subj = rbind(df.figure.subj,
                         data.frame(subject = subj,
                                    rank = 0,
                                    type = 1,
                                    considered = p,
                                    chosen = ifelse(p == 0, NA, mean(df2.temp.chosen$val.spec <= cutoff) / p / cs.size)))
  p = mean(df2.temp$val.spec > cutoff)
  df.figure.subj = rbind(df.figure.subj,
                         data.frame(subject = subj,
                                    rank = 1,
                                    type = 1,
                                    considered = p,
                                    chosen = ifelse(p == 0, NA, mean(df2.temp.chosen$val.spec > cutoff) / p / cs.size)))
}

df.figure.subj = df.figure.subj %>%
  mutate(type.fac = factor(type, c(0,1), c('MF', 'MB')),
         rank.fac = factor(rank, c(0,1), c('Low', 'High')))

df.figure.subj.coll = df.figure.subj %>% group_by(type.fac, rank.fac) %>%
  summarize(considered.m = mean(considered, na.rm = T), considered.se = se(considered),
            chosen.m = mean(chosen, na.rm = T), chosen.se = se(chosen))
ggplot(df.figure.subj.coll, aes(x = rank.fac, y = considered.m, group = type.fac, color = type.fac)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = considered.m - considered.se, ymax = considered.m + considered.se), width = .2)
ggplot(df.figure.subj.coll, aes(x = rank.fac, y = chosen.m, group = type.fac, color = type.fac)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m + chosen.se), width = .2)

m.test = lmer(considered ~ rank.fac * type.fac + (1 | subject), data = df.figure.subj)
summary(m.test)



## logit
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

m.selection = mlogit(chosen ~ val + val.spec | -1, df2.logit2, panel = T,
                     rpar = c(val = "n", val.spec = "n"), correlation = F, halton = NA, R = 1000, tol = .001)
summary(m.selection)

ll1 = logLik(m.selection)
BIC1 = attr(ll1, 'df') * log(length(m.selection$fitted.values)) - 2 * as.numeric(ll1)

m.selection.nogeneral = mlogit(chosen ~ val.spec | -1, df2.logit2, panel = T,
                     rpar = c(val.spec = "n"), correlation = F, halton = NA, R = 1000, tol = .001)
summary(m.selection.nogeneral)

ll0 = logLik(m.selection.nogeneral)
BIC0 = attr(ll0, 'df') * log(length(m.selection.nogeneral$fitted.values)) - 2 * as.numeric(ll0)

BFnull = exp((BIC1 - BIC0) / 2)
BFnull

# save --------------------------------------------------------------------


save.image(paste0(expt, '/analysis.RData'))

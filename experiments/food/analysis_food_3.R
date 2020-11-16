# setup -------------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(ggExtra)
require(lme4)
require(lmerTest)
require(mlogit)
require(stringdist)
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

# import data -------------------------------------------------------------

df = read.csv('3/data.csv') %>% filter(DistributionChannel != 'preview', Status == 0) %>%
  dplyr::select(StartDate, ResponseId, movieChoice, choiceTime_Page.Submit,
                cs_1, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_13,
                cs_time_Page.Submit,
                val_q_1, val_q_2, val_q_3, val_q_4, val_q_5, val_q_6, val_q_7, val_q_8,
                val_time_Page.Submit,
                freq_q_1, freq_q_2, freq_q_3, freq_q_4, freq_q_5, freq_q_6, freq_q_7, freq_q_8,
                freq_time_Page.Submit)


df2 = data.frame(subject = character(), food = character(), val = numeric(), freq = numeric())

cs.names = c('cs_1', 'cs_3', 'cs_4', 'cs_5', 'cs_6', 'cs_7', 'cs_13')
val.names = c('val_q_1', 'val_q_2', 'val_q_3', 'val_q_4', 'val_q_5', 'val_q_6', 'val_q_7', 'val_q_8')
freq.names = c('freq_q_1', 'freq_q_2', 'freq_q_3', 'freq_q_4', 'freq_q_5', 'freq_q_6', 'freq_q_7', 'freq_q_8')
for (i in 1:nrow(df)) {
  subj = df$ResponseId[i]
  for (j in 1:length(cs.names)) {
    food = as.character(df[i, cs.names[j]])
    if (nchar(food) > 0) {
      val = as.numeric(as.character(df[i, val.names[j]]))
      freq = as.numeric(as.character(df[i, freq.names[j]]))
      df2 = rbind(df2, data.frame(subject = subj, food = food, val = val, freq = freq))
    }
  }
}


# analyze results ---------------------------------------------------------

# plot data
p = ggplot(data = df2, aes(x = val, y = freq)) +
  geom_jitter(size = 2) +
  geom_vline(xintercept = 4) +
  geom_hline(yintercept = 4) +
  labs(x = '', y = '')
ggMarginal(p, type='histogram', xparams = list(bins=8), yparams = list(bins = 8))

cor.test(df2$val, df2$freq)

# test whether general value > midpoint (4), controlling for freq
m1 = lmer(val ~ freq + (freq | subject), data = df2 %>% mutate(val = val - 4, freq = freq - 4))
summary(m1)

# test whether freq > midpoint (4), controlling for value
m2 = lmer(freq ~ val + (val | subject), data = df2 %>% mutate(val = val - 4, freq = freq - 4))
summary(m2)

m2.null = lmer(freq ~ 0+val + (val | subject), data = df2 %>% mutate(val = val - 4, freq = freq - 4))
BFnull.freq = exp((BIC(m2) - BIC(m2.null)) / 2)
BFnull.freq

## do proper "inferred probability" analysis
df.subj = df2 %>% group_by(subject) %>%
  summarize(cs.size = n(), nas = any(is.na(val) | is.na(freq)),
            var.choice = sd(c(val, freq)),
            val = mean(val), freq = mean(freq))

df2.filt = df2 %>% filter(!(subject %in% df.subj$subject[df.subj$nas == T | df.subj$var.choice == 0]))
df2.subj.filt = df.subj %>% filter(!nas & var.choice != 0)


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
  
  if (df2.filt$freq[i] == 4) {
    df2.filt$freq.four[i] = c2
    c2 = ifelse(c2 == 0, 1, 0)
  } else {
    df2.filt$freq.four[i] = -1
  }
}


df2.filt = df2.filt %>%
  mutate(val.adj = ifelse(val == 4, ifelse(val.four == 0, 3.9, 4.1), val),
         freq.adj = ifelse(val == 4, ifelse(freq.four == 0, 3.9, 4.1), freq),
         val.fac = factor(val.adj > 4, c(F,T), c('Low', 'High')),
         freq.fac = factor(freq.adj > 4, c(F,T), c('Low', 'High')))

df.model.cs = df2.filt %>% group_by(subject) %>%
  summarize(ll = mean(val.fac == 'Low' & freq.fac == 'Low'),
            lh = mean(val.fac == 'Low' & freq.fac == 'High'),
            hl = mean(val.fac == 'High' & freq.fac == 'Low'),
            hh = mean(val.fac == 'High' & freq.fac == 'High')) %>%
  gather("type", "prob", -subject) %>%
  mutate(val.high = type %in% c('hl', 'hh'), freq.high = type %in% c('lh', 'hh'))

m.cs = lmer(prob ~ val.high + freq.high + (1 | subject), data = df.model.cs)
summary(m.cs)

m.cs.null = lmer(prob ~ val.high + (1 | subject), data = df.model.cs)
exp((BIC(m.cs) - BIC(m.cs.null)) / 2)

graph.cs = df2.filt %>% select(subject, val.fac, freq.fac) %>%
  rename(val = val.fac, freq = freq.fac) %>%
  gather("type", "rank", -subject) %>%
  group_by(type, subject) %>%
  summarize(low = mean(rank == 'Low'), high = mean(rank == 'High')) %>%
  gather("rank", "prob", -c(subject,type)) %>%
  group_by(type, rank) %>%
  summarize(prob.m = mean(prob), prob.se = se(prob)) %>%
  group_by() %>%
  mutate(type = factor(type, c('val', 'freq')), rank = factor(rank, c('low', 'high')))

ggplot(graph.cs, aes(x = rank, y = prob.m, group = type, color = type)) +
  geom_point(size = 4) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = prob.m - prob.se, ymax = prob.m + prob.se), width = .2) +
  xlab('') + ylab('') +
  scale_x_discrete(labels = NULL) +
  #scale_y_continuous(limits = c(.15,.85), breaks = c(.2, .8)) +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#d11a02', '#105db0'))

# save --------------------------------------------------------------------

save.image('3/analysis.RData')

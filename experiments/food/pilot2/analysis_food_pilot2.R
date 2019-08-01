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

df = read.csv('data.csv') %>% filter(DistributionChannel != 'preview', Status == 0) %>%
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

# save --------------------------------------------------------------------

save.image('analysis.RData')

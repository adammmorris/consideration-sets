# setup -------------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(ggExtra)
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

se = function(x) {return(sd(x, na.rm = T) / sqrt(length(x)))}
dodge <- position_dodge(width=0.9)

betterLine = function(data, formula, color = '#105db0') {
  lg = summary(lm(formula, data))$coefficients
  return(c(geom_abline(intercept = lg[1,1], slope = lg[2,1], color = color, size = 1.25),
           geom_abline(intercept = lg[1,1] - lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8),
           geom_abline(intercept = lg[1,1] + lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8)))
}

# import data -------------------------------------------------------------

df = read.csv('pilot/data.csv') %>% filter(DistributionChannel != 'preview', Status == 0) %>%
  dplyr::select(ResponseId, choice_5, choice_10, choice_15, choice_20, choice_30,
                choicetime_5_Page.Submit, choicetime_10_Page.Submit, choicetime_15_Page.Submit, choicetime_20_Page.Submit, choicetime_30_Page.Submit,
                cs_1, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_13,
                cs_time_Page.Submit,
                val_q_1, val_q_2, val_q_3, val_q_4, val_q_5, val_q_6, val_q_7, val_q_8,
                val_time_Page.Submit,
                freq_q_1, freq_q_2, freq_q_3, freq_q_4, freq_q_5, freq_q_6, freq_q_7, freq_q_8,
                freq_time_Page.Submit)

df2 = data.frame(subject = character(), food = character(), val = numeric(), freq = numeric(), rt = numeric(), choice = character())

choice.names = c('choice_5', 'choice_10', 'choice_15', 'choice_20', 'choice_30')
cond.names = c('choicetime_5_Page.Submit', 'choicetime_10_Page.Submit', 'choicetime_15_Page.Submit', 'choicetime_20_Page.Submit', 'choicetime_30_Page.Submit')
cs.names = c('cs_1', 'cs_3', 'cs_4', 'cs_5', 'cs_6', 'cs_7', 'cs_8', 'cs_13')
val.names = c('val_q_1', 'val_q_2', 'val_q_3', 'val_q_4', 'val_q_5', 'val_q_6', 'val_q_7', 'val_q_8')
freq.names = c('freq_q_1', 'freq_q_2', 'freq_q_3', 'freq_q_4', 'freq_q_5', 'freq_q_6', 'freq_q_7', 'freq_q_8')
for (i in 1:nrow(df)) {
  subj = df$ResponseId[i]
  choice = NULL
  rt = NULL
  for (j in 1:length(cond.names)) {
    cond = as.character(df[i, cond.names[j]])
    if (nchar(cond) > 0) {
      choice = as.character(df[i, choice.names[j]])
      rt = j
      break
    }
  }
  for (j in 1:length(cs.names)) {
    food = as.character(df[i, cs.names[j]])
    if (nchar(food) > 0) {
      val = as.numeric(as.character(df[i, val.names[j]]))
      freq = as.numeric(as.character(df[i, freq.names[j]]))
      df2 = rbind(df2, data.frame(subject = subj, food = food, val = val, freq = freq, rt = rt, choice = choice))
    }
  }
}

df2$choice = as.character(df2$choice)
times = c(5,10,15,20,25)
df2 = df2 %>% mutate(rt.true = times[rt])


# graphs and analysis -----------------------------------------------------

# test for effect of time condition
df.subj = df2 %>% filter(nchar(choice) > 0) %>%
  group_by(rt.true, subject) %>% summarize(cs.size = n())
df.graph = df.subj %>% group_by(rt.true) %>%
  summarize(cs.size.m = mean(cs.size), cs.size.se = se(cs.size), cs.size.med = median(cs.size))
ggplot(df.graph, aes(x = rt.true, y = cs.size.m)) +
  geom_point(size = 5) +
  betterLine(df.graph, cs.size.m ~ rt.true) +
  geom_errorbar(aes(ymin = cs.size.m - cs.size.se, ymax = cs.size.m + cs.size.se), width = .2) +
  labs(x = '', y = '')

m = lm(cs.size ~ rt.true, data = df.subj)
summary(m)

# look at distribution of choice set sizes
ggplot(df.subj, aes(x=cs.size)) +
  geom_histogram() +
  labs(x = '', y = '') +
  scale_x_continuous(breaks = c(1,4,8), labels = c(1,4,8)) +
  geom_vline(aes(xintercept = mean(df.subj$cs.size)), color = 'red', linetype = 'dashed', size = 1.2)
mean(df.subj$cs.size)
se(df.subj$cs.size)

# save --------------------------------------------------------------------

save.image('pilot/analysis.RData')

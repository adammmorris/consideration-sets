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
             axis.ticks = element_line(size = 3),
             axis.ticks.length = unit(.25, 'cm')
             )
betterLine = function(data, formula, color = '#105db0') {
  lg = summary(lm(formula, data))$coefficients
  return(c(geom_abline(intercept = lg[1,1], slope = lg[2,1], color = color, size = 1.25),
           geom_abline(intercept = lg[1,1] - lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8),
           geom_abline(intercept = lg[1,1] + lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8)))
}

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

se = function(x) {return(sd(x, na.rm = T) / sqrt(length(x)))}
dodge <- position_dodge(width=0.9)



# import data -------------------------------------------------------------

df = read.csv('data.csv', stringsAsFactors = F) %>% filter(DistributionChannel != 'preview', Status == 0) %>%
  dplyr::select(ResponseId, choice_30, choicetime_30_Page.Submit, Q164, Q165_Page.Submit,
                cs_1, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_13,
                cs_time_Page.Submit,
                val_q_1, val_q_2, val_q_3, val_q_4, val_q_5, val_q_6, val_q_7, val_q_8, val_q_9, val_q_10, val_q_11,
                val_time_Page.Submit)


df2 = data.frame(subject = character(), food = character(), chosen = logical(), val = numeric(), choice = character(), cond = numeric(), order = numeric())

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


# graph -------------------------------------------------------------------

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
            val.opp = mean(val.opp), val.opp.se = sqrt(val.opp * (1 - val.opp) / n()),
            val.rank.m = mean(val.rank), val.rank.se = se(val.rank))

ggplot(df2.test, aes(x = cond, y = val.opp,fill = cond)) +
  geom_bar(stat = 'identity', position = dodge) +
  geom_errorbar(aes(ymin = val.opp - val.opp.se, ymax = val.opp + val.opp.se), width = .2, position = dodge, color = 'black') +
  ylab('') +
  xlab('') +
  scale_x_discrete(labels = c('', '')) +
  scale_y_continuous(breaks = c(0,.2,.4)) +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('#8E44AD', '#58D68D'))

# stats
m = glmer(val.opp ~ cond + (1 | subject), family = 'binomial', data = df2.filt)
summary(m)

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

df2.filt %>% filter(chosen == 'Chosen') %>% group_by(cond) %>%
  summarize(best.in.cs = mean(val.rank == 1))


# save --------------------------------------------------------------------

save.image('analysis.RData')

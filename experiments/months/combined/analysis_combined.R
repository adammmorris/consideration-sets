
# setup -------------------------------------------------------------------



require(dplyr)
require(ggplot2)
require(lme4)
require(lmerTest)

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

se = function(x) {return(sd(x, na.rm = T) / sqrt(length(x)))}
dodge <- position_dodge(width=0.9)

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

betterLine = function(data, formula, color = '#105db0') {
  lg = summary(lm(formula, data))$coefficients
  return(c(geom_abline(intercept = lg[1,1], slope = lg[2,1], color = color, size = 1.25),
           geom_abline(intercept = lg[1,1] - lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8),
           geom_abline(intercept = lg[1,1] + lg[1,2], slope = lg[2,1], color = color, linetype = 'dashed', size = 1, alpha = .8)))
}

# load data ----------------------------------------------------------------


load('../1/analysis.rdata')

df1 = df.words.filt %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1value.fac, s2value.fac, in.cs, subject, chosen)

load('../2/analysis.rdata')

df2 = df.words.filt %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1value.fac, s2value.fac, in.cs, subject, chosen)

load('../3/analysis.rdata')

df3 = df.words.filt %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1value.fac, s2value.fac, in.cs, subject, chosen)

df3.all = df.words.all %>%
  mutate(s1value.fac = cut(s1_value, 6, 1:6), s2value.fac = cut(s2_value, 6, 1:6)) %>%
  dplyr::select(s1_value, s2_value, s1value.fac, s2value.fac, in.cs, subject, chosen)

df = rbind(df1,df2,df3) %>% mutate(s1value.num = as.numeric(s1value.fac), s2value.num = as.numeric(s2value.fac))

df.selection = rbind(df1,df2,df3.all) %>% mutate(s1value.num = as.numeric(s1value.fac), s2value.num = as.numeric(s2value.fac))


# make df for logistic
df.logit = df.selection %>% filter(in.cs == T) %>% select(subject, s1value.num, s2_value, chosen) %>%
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

m.s1.cs = glmer(in.cs ~ s1value.num + (s1value.num || subject), data = df, family = 'binomial')
summary(m.s1.cs)

# on choice (for main text)
ggplot(df.s1, aes(x = s1value.num, y = chosen.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0,.8), limits = c(0,.8)) + 
  scale_x_continuous(breaks = c(1,6), labels = c('', ''))
# (for supplement)
ggplot(df.s1, aes(x = s1value.num, y = chosen.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(.15, .3), limits = c(.13, .3)) + 
  scale_x_continuous(breaks = c(1,6), labels = c('1', '6'))+
  betterLine(df.s1, chosen.m ~ s1value.num)

m.s1.choice = mlogit(chosen ~ s1value.num | -1, df.logit, panel = T,
                        rpar = c(s1value.num = "n"), halton = NA, R = 1000, tol = .001)
summary(m.s1.choice)

m.s1.choice.null = mlogit(chosen ~ intercept | -1, df.logit)
summary(m.s1.choice.null)

ll1 = logLik(m.s1.choice)
BIC1 = attr(ll1, 'df') * log(length(m.s1.choice$fitted.values)) - 2 * as.numeric(ll1)

ll0 = logLik(m.s1.choice.null)
BIC0 = attr(ll0, 'df') * log(length(m.s1.choice.null$fitted.values)) - 2 * as.numeric(ll0)

BFnull = exp((BIC1 - BIC0) / 2)
BFnull

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

m.s2.cs = glmer(in.cs ~ s2value.num + (s2value.num || subject), data = df, family = 'binomial')
summary(m.s2.cs)

# on choice
ggplot(df.s2, aes(x = s2_value, y = chosen.m)) +
  geom_point(size = 5, color = 'black') + geom_line(color = 'black') +
  geom_errorbar(aes(ymin = chosen.m - chosen.se, ymax = chosen.m+chosen.se), width = .02, color = 'black') +
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0,.8), limits = c(0,.82)) + 
  scale_x_continuous(breaks = c(2,25), labels = c(2,25))+
  betterLine(df.s2, chosen.m ~ s2_value)

m.s2.choice = mlogit(chosen ~ s2_value | -1, df.logit, panel = T,
                        rpar = c(s2_value = "n"), halton = NA, R = 1000, tol = .001)
summary(m.s2.choice)


# save --------------------------------------------------------------------


save.image('analysis.RData')

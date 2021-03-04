## Generating options and choosing between them rely on distinct forms of value representation
## Adam Morris, Jonathan Phillips, Karen Huang, Fiery Cushman
## Graphs simulation results
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
pkgs = c('dplyr', 'ggplot2')
groundhog.library(pkgs, '2020-03-01')

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


# make graph --------------------------------------------------------------


# which parameters do you want?
n = 20 # number of options (20, 100, 1000, or 10000)
a = 10 # action cost (2, 5, 10, 25, 50, or 100)

path = paste0(n,'_',a)
df = read.csv(paste0('data/', path, '.csv'), header=F)
colnames(df) = c('Output','Correlation','CS.Size', 'Random')
df = df %>% mutate(Correlation = factor(Correlation), Random = factor(Random), CS.Size = factor(CS.Size))

if (n == 20) {
  sizes = c(1:5, 7, 9, 12, 15, 20)
} else {
  sizes = c(1:5, 7, 9, 12, 15, 20, 40, 80, max(as.numeric(as.character(df$CS.Size))))
}
df.filt = df %>% filter(CS.Size %in% sizes)

ggplot(mapping = aes(x = CS.Size, y = Output, color = Correlation, group = Correlation)) +
  stat_summary(data = df.filt %>% filter(Random == 0), fun.y = sum, geom = "line", size = 2) +
  stat_summary(data = df.filt %>% filter(Random == 1), fun.y = sum, geom = "line", linetype = 'dashed', size = 2) +
  xlab('') + ylab('') + theme(legend.position = 'none')

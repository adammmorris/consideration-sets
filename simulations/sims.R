require(dplyr)
require(ggplot2)

# which parameters do you want?
n = 10000
a = 100

path = paste0(n,'_',a)
df = read.csv(paste0('data/', path, '.csv'), header=F)
colnames(df) = c('Output','Correlation','CS.Size', 'Random')
df = df %>% mutate(Correlation = factor(Correlation), Random = factor(Random), CS.Size = factor(CS.Size))#, labels = c('Random sampling', '0.25', '0.50', '0.75', 'Perfect')))

sizes = c(1:5, 7, 9, 12, 15, 20, 40, 80, max(as.numeric(as.character(df$CS.Size))))
df.filt = df %>% filter(CS.Size %in% sizes)

ggplot(mapping = aes(x = CS.Size, y = Output, color = Correlation, group = Correlation)) +
  stat_summary(data = df.filt %>% filter(Random == 0), fun.y = sum, geom = "line", size = 2) +
  stat_summary(data = df.filt %>% filter(Random == 1), fun.y = sum, geom = "line", linetype = 'dashed', size = 2) +
  xlab('') + ylab('') + theme(legend.position = 'none')
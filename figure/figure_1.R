library(gplots)
library(ggplot2)
source('multiplot.R')

cov <- read.table("cov.tsv", stringsAsFactors = F, header=T)
cov$sample <- gsub(pattern = ".dedup.realign.base.txt.gz", replacement = "", x = cov$sample)

p1 <- ggplot(cov, aes(x=sample, y=X20x)) + geom_bar(stat="identity") +
  geom_abline(intercept=80, colour="red", slope=0) +
  ggtitle('20X coverage') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample") +
  ylab("Percentage") +
  scale_y_continuous(limit = c(0, 100)) +
  coord_cartesian(ylim = c(50, 100)) +
  theme(axis.text=element_text(size=18),
       axis.title=element_text(size=24,face="bold"),
       plot.title=element_text(size=30),
       axis.ticks = element_blank(),
       axis.text.x = element_blank())

p2 <- ggplot(cov, aes(x=sample, y=X30x)) + geom_bar(stat="identity") +
  geom_abline(intercept=60, colour="red", slope=0) +
  ggtitle('30X coverage') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample") +
  ylab("Percentage") +
  scale_y_continuous(limit = c(0, 100)) +
  coord_cartesian(ylim = c(50, 100)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"),
        plot.title=element_text(size=30),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

multiplot(p1, p2, cols=2)


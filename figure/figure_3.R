library(ggplot2)
source('multiplot.R')
options(scipen=10000)

class <- c('intronic',
           'UTR3',
           'exonic',
           'ncRNA_exonic',
           'UTR5',
           'ncRNA_intronic',
           'upstream',
           'downstream',
           'intergenic',
           'splicing',
           'ncRNA_splicing')

class <- factor(x = class, levels = rev(class))
df <- data.frame(class = class, count = count)

count <- c(142280, #intronic
           90117, #UTR3
           87226, #exonic
           18049, #ncRNA_exonic
           17311, #UTR5
           12900, #ncRNA_intronic
           6376, #upstream
           5917, #downstream
           1366, #intergenic
           607, #splicing
           95) #ncRNA_splicing

class2 <- c('nonsynonymous SNV',
            'synonymous SNV',
            'unknown',
            'nonframeshift deletion',
            'frameshift deletion',
            'stopgain',
            'nonframeshift insertion',
            'frameshift insertion',
            'stoploss')

class2 <- factor(x = class2, levels = rev(class2))
count2 <- c(43756, 38910, 1303, 926, 708, 565, 542, 495, 49)
df2 <- data.frame(class = class2, count = count2)

p1 <- ggplot(df, aes(x=class, y=count)) + geom_bar(stat="identity") +
  ggtitle("Genomic annotation") +
  ylim(0,160000) +
  coord_flip() +
  geom_text(aes(label=count), size=9, hjust=-0.2) +
  theme(axis.text=element_text(size=28),
        plot.title=element_text(size=34),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p2 <- ggplot(df2, aes(x=class, y=count)) + geom_bar(stat="identity") +
  ggtitle("Exonic annotation") +
  ylim(0,50000) +
  coord_flip() +
  geom_text(aes(label=count), hjust=-0.2, size=10) +
  theme(axis.text=element_text(size=28),
        plot.title=element_text(size=34),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

multiplot(p1, p2, cols = 2)


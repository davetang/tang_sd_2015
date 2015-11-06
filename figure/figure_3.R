library(ggplot2)

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

class <- factor(x = class, levels = rev(class))
df <- data.frame(class = class, Count = count)

p1 <- ggplot(df, aes(x=class, y=count)) + geom_bar(stat="identity") +
ggtitle("Genomic annotation") +
ylim(0,175000) +  theme_bw() +
coord_flip() + geom_text(aes(label=count), hjust=-0.2) + labs(y = "Count") +
theme(title =  element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),axis.ticks=element_blank())

ggsave("Genomic.pdf",p1, width = 16, height = 12, units = c("cm"), dpi = 300)

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

p2 <- ggplot(df2, aes(x=class, y=count)) + geom_bar(stat="identity") +
ylim(0,50000) +  theme_bw() +
coord_flip() + geom_text(aes(label=count), hjust=-0.2) + labs(y = "Count") +
theme(title =  element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),axis.ticks=element_blank())

ggsave("Exonic.pdf",p2, width = 16, height = 12, units = c("cm"), dpi = 300)

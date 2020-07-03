library(survival)
data(pbc)
?pbc

View(pbc)
summary(pbc)

install.packages("tidyverse")
library(tidyverse)

pbc_select <- select(pbc, sex, stage, age)
View(pbc_select)
pbc_select <- select(pbc, -sex, -stage, -age)
View(pbc_select)

pbc_filter <- filter(pbc, sex=="m") 
View(pbc_filter)
pbc_filter <- filter(pbc, age > 70)
View(pbc_filter)

pbc_mutate <- mutate(pbc, new_val=ast/1000) 
View(pbc_mutate)

write.table(pbc_mutate,"pbc_mutate.txt",row.names=F,sep="\t")


ave_age <-summarise(pbc, mean = mean(age))
View(ave_age)

group_by_sex <- group_by(pbc,sex)
View(group_by_sex)
groups(group_by_sex)
groups(pbc)

pbc_final <- pbc %>% group_by(sex) %>% summarise(mean = mean(age))
View(pbc_final)



ggplot(pbc, aes(x=sex)) + geom_bar()

ggplot(pbc, aes(x=age)) + geom_histogram()
ggplot(pbc, aes(x=age)) + geom_histogram(binwidth=5)

ggplot(pbc, aes(x=age)) + geom_density()

ggplot(pbc, aes(x= sex, y=age)) + geom_boxplot()

ggplot(pbc, aes(x = age, y = platelet)) + geom_point()

ggplot(pbc, aes(x = age, y = platelet)) + geom_point() + geom_smooth(method=lm)

ggplot(pbc, aes(x=age, color=sex)) + geom_density()

ggplot(pbc, aes(x=age, y=platelet, color=sex)) + geom_point()
ggplot(pbc, aes(x=age, y=platelet, size=copper)) + geom_point()
ggplot(pbc, aes(x=age, y=platelet, color=sex, size=copper, shape=status)) + geom_point()
ggplot(pbc, aes(x=age, y=platelet, color=sex, size=copper, shape=as.factor(status))) + geom_point()

ggplot(pbc, aes(x=age, y=platelet, color=sex)) + geom_point() + geom_smooth(method=lm)
ggplot(pbc, aes(x=age, y=platelet)) + geom_point() + geom_smooth(method = lm) + facet_wrap(vars(sex), nrow=1)

pbc_long <- pbc %>% select(1,6,11:18) %>% gather(key, value, -sex, -id)
View(pbc_long)
ggplot(pbc_long , aes(value)) + geom_density() + facet_wrap(vars(key), scales = "free")
ggplot(pbc_long , aes(value, color=sex)) + geom_density() +  facet_wrap(vars(key), scales = "free")





exp <- read.table("AMvsEM_deseq2_results.tabular") 
View(exp)

colnames(exp) <- c("GeneID",	"Basemean",	"log2FC",	"StdErr", "WaldStats", "Pvalue", "Padj")
View(exp)

ggplot(data=exp, aes(x=log2FC, y=-log10(Padj), color=log10(Basemean)))  + geom_point(shape=1, size=1.5) + scale_color_gradient(low="green",high="red")

ggsave("volcano.png", width = 6, height = 4)



ADAR2KO <-read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70588/suppl/GSE70588_ADAR2KO_DESEQ_coutfilter.csv.gz")
View(ADAR2KO)
ggplot(data=ADAR2KO, aes(x=log2FoldChange, y=-log10(padj), color=log10(baseMean)))  + geom_point(shape=1, size=1.5) + scale_color_gradient(low="green",high="red")

install.packages("esquisse")
library(esquisse)





install.packages("ggThemeAssist") 
library(ggThemeAssist) 
ggplot(pbc, aes(x= sex, y=age)) + geom_boxplot()



install.packages("plotly") 
library(plotly) 
plot_ly(data = exp, x = ~log2FC, y = ~-log10(Padj), color=~log10(Basemean), type = 'scatter', mode='markers', text= ~GeneID)


install.packages("pheatmap")
library(pheatmap)

counts <- read.table("AMvsEM_deseq2_counts.tabular", sep="\t", header=TRUE)
View(counts)

gene_select <- as.data.frame(filter(counts, X %in% c("75426","71951","53321", "66425")))
View(gene_select)
rownames(gene_select) <- gene_select$X
pheatmap(gene_select[2:8])

library(tidyverse)
library(ggplot2)
library(readxl)
library(pheatmap)
library(ComplexHeatmap)
library(rjson)

#Import of T-Rip data
genes <- read_excel("Enriched transcripts.xlsx", skip = 1)
#Getting the Enriched genes for mock
mock_genes <- genes %>% filter(`Log ratio mock`>1) %>% pull(`geneID`)

#RPM values for the heatmaps
nRPMs <- read_excel("pnas.1712312115.sd01.xlsx", sheet=1, skip =4)
#Alpha values for the scatterplots of UTR vs Half Life
results_data <- read_excel("pnas.1712312115.sd01.xlsx", sheet=2, skip =4)

mock_nRPMs <- nRPMs %>% filter(`geneID` %in% mock_genes)
mock_alphas <- results_data %>% filter(`...1` %in% mock_genes)


### Heatmap code ###

#initial Filtering
a <- data.frame(drop_na(mock_nRPMs)[,c(1,13,14,16)])
rownames(a) <- a$geneID
a <- a[,c(2,3,4)]
a <- a%>% filter_all(all_vars(. < 1))
a <- a[order(rowSums(a), decreasing = TRUE),]

#Ploting the Heatmap
pheatmap::pheatmap(a, 
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   labels_col = c("30 min","60 min","240 min"),
                   legend_breaks = c(0.99,0.5,0),
                   legend_labels = c("1","0.5","0"),
                   border_color = NA,
                   show_rownames = F,
                   main="PB-depleted Mock")



### Scatterplot Figures ###

#Import of a dictionary of PB-enriched primary transcript sequences. Keys are the gene names and values are the sequences
#Could be used twice to import the same dictionary for mock and heat conditions
UTR_5 <- as.data.frame(fromJSON(file = 'relevant_genes_5_mock.json'))
UTR_3 <- as.data.frame(fromJSON(file = 'relevant_genes_3_mock.json'))


conts <- data.frame('geneID'=strtrim(names(UTR_5),9), 
                    'A cont'=unlist(lapply(UTR_5, function (x) (str_count(x, 'A')/nchar(x)))),
                    'T cont'=unlist(lapply(UTR_5, function (x) (str_count(x, 'T')/nchar(x)))),
                    'C cont'=unlist(lapply(UTR_5, function (x) (str_count(x, 'C')/nchar(x)))),
                    'G cont'=unlist(lapply(UTR_5, function (x) (str_count(x, 'G')/nchar(x)))))

conts_3 <- data.frame('geneID'=strtrim(names(UTR_3),9), 
                      'A cont'=unlist(lapply(UTR_3, function (x) (str_count(x, 'A')/nchar(x)))),
                      'T cont'=unlist(lapply(UTR_3, function (x) (str_count(x, 'T')/nchar(x)))),
                      'C cont'=unlist(lapply(UTR_3, function (x) (str_count(x, 'C')/nchar(x)))),
                      'G cont'=unlist(lapply(UTR_3, function (x) (str_count(x, 'G')/nchar(x)))))




t <- mock_alphas[mock_alphas$`...1` %in% conts$geneID, ]%>%arrange(`...1`)%>%pull(alpha_WT, `...1`)
conts <- conts[conts$geneID%in%names(t),]
conts$alpha <- log(2)/unlist(t)
conts<- conts[is.finite(conts$alpha),]

t <- mock_alphas[mock_alphas$`...1` %in% conts_3$geneID, ]%>%arrange(`...1`)%>%pull(alpha_WT, `...1`)
conts_3 <- conts_3[conts_3$geneID%in%names(t),]
conts_3$alpha <- log(2)/unlist(t)
conts_3<- conts_3[is.finite(conts_3$alpha),]


#5 UTR
A <- ggplot(conts,aes(log10(conts$alpha), conts$A.cont))+geom_point(colour="#58a773")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("A content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()

L <- ggplot(conts,aes(log10(conts$alpha), conts$T.cont))+geom_point(colour="#cd3332")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("T content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()

G <- ggplot(conts,aes(log10(conts$alpha), conts$G.cont))+geom_point(colour="#e5c41a")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("G content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()

C <- ggplot(conts,aes(log10(conts$alpha), conts$C.cont))+geom_point(colour="#1ba8e4")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("C content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()


grid.arrange(A,L,C,G, ncol=4, widths=c(30,30,30,30),
             bottom=textGrob('WT Half Life (log10)',
                             gp = gpar(fontsize = 12)),
             top=textGrob("heat 5' UTR",
                          gp = gpar(fontsize = 12)))


#3 UTR
A <- ggplot(conts_3,aes(log10(conts_3$alpha), conts_3$A.cont))+geom_point(colour="#58a773")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("A content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()

L <- ggplot(conts_3,aes(log10(conts_3$alpha), conts_3$T.cont))+geom_point(colour="#cd3332")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("T content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()

G <- ggplot(conts_3,aes(log10(conts_3$alpha), conts_3$G.cont))+geom_point(colour="#e5c41a")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("G content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()

C <- ggplot(conts_3,aes(log10(conts_3$alpha), conts_3$C.cont))+geom_point(colour="#1ba8e4")+
  scale_x_log10()+
  xlab(NULL)+
  ylab("C content (%)")+
  ylim(c(0,1))+
  stat_smooth(method = "loess", formula = 'y~x')+
  theme_bw()


grid.arrange(A,L,C,G, ncol=4, widths=c(30,30,30,30),
             bottom=textGrob('WT Half Life (log10)',
                             gp = gpar(fontsize = 12)),
             top=textGrob("heat 3' UTR",
                          gp = gpar(fontsize = 12)))

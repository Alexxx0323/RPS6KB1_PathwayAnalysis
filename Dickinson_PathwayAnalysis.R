#source("install_and_load_packages.R")

library(plyr)
library(dplyr)
library(readr)
library(ggplot2)
library(reshape)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)
library(hexbin)
library(vsn)

setwd("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL")
source("support_functions.R")

#==================================================================================================================================================================
#==================================================================================================================================================================
#==================================================================================================================================================================
print("### 1: Data importing and manipulation in r")
setwd("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/GSE107934_RAW")

fileList = list.files(pattern = ".*.txt.gz") 

data_all = list.files(pattern = ".*.txt.gz") %>%
  lapply(read.table) %>%
  bind_cols

count_df = data_all[,c(seq(2,58,2))]
count_df

row.names(count_df) = data_all$V1...1
count_df

dataSet_names = c("basal.P1", "AT.1h.P1", "AT.4h.P1", "RT.1h.P1", "RT.4h.P1",
                  "basal.P2", "AT.1h.P2", "AT.4h.P2", "RT.1h.P2", "RT.4h.P2",
                  "basal.P3", "AT.1h.P3", "AT.4h.P3", "RT.1h.P3", "RT.4h.P3",
                  "basal.P4", "AT.1h.P4",             "RT.1h.P4", "RT.4h.P4",
                  "basal.P5", "AT.1h.P5", "AT.4h.P5", "RT.1h.P5", "RT.4h.P5",
                  "basal.P7", "AT.1h.P7", "AT.4h.P7", "RT.1h.P7", "RT.4h.P7")

names(count_df) = dataSet_names
head(count_df)

targets <- as.data.frame(matrix(NA,length(names(count_df)),2))

names(targets) <- c("sample","condition")

targets$sample <- names(count_df)

targets$condition <- gsub(".P[1-7]$","",targets$sample)

head(targets)

count_df <- count_df[rowSums(count_df) > 0,]

count_df[count_df == 0] <- NA

plots <- magicPlotMakerLight(df = log2(count_df), targets = targets)
plot(plots[[1]] + geom_hline(yintercept=0.5)) #violin plot gemo_hline() to see where bimodal starts
plot(plots[[2]]) 
count_df[log2(count_df) < 0.5 ] <- NA
count_df <- count_df[rowSums(is.na(count_df[,c(1:3)])) < 2,] #c(1:3): basal.P1 AT.1h.P1 AT.4h.P1
count_df <- count_df[rowSums(is.na(count_df[,c(4:6)])) < 2,] #c(4:6): RT.1h.P1 RT.4h.P1 basal.P2

fit <- vsn::vsnMatrix(as.matrix(count_df)) #train vsn parameters
meanSdPlot(fit)

count_df_vsn <- as.data.frame(vsn::predict(fit, as.matrix(count_df)))

plots <- magicPlotMakerLight(df = log2(count_df_vsn), targets = targets)
plot(plots[[1]]) #violins
plot(plots[[2]]) #PCA

gene_id_mapping_from_uniprot <- as.data.frame(
  read_delim("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/support/gene_id_mapping_from_uniprot.tab", 
             "\t", escape_double = FALSE, trim_ws = TRUE))
gene_id_mapping_from_uniprot <- gene_id_mapping_from_uniprot[!is.na(gene_id_mapping_from_uniprot$`Gene names`),]

ensembl_to_symbol <- gsub(" .*","",gene_id_mapping_from_uniprot$`Gene names`)

names(ensembl_to_symbol) <- gene_id_mapping_from_uniprot[,1]

row.names(count_df_vsn) <- gsub("[.][1-7]*","",row.names(count_df_vsn))

count_df_vsn <- count_df_vsn[row.names(count_df_vsn) %in% names(ensembl_to_symbol),]

for(i in 1:length(count_df_vsn[,1]))
{
  row.names(count_df_vsn)[i] <- ensembl_to_symbol[row.names(count_df_vsn)[i]]
}

to_write <- count_df_vsn
to_write$gene <- row.names(to_write)

to_write <- to_write[,c(length(to_write[1,]),1:(length(to_write[1,])-1))]
write_csv(to_write, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/count_df_vsn.csv")
write_csv(targets, "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/support/targets.csv")






#==================================================================================================================================================================
#==================================================================================================================================================================
#==================================================================================================================================================================
print("### 2: Differential analysis")
library(limma)
unique(targets$condition)
comparisons = list("B_vs_AT.1h"=c(2,-1), "B_vs_AT.4h"=c(3,-1), "AT.1h_vs_AT.4h"=c(3,-2), 
                   "B_vs_RT.1h"=c(4,-1), "B_vs_RT.4h"=c(5,-1), "RT.1h_vs_RT.4h"=c(5,-4))  

limmaRes = runLimma(measurements = count_df_vsn,
                    targets = targets,
                    comparisons = comparisons)

ttop_B_vs_AT.1h <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

null_model <- pnorm(rnorm(length(ttop_B_vs_AT.1h[,1])))
plot(sort(null_model), sort(ttop_B_vs_AT.1h$P.Value) ,
     xlim = c(1,0), ylim = c(1,0),
     title("AT 1 hour versus basal")) #not bad, not great, let's proceed
abline(coef = c(0,1))

ttop_B_vs_AT.4h <- ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

null_model <- pnorm(rnorm(length(ttop_B_vs_AT.4h[,1])))
plot(sort(null_model), sort(ttop_B_vs_AT.4h$P.Value), 
     xlim = c(1,0), ylim = c(1,0), 
     title("AT 4 hour versus basal"))
abline(coef = c(0,1))

#AT 1 hour versus AT 4 hour
ttop_AT.1h_vs_AT.4h <- ttopFormatter(topTable(limmaRes[[1]], coef = 3, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

null_model <- pnorm(rnorm(length(ttop_B_vs_AT.4h[,1])))
plot(sort(null_model), sort(ttop_B_vs_AT.4h$P.Value), 
     xlim = c(1,0), ylim = c(1,0),
     title("AT 4 hour versus AT 1 hour"))
abline(coef = c(0,1))

#Basal versus RT 1 hour
ttop_B_vs_RT.1h <- ttopFormatter(topTable(limmaRes[[1]], coef = 4, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

null_model <- pnorm(rnorm(length(ttop_B_vs_RT.1h[,1])))
plot(sort(null_model), sort(ttop_B_vs_RT.1h$P.Value), 
     xlim = c(1,0), ylim = c(1,0),
     title("RT 1 hour versus basal")) 
abline(coef = c(0,1))

#Basal versus RT 4 hour
ttop_B_vs_RT.4h <- ttopFormatter(topTable(limmaRes[[1]], coef = 5, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

null_model <- pnorm(rnorm(length(ttop_B_vs_RT.4h[,1])))
plot(sort(null_model), sort(ttop_B_vs_RT.4h$P.Value), 
     xlim = c(1,0), ylim = c(1,0),
     title("RT 4 hour versus basal"))
abline(coef = c(0,1))


#RT 1 hour versus RT 4 hour
ttop_RT.1h_vs_RT.4h <- ttopFormatter(topTable(limmaRes[[1]], coef = 6, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

null_model <- pnorm(rnorm(length(ttop_RT.1h_vs_RT.4h[,1])))
plot(sort(null_model), sort(ttop_RT.1h_vs_RT.4h$P.Value),
     xlim = c(1,0), ylim = c(1,0),
     title("RT 4 hour versus RT 1 hour")) 
abline(coef = c(0,1))

write_csv(ttop_B_vs_AT.1h, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_B_vs_AT.1h.csv") #no diff
write_csv(ttop_B_vs_AT.4h, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_B_vs_AT.4h.csv")
write_csv(ttop_AT.1h_vs_AT.4h, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_AT.1h_vs_AT.4h.csv")
write_csv(ttop_B_vs_RT.1h, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_B_vs_RT.1h.csv") #no diff
write_csv(ttop_B_vs_RT.4h, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_B_vs_RT.4h.csv") #high
write_csv(ttop_RT.1h_vs_RT.4h, 
          file = "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_RT.1h_vs_RT.4h.csv")






#==================================================================================================================================================================
#==================================================================================================================================================================
#==================================================================================================================================================================
print("### 3: Pathway activity (PROGENy)")
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)

Normalised_counts = count_df_vsn
Experimental_design = targets

Normalised_counts_matrix <- Normalised_counts %>% 
    dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
    # tibble::column_to_rownames(var = "gene") %>% 
    as.matrix()

#Comparison of each condition
ttop_B.vs.AT.1h_matrix <- ttop_B_vs_AT.1h %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    dplyr::select(-"ID") %>%
    as.matrix()

ttop_B_vs_AT.4h_matrix <- ttop_B_vs_AT.4h %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    dplyr::select(-"ID") %>%
    as.matrix()

ttop_AT.1h_vs_AT.4h_matrix <- ttop_AT.1h_vs_AT.4h %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    dplyr::select(-"ID") %>%
    as.matrix()

ttop_B_vs_RT.1h_matrix <- ttop_B_vs_RT.1h %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    dplyr::select(-"ID") %>%
    as.matrix()

ttop_B_vs_RT.4h_matrix <- ttop_B_vs_RT.4h %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    dplyr::select(-"ID") %>%
    as.matrix()

ttop_RT.1h_vs_RT.4h_matrix <- ttop_RT.1h_vs_RT.4h %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    dplyr::select(-"ID") %>%
    as.matrix()




#==================================================================================================================================================================
#==================================================================================================================================================================
#==================================================================================================================================================================
print("## Pathway activity with Progeny (pathway activity estimator)")
PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)

paletteLength <- 100
myColor <- colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
    length.out=ceiling(paletteLength/2) + 1),
    seq(max(Activity_counts)/paletteLength, 
    max(Activity_counts), 
    length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
    fontsize_row = 10, fontsize_col = 10, 
    color=myColor, breaks = progenyBreaks, 
    main = "PROGENy (100)", angle_col = 45,
    treeheight_col = 0,  border_color = NA)

PathwayActivity_zscore <- progeny(ttop_B.vs.AT.1h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
    t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>% #add column name to all rows
    dplyr::arrange(NES) %>% #sort pathway by ascending NES
    dplyr::mutate(Pathway = factor(Pathway))

#reorder the levels of Pathway variables based on the magnitude of its second arg NES
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    ggtitle("AT 1 hour versus basal") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Pathways")

# getting a model matrix with 100 top significant genes for each pathway and converting to df
# https://saezlab.github.io/progeny/reference/index.html
prog_matrix <- getModel("Human", top=100) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("GeneID")

ttop_B_vs_AT.1h_df <- ttop_B.vs.AT.1h_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneID")

scat_plots <- progeny::progenyScatter(df = ttop_B_vs_AT.1h_df, #y, gene level stat t values
    weight_matrix = prog_matrix, #x
    statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`MAPK`) # red/blue = positive/negative contributions of genes to the progeny pathway score

PathwayActivity_CARNIVALinput_B_vs_AT.1h <- progeny(ttop_B.vs.AT.1h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
    t () %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput_B_vs_AT.1h)[2] <- "score"
write_csv(PathwayActivity_CARNIVALinput_B_vs_AT.1h, 
    "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/PathwayActivity_CARNIVALinput_B_vs_AT.1h.csv")


## AT 4 hour vs basal
PathwayActivity_zscore <- progeny(ttop_B_vs_AT.4h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
    t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    ggtitle("AT 4 hour versus basal") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Pathways")

prog_matrix <- getModel("Human", top=100) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("GeneID")

ttop_B_vs_AT.4h_df <- ttop_B_vs_AT.4h_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneID")

scat_plots <- progeny::progenyScatter(df = ttop_B_vs_AT.4h_df, 
    weight_matrix = prog_matrix, 
    statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`NFkB`) 

PathwayActivity_CARNIVALinput_B_vs_AT.4h <- progeny(ttop_B_vs_AT.4h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
    t () %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput_B_vs_AT.4h)[2] <- "score"
write_csv(PathwayActivity_CARNIVALinput_B_vs_AT.4h, 
    "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/PathwayActivity_CARNIVALinput_B_vs_AT.4h.csv")


## AT 1 hour versus AT 4 hour

PathwayActivity_zscore <- progeny(ttop_AT.1h_vs_AT.4h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
    t()
colnames(PathwayActivity_zscore) <- "NES"


PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    ggtitle("AT 1 hour versus AT 4 hour") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Pathways")


print("## RT 1 hour versus basal")
PathwayActivity_zscore <- progeny(ttop_B_vs_RT.1h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
    t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    ggtitle("Basal versus RT 1 hour") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Pathways")

print("## RT 4 hour versus basal")

PathwayActivity_zscore <- progeny(ttop_B_vs_RT.4h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
    t()
colnames(PathwayActivity_zscore) <- "NES"


PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    ggtitle("RT 4 hour versus basal") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Pathways")

prog_matrix <- getModel("Human", top=100) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("GeneID")

ttop_B_vs_RT.4h_df <- ttop_B_vs_RT.4h_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneID")

scat_plots <- progeny::progenyScatter(df = ttop_B_vs_RT.4h_df, 
    weight_matrix = prog_matrix, 
    statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`NFkB`) 
plot(scat_plots[[1]]$`TNFa`)
plot(scat_plots[[1]]$`MAPK`) 
plot(scat_plots[[1]]$`EGFR`) 
plot(scat_plots[[1]]$`JAK-STAT`) 


PathwayActivity_CARNIVALinput_B_vs_RT.4h <- progeny(ttop_B_vs_RT.4h_matrix, 
    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
    t () %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput_B_vs_RT.4h)[2] <- "score"
write_csv(PathwayActivity_CARNIVALinput_B_vs_RT.4h, 
    "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/PathwayActivity_CARNIVALinput_B_vs_RT.4h.csv")

print("## RT 1 hour versus RT 4 hour")
PathwayActivity_zscore <- progeny(ttop_RT.1h_vs_RT.4h_matrix, 
    scale=TRUE, organism="Human", top = 200, perm = 10000, z_scores = TRUE) %>%
    t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    ggtitle("RT 1 hour versus RT 4 hour") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Pathways")




#==================================================================================================================================================================
#==================================================================================================================================================================
#==================================================================================================================================================================
print("### 4: Transcription Factor activity (DoRothEA)")
# Normalised_counts_dorothea = count_df_vsn
# Experimental_design_dorothea = targets
Normalised_counts_dorothea <- read_csv("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/count_df_vsn.csv")
Experimental_design_dorothea <- read_csv("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/support/targets.csv")
# ttop_B_vs_AT.4h_dorothea = ttop_B_vs_AT.4h
# ttop_B_vs_RT.4h_dorothea = ttop_B_vs_RT.4h
ttop_B_vs_AT.4h_dorothea <- read_csv("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_B_vs_AT.4h.csv")
ttop_B_vs_RT.4h_dorothea <- read_csv("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/data/ttop_B_vs_RT.4h.csv")

#replace NA entries with 0
Normalised_counts_matrix_dorothea <- Normalised_counts_dorothea %>% 
    dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
    tibble::column_to_rownames(var = "gene") %>% 
    as.matrix()

ttop_B_vs_AT.4h_dorothea_matrix <- ttop_B_vs_AT.4h_dorothea %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    column_to_rownames(var = "ID") %>%
    as.matrix()

ttop_B_vs_RT.4h_dorothea_matrix <- ttop_B_vs_RT.4h_dorothea %>% 
    dplyr::select(ID, t) %>% 
    dplyr::filter(!is.na(t)) %>% 
    column_to_rownames(var = "ID") %>%
    as.matrix()


print("## Transcription Factor activity with DoRothEA")
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat_B_vs_AT.4h <- dorothea::run_viper(ttop_B_vs_AT.4h_dorothea_matrix, regulons,
    options =  list(minsize = 5, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_B_vs_RT.4h <- dorothea::run_viper(ttop_B_vs_RT.4h_dorothea_matrix, regulons,
    options =  list(minsize = 5, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, nes = TRUE))


print("## AT 4 hour versus basal")
tf_activities_stat_top25_B_vs_AT.4h <- tf_activities_stat_B_vs_AT.4h %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "t") %>%
    dplyr::top_n(25, wt = abs(NES)) %>% #top 25, ordered by absolute values of NES, select both - && +ve
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25_B_vs_AT.4h, aes(x = reorder(GeneID, NES), y = NES)) + 
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Transcription Factors")

targets_NFKB1 <- regulons$target[regulons$tf == "NFKB1"] #targets of "TCF4"

volcano_nice(as.data.frame(ttop_B_vs_AT.4h_dorothea[ttop_B_vs_AT.4h_dorothea$ID %in% targets_NFKB1,]), 
    FCIndex = 2, pValIndex = 5, IDIndex = 1,nlabels = 20, label = TRUE, 
    straight = FALSE) 

tf_activities_CARNIVALinput_AT.4h_vs_b<- tf_activities_stat_B_vs_AT.4h %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "TF") 
write_csv(tf_activities_CARNIVALinput_AT.4h_vs_b, "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/TFActivity_CARNIVALinput_AT4h_vs_B.csv")


print("## RT 4 hour versus basal")
tf_activities_stat_top25_B_vs_RT.4h <- tf_activities_stat_B_vs_RT.4h %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "t") %>%
    dplyr::top_n(25, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25_B_vs_RT.4h, aes(x = reorder(GeneID, NES), y = NES)) + 
    ggtitle("RT 4 hour versus basal: TF activity") +
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
    xlab("Transcription Factors")


targets_STAT3 <- regulons$target[regulons$tf == "STAT3"]
volcano_nice(as.data.frame(ttop_B_vs_RT.4h_dorothea[ttop_B_vs_RT.4h_dorothea$ID %in% targets_STAT3,]), 
    FCIndex = 2, pValIndex = 5, IDIndex = 1,nlabels = 20, label = TRUE, 
    straight = FALSE) 


tf_activities_CARNIVALinput_RT.4h_vs_b<- tf_activities_stat_B_vs_RT.4h %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "TF") 
write_csv(tf_activities_CARNIVALinput_RT.4h_vs_b, "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/TFActivity_CARNIVALinput_RT4h_vs_B.csv")

tf_activities_counts <- 
    dorothea::run_viper(Normalised_counts_matrix_dorothea, regulons,
    options =  list(minsize = 5, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::filter(GeneID %in% tf_activities_stat_top25_B_vs_AT.4h$GeneID) %>%
    column_to_rownames(var = "GeneID") %>%
    as.matrix()

tf_activities_vector <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
    colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
    length.out=ceiling(paletteLength/2) + 1),
    seq(max(tf_activities_vector)/paletteLength, 
    max(tf_activities_vector), 
    length.out=floor(paletteLength/2)))

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
    fontsize=14, fontsize_row = 8, fontsize_col = 8, 
    color=myColor, breaks = dorotheaBreaks,
    main = "Dorothea ABC: AT 4h vs basal", angle_col = 45,
    border_color = NA)


tf_activities_counts_filter <- tf_activities_counts %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::filter(GeneID %in% tf_activities_stat_top25_B_vs_RT.4h$GeneID) %>%
    column_to_rownames(var = "GeneID") %>%
    as.matrix()

tf_activities_vector <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
    colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
    length.out=ceiling(paletteLength/2) + 1),
    seq(max(tf_activities_vector)/paletteLength, 
    max(tf_activities_vector), 
    length.out=floor(paletteLength/2)))
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
    fontsize=14, fontsize_row = 8, fontsize_col = 8, 
    color=myColor, breaks = dorotheaBreaks,
    main = "Dorothea ABC: RT 4h vs basal", angle_col = 45,
    treeheight_col = 0,  border_color = NA)

### 5: Network reconstrcution with CARNIVAL
# devtools::install_github("saezlab/CARNIVAL@v1.3")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("OmnipathR", force=TRUE)
# BiocManager::install("CARNIVAL", force=TRUE)
library(OmnipathR)
library(CARNIVAL)

## We also load the support functions
source("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/assignPROGENyScores.r")
source("/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/generateTFList.r")

## We read the normalised counts and the experimental design 
# tf_activities <- tf_activities_CARNIVALinput_RT.4h_vs_b
# PathwayActivity <- PathwayActivity_CARNIVALinput_B_vs_RT.4h
tf_activities <- tf_activities_CARNIVALinput_RT.4h_vs_b
PathwayActivity <- PathwayActivity_CARNIVALinput_B_vs_RT.4h


print("## Getting the scaffold network from Omnipath")
omniR <- import_omnipath_interactions()
  
# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 & #directed with either inhibition or stimulation
                                (consensus_stimulation == 1 | 
                                 consensus_inhibition == 1
                                 ))
  
# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1 
# (converting *both* columns to 1 or -1 - either column indicated stimulation or inhibition) now all inhibition denoted by -1
# stimulation denoted by 1 regardless of column types
omnipath_sd$consensus_stimulation[which(omnipath_sd$consensus_stimulation == 0)] = -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 1)] = -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 0)] = 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
      dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
      unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c('source', 'interaction', 'target')

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)

#save SIF
write_tsv(sif, "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/omnipath_carnival.tsv")



#==================================================================================================================================================================
#==================================================================================================================================================================
#==================================================================================================================================================================
print("### Part 6: Transcription Factor and pathway activities for CARNIVAL")
# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$TF
tf_activities_carnival$TF <- NULL

# generate a list with each TF and its corresponding t value
tfList = generateTFList(tf_activities_carnival, top=100, access_idx = 1)

# progeny for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
PathwayActivity_carnival$Pathway <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                            progenyMembers = progenyMembers, 
                                            id = "gene", 
                                            access_idx = 1)
saveRDS(progenylist, "/home/rza104/scratch/usra/BPK_Work_Study/CARNIVAL/Results/myProgenyList.rds")

current_path <- Sys.getenv("PATH")
new_paths <- "/home/rza104/scratch/cplex/bin/x86-64_linux"
updated_path <- paste(current_path, new_paths, sep=":")
Sys.setenv(PATH = updated_path)
print(Sys.getenv("PATH"))

carnival_result = runCARNIVAL( measObj = data.frame(tfList$t),
                               netObj =  sif, 
                               weightObj = progenylist$score, 
                               solverPath = "/home/rza104/scratch/cplex/bin/x86-64_linux/cplex",
                               solver = "cplex",
                               timelimit = 30000,
                               betaWeight = 0.2)

getwd()
timeStamp = Sys.time()
filename = paste0(format(timeStamp, "%Y%m%d_%H%M%S_"), "carnival.rds")
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)
carnival_result$weightedSIF = carnival_result$weightedSIF[carnival_result$weightedSIF$Weight !=0,]

carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)

saveRDS(carnival_result, file = filename)

filename_weightedSIF = paste0(format(timeStamp, "%Y%m%d_%H%M%S_"), "carnival_weightedSIF.csv")
filename_nodesAttributes = paste0(format(timeStamp, "%Y%m%d_%H%M%S_"), "carnival_nodeAttributes.csv")

write.csv(carnival_result$weightedSIF, filename_weightedSIF)
write.csv(carnival_result$nodesAttributes, filename_nodesAttributes)


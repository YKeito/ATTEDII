####all####
object <- list(ATTEDII_CY20_1hTFDEGs_sort, ATTEDII_CY20_3hTFDEGs_sort, ATTEDII_CY20_12hTFDEGs_sort, ATTEDII_CY20_24hTFDEGs_sort, ATTEDII_CY20_TFDEGsintersection)
Figtitle <- c("ATTEDII_CY20_Table_1h_TFDEGs", "ATTEDII_CY20_Table_3h_TFDEGs", "ATTEDII_CY20_Table_12h_TFDEGs", "ATTEDII_CY20_Table_24h_TFDEGs", "ATTEDII_CY20_Table_TFDEGs_intersection")
g <-c()
library(ggplot2)
i <- 1
for(i in i:length(object)){
  AGI <- rep(object[[i]][, "AGI"], times = 4)
  timecourse <- rep(c("01h", "03h", "12h", "24h"), each = nrow(object[[i]]))
  normalize_degree <- c(object[[i]][, "CY20_1h_normalize_degree"],
                        object[[i]][, "CY20_3h_normalize_degree"],
                        object[[i]][, "CY20_12h_normalize_degree"],
                        object[[i]][, "CY20_24h_normalize_degree"]
  )
  recode <- rep(1:length(object[[i]][, "AGI"]), times = 4)
  data <- data.frame(AGI = AGI,
                     timecourse = timecourse,
                     normalize_degree = normalize_degree,
                     recode = recode
  )
  levels(data$timecourse) <- c("1h", "3h", "12h", "24h")
  levels(data$AGI) <- AGI
  g <- ggplot(data, aes(x = timecourse, y = reorder(AGI, recode), fill = normalize_degree))
  g <- g + geom_tile()
  g <- g + theme_bw()
  g <- g + scale_fill_gradient2(low = "blue", high = "red")
  Figtitle1 <- paste0(Figtitle[i], "NodeNum:", nrow(object[[i]]))
  g <- g + ggtitle(Figtitle1)
  plot(g)
  output <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/network_degree/ATTEDII/heatmap/CY20/all_TFDEGs/", Figtitle1, "_heatmap", ".png")
  ggsave(file = output, plot = g)
  i <- i+1
}
####top5####
object <- list(ATTEDII_CY20_1hTFDEGs_sort[1:5, ], ATTEDII_CY20_3hTFDEGs_sort[1:5, ], ATTEDII_CY20_12hTFDEGs_sort[1:5, ], ATTEDII_CY20_24hTFDEGs_sort[1:5, ], ATTEDII_CY20_TFDEGsintersection[1:5, ])
Figtitle <- c("ATTEDII_CY20_Table_1h_TFDEGs", "ATTEDII_CY20_Table_3h_TFDEGs", "ATTEDII_CY20_Table_12h_TFDEGs", "ATTEDII_CY20_Table_24h_TFDEGs", "ATTEDII_CY20_Table_TFDEGs_intersection")
g <-c()
library(ggplot2)
i <- 1
for(i in i:length(object)){
  AGI <- rep(object[[i]][, "AGI"], times = 4)
  timecourse <- rep(c("01h", "03h", "12h", "24h"), each = nrow(object[[i]]))
  normalize_degree <- c(object[[i]][, "CY20_1h_normalize_degree"],
                        object[[i]][, "CY20_3h_normalize_degree"],
                        object[[i]][, "CY20_12h_normalize_degree"],
                        object[[i]][, "CY20_24h_normalize_degree"]
  )
  recode <- rep(1:length(object[[i]][, "AGI"]), times = 4)
  data <- data.frame(AGI = AGI,
                     timecourse = timecourse,
                     normalize_degree = normalize_degree,
                     recode = recode
  )
  levels(data$timecourse) <- c("1h", "3h", "12h", "24h")
  levels(data$AGI) <- AGI
  g <- ggplot(data, aes(x = timecourse, y = reorder(AGI, recode), fill = normalize_degree))
  g <- g + geom_tile()
  g <- g + theme_bw()
  g <- g + scale_fill_gradient2(low = "blue", high = "red")
  Figtitle1 <- paste0(Figtitle[i], "NodeNum:", nrow(object[[i]]))
  g <- g + ggtitle(Figtitle1)
  plot(g)
  output <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/network_degree/ATTEDII/heatmap/CY20/TFDEGs_top5/", Figtitle1, "_heatmap", ".png")
  ggsave(file = output, plot = g)
  i <- i+1
}

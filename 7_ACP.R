library(dplyr)
library(ggplot2)
library(mclust)
library(maps)
library(FactoMineR)
library(dplyr)
library(tidyr)
library(pracma)
library(here)

# Chemin relatif à la racine du projet
all_seals_data_path <- file.path(here(), "oceanographic_data_corrected", "all_deployment_filtered.rds")

# Charger les données
all_seals_data <- readRDS(all_seals_data_path)

# Maximum Chl-a per profile
df_max_chla <- all_seals_data %>%
  filter(!is.na(chla_corrected)) %>%
  group_by(seal_id, profile_id) %>%
  slice_max(order_by = chla_corrected, n = 1, with_ties = FALSE) %>%  
  ungroup() %>%
  dplyr::select(chla_corrected, light_corrected, temperature, salinity, 
                MLD, SIC_mean, SIC_sd, seal_id, profile_id, date_simple, depth, lat, lon) %>% 
  mutate(SIC_mean = replace_na(SIC_mean, 0), SIC_sd = replace_na(SIC_sd, 0)) %>% 
  na.omit()

head(df_max_chla)

# Computation of integrated mean per profile
df_avg <- all_seals_data %>%
  filter(!is.na(chla_corrected)) %>%
  group_by(seal_id, profile_id, date_simple) %>%
  summarise(
    mean_chla = if(n() > 1) {
      trapz(depth, chla_corrected) / max(depth)
    } else {
      chla_corrected
    },
    .groups = "drop"
  )

head(df_avg)

# Merge max Chl-a and integrated mean Chl-a
df <- df_max_chla %>%
  left_join(df_avg, by = c("seal_id", "profile_id", "date_simple"))

head(df)

# PCA
res_pca=PCA(df[, -c(8:13)],scale=TRUE,ncp=4)
res_pca$var$cos2
res_pca$var$contrib

# > res_pca$var$cos2
#                 Dim.1        Dim.2        Dim.3       Dim.4
# chla_corrected  0.73689990 7.030468e-05 1.489218e-02 0.054463120
# light_corrected 0.08167152 2.175912e-02 5.337051e-01 0.003815076
# temperature     0.33418147 1.435627e-02 4.060742e-05 0.465959517
# salinity        0.48041969 1.716618e-02 1.798463e-02 0.005811424
# MLD             0.10131755 2.024588e-01 4.137989e-01 0.001380953
# SIC_mean        0.01727017 6.358089e-01 3.563688e-03 0.130991677
# SIC_sd          0.07742498 3.659444e-01 1.806421e-01 0.077628543
# mean_chla       0.74320481 8.512698e-03 1.299504e-02 0.060088735

# > res_pca$var$contrib
#                 Dim.1        Dim.2        Dim.3      Dim.4
# chla_corrected  28.6465067  0.005552955  1.264597513  6.8067070
# light_corrected  3.1749275  1.718625623 45.320565822  0.4768016
# temperature     12.9910885  1.133917712  0.003448255 58.2348181
# salinity        18.6760046  1.355856500  1.527198809  0.7263017
# MLD              3.9386541 15.991035032 35.138508750  0.1725891
# SIC_mean         0.6713668 50.218831328  0.302617209 16.3711143
# SIC_sd           3.0098460 28.903812597 15.339565109  9.7018816
# mean_chla       28.8916058  0.672368254  1.103498533  7.5097867

# Visualization
library(factoextra)

fviz_pca_biplot(res_pca, label="var")

fviz_pca_var(res_pca,axes=c(1,2),labelsize = 2,
             col.circle = "grey70",col.var="contrib",repel="TRUE")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=5, space ="Lab") +
  theme_minimal()
#ggsave("plot_2a.png",height=15,width=16,units=c("cm"),dpi=300)

fviz_pca_var(res_pca,axes=c(3,4),labelsize = 2,
             col.circle = "grey70",col.var="contrib",repel="TRUE")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=5, space ="Lab") +
  theme_minimal()

# Clustering
pc1 = res_pca$ind$coord
pc1.mbc = Mclust(pc1, G = 5, modelNames = "EEE")
summary(pc1.mbc)
# Clustering table:
# 1    2   3    4   5 
# 115  55  79  30  252 

# Create pc_scores with identifiers and clusters
pc_scores <- as.data.frame(res_pca$ind$coord) %>%
  mutate(
    seal_id = df$seal_id,
    profile_id = df$profile_id,
    cluster = as.factor(pc1.mbc$classification),
    date_simple = as.Date(df$date_simple),
    month_name = format(date_simple, "%b"),
    month_num  = as.numeric(format(date_simple, "%m")),
    year = format(date_simple, "%Y") 
  )

pc_scores$month_name <- factor(pc_scores$month_name, levels = month.abb, ordered = TRUE)

colnames(pc_scores)[1:2] <- c("PC1", "PC2")

# Create loadings (variable vectors)
loadings <- as.data.frame(res_pca$var$coord[, 1:2])
loadings$varnames <- rownames(loadings)
colnames(loadings)[1:2] <- c("PC1", "PC2")

labels_remanies <- c(
  "chla_corrected" = "FC",
  "light_corrected" = "Light",
  "temperature" = "T°",
  "salinity" = "Sal",
  "MLD" = "MLD",
  "SIC_mean" = "SIC m",
  "SIC_sd" = "SIC sd",
  "mean_chla" = "FC m"
)

# Add a column with modified labels
loadings$label_affiche <- labels_remanies[loadings$varnames]

# Factor to scale up arrows
mult <- 2

# Biplot with clusters and ellipses
ggplot(pc_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +  # transparent points
  # Variable arrows
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  # Variable names
  geom_text(data = loadings, aes(x = PC1 * mult * 1.2, y = PC2 * mult * 1.2, label = label_affiche), color = "black", size = 2) +
  # 95% confidence ellipses
  stat_ellipse(level = 0.95, size = 1) +
  labs(
    title = "Mclust clustering on PCA",
    x = paste0("PC1 (", round(summary(res_pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(res_pca)$importance[2,2]*100, 1), "%)")
  ) +
  theme_minimal()

# PCA biplot with monthly ellipses
ggplot(pc_scores, aes(x = PC1, y = PC2, color = month_name)) +
  geom_point(size = 0.7, alpha = 0.5) +
  stat_ellipse(aes(group = month_name, color = month_name), level = 0.95, size = 1) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text(data = loadings, 
            aes(x = PC1 * mult * 1.2, y = PC2 * mult * 1.2, label = label_affiche), 
            color = "black", size = 2, inherit.aes = FALSE) +
  labs(
    title = "PCA with ellipses by month (seasonal cycle)",
    x = paste0("PC1 (", round(summary(res_pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(res_pca)$importance[2,2]*100, 1), "%)")
  ) +
  theme_minimal()

# PCA biplot with yearly ellipses
ggplot(pc_scores, aes(x = PC1, y = PC2, color = year)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text(data = loadings, aes(x = PC1 * mult * 1.2, y = PC2 * mult * 1.2, label = label_affiche), 
            color = "black", size = 2, inherit.aes = FALSE) +
  stat_ellipse(aes(group = year, color = year), level = 0.95, size = 1) +
  labs(
    title = "PCA with ellipses by year",
    x = paste0("PC1 (", round(summary(res_pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(res_pca)$importance[2,2]*100, 1), "%)")
  ) +
  theme_minimal()

# Geographical projection of clusters
df_map <- df %>%
  select(seal_id, profile_id, lat, lon, depth) %>% distinct() %>%
  left_join(pc_scores %>% select(seal_id, profile_id, cluster, month_name), by = c("seal_id","profile_id"))

min_lon <- 136
max_lon <- 146
min_lat <- -67.3
max_lat <- -65.4

world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray90", color = "gray70") +
  geom_point(data = df_map, aes(x = lon, y = lat, color = cluster), size = 0.5, alpha = 0.8) +
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) +
  labs(title = "Spatial projection of clusters", x = "Longitude", y = "Latitude") +
  theme_minimal()

### ---------------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(maps)

df_map <- df_map %>%
  mutate(cluster = factor(cluster))

world <- map_data("world")

# Plot with facets by cluster
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray90", color = "gray70") +
  geom_point(data = df_map, aes(x = lon, y = lat, color = cluster), size = 0.5, alpha = 0.8) +
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) +
  facet_wrap(~cluster, nrow = 2, ncol = 3) +
  labs(title = "Spatial projection of clusters", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "bottom")

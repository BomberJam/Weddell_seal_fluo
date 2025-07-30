library(dplyr)
library(ggplot2)
library(viridis)


# Charger les données
df <- readRDS("/home/cactus/Documents/Oceano/M2/LOCEAN/Weddell_seal_fluo/oceanographic_data_corrected/all_deployment_filtered.rds")

df <- df %>%
  mutate(seal_id = ifelse(is.na(seal_id) | seal_id == "", SEAL_NAME, seal_id))
df <- df %>% dplyr::select(-SEAL_NAME)
df <- df[!apply(is.na(df), 1, all), ]

# Sélectionner les variables utiles + seal_id (SIC et Fast-ice)
df_acp <- df %>% dplyr::select(seal_id, temperature, salinity, MLD, chla_corrected, depth, light_corrected) %>%
  na.omit()

# Garder seal_id à part et centrer les données
seal_names <- df_acp$seal_id
df_numeric <- df_acp %>% dplyr::select(-seal_id)

# Appliquer l'ACP
res_pca <- prcomp(df_numeric, scale. = TRUE)

# Créer un data frame avec les coordonnées des individus
pca_df <- as.data.frame(res_pca$x)  # contient PC1, PC2, etc.
pca_df$seal_id <- seal_names      # ajouter seal_id

# Données pour les vecteurs des variables
loadings <- as.data.frame(res_pca$rotation[, 1:2])  # PC1 et PC2
loadings$varnames <- rownames(loadings)

# Facteur pour bien voir les flèches
mult <- 3  

centroids <- pca_df %>%
  group_by(seal_id) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# Visualiser avec ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = seal_id)) +
  
  stat_ellipse(aes(linetype = seal_id, color = seal_id), type = "norm", size = 1) +
  # stat_ellipse(type = "norm", linetype = 2) +
  
  geom_point(data = centroids, aes(x = PC1, y = PC2, fill = seal_id), size = 4, shape = 21, color = "black")  +
  # geom_point(alpha = 0.7, size = 1) +

  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  # Noms des variables au bout des flèches
  geom_text(data = loadings, aes(x = PC1 * mult * 1.1, y = PC2 * mult * 1.1, label = varnames),
            color = "black", size = 4) +
  labs(title = "ACP wd",
       x = paste0("PC1 (", round(summary(res_pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(res_pca)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# -------------------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(lubridate)  # pour month()

# Charger les données
df <- readRDS("/home/cactus/Documents/Oceano/M2/LOCEAN/Weddell_seal_fluo/oceanographic_data_corrected/all_deployment_filtered.rds")

df <- df %>%
  mutate(seal_id = ifelse(is.na(seal_id) | seal_id == "", SEAL_NAME, seal_id))
df <- df %>% dplyr::select(-SEAL_NAME)
df <- df[!apply(is.na(df), 1, all), ]

# Vérifier que la colonne date existe et est bien au format Date
df$date <- as.Date(df$date)  # ou POSIXct si nécessaire

# Sélectionner les colonnes utiles pour l’ACP
df_acp <- df %>%
  dplyr::select(seal_id, temperature, salinity, MLD, chla_corrected, depth, light_corrected, date) %>%
  na.omit()

# Extraire seal_id et mois
seal_names <- df_acp$seal_id
months <- month(df_acp$date, label = TRUE, abbr = TRUE)  # factor avec labels "Jan", "Feb", etc.
df_numeric <- df_acp %>% dplyr::select(-seal_id, -date)

# ACP
res_pca <- prcomp(df_numeric, scale. = TRUE)

# Résultats PCA + mois
pca_df <- as.data.frame(res_pca$x)
pca_df$seal_id <- seal_names
pca_df$month <- months  # facteur des mois

# vecteurs
loadings <- as.data.frame(res_pca$rotation[, 1:2])
loadings$varnames <- rownames(loadings)
mult <- 3

centroids <- pca_df %>%
  group_by(month) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = month)) +
  
  stat_ellipse(aes(linetype = month, color = month), type = "norm", size = 1) +

  geom_point(data = centroids, aes(x = PC1, y = PC2, fill = month), size = 4, shape = 21, color = "black")  +
  # geom_point(alpha = 0.6, size = 0.5) +
  # scale_color_viridis_d(name = "Mois", option = "A") +
  # stat_ellipse(type = "norm", linetype = 2, color = "grey40") +
  geom_segment(data = loadings, aes(x = 0, y = 0,
                                    xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text(data = loadings, aes(x = PC1 * mult * 1.1,
                                 y = PC2 * mult * 1.1,
                                 label = varnames), color = "black", size = 4) +
  labs(title = "ACP (monthly)",
       x = paste0("PC1 (", round(summary(res_pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(res_pca)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9))


# -------------------------------------------------------------------------------------------

df <- readRDS("/home/cactus/Documents/Oceano/M2/LOCEAN/Weddell_seal_fluo/oceanographic_data_corrected/all_deployment_filtered.rds")

df <- df %>%
  mutate(seal_id = ifelse(is.na(seal_id) | seal_id == "", SEAL_NAME, seal_id))
df <- df %>% dplyr::select(-SEAL_NAME)
df <- df[!apply(is.na(df), 1, all), ]

# Extraire la valeur max(Chla), extraire la profondeur, Temp, sal, SIC (moyenne)
df_max_chla <- df %>%
  filter(!is.na(chla_corrected)) %>%  # on évite les NA
  group_by(seal_id, profile_id) %>%
  slice_max(chla_corrected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(seal_id, profile_id, chla_corrected, light_corrected, depth, temperature, salinity, MLD)

# -------------------------------------------------------------------------------------------
#Table de résultats
library(knitr)

kable(df_max_chla, digits = 3)

df_max_chla %>%
  arrange(desc(chla_corrected)) %>%
  head(10) %>%
  kable(digits = 3)

# -------------------------------------------------------------------------------------------
# ACP sur max(chla)
library(dplyr)
library(ggplot2)

# Sélection des colonnes numériques utiles
df_acp <- df_max_chla %>%
  dplyr::select(seal_id, chla_corrected, light_corrected, depth, temperature, salinity, MLD) %>%
  na.omit()  # ACP ne supporte pas les NA

# Vérifier et retirer colonnes constantes
seal_names <- df_acp$seal_id
df_numeric <- df_acp %>% dplyr::select(-seal_id)

# Lancer l’ACP
res_pca <- prcomp(df_numeric, scale. = TRUE)
pca_df <- as.data.frame(res_pca$x)  # contient PC1, PC2, etc.
pca_df$seal_id <- seal_names

# Résumé de l’ACP
summary(res_pca)

# Données pour les vecteurs des variables
loadings <- as.data.frame(res_pca$rotation[, 1:2])  # PC1 et PC2
loadings$varnames <- rownames(loadings)

# Facteur pour bien voir les flèches
mult <- 3  

# Graphe PCA avec points et vecteurs
ggplot(pca_df, aes(x = PC1, y = PC2, color = seal_id)) +
  geom_point(size = 0.5, alpha = 0.7) +
  # Flèches des variables
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  # Noms des variables au bout des flèches
  geom_text(data = loadings, aes(x = PC1 * mult * 1.1, y = PC2 * mult * 1.1, label = varnames),
            color = "black", size = 4) +
  labs(title = "ACP max(Chla)",
       x = paste0("PC1 (", round(summary(res_pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(res_pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()

# -------------------------------------------------------------------------------------------





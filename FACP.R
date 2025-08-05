library(dplyr)
library(tidyr)
library(fdapace)
library(plot3D)

df <- readRDS("/home/cactus/Documents/Oceano/M2/LOCEAN/Weddell_seal_fluo/oceanographic_data_corrected/all_deployment_filtered.rds")
# Supprimer lignes entièrement vides
df <- df[!apply(is.na(df), 1, all), ]
# Variables qu'on garde pour la FACP
varname <- "temperature"  # ou "salinity", etc.

# On garde juste les colonnes utiles
df <- df %>%
  select(profile_id, depth, temperature) %>%
  filter(!is.na(depth), !is.na(temperature)) %>%
  arrange(depth, profile_id)

# Identifier les profondeurs et profils
depths <- sort(unique(df$depth))
ndepth <- length(depths)
profiles <- unique(df$profile_id)
nprofiles <- length(profiles)

# Construire la matrice : lignes = profondeurs, colonnes = profils
mat <- matrix(NA, nrow = ndepth, ncol = nprofiles)
for (i in seq_along(profiles)) {
  sel <- df %>% filter(profile_id == profiles[i]) %>% arrange(depth)
  mat[match(sel$depth, depths), i] <- sel[[varname]]
}
#mat_debug <- cbind(depth = depths, mat) # pour le debug
matplot(depths,mat,type="l",xlab="Depth",ylab="T°")

# ---- Normalisation par profondeur ----
nb.obs.na <- apply(is.na(mat), 1, sum)
mat.m <- apply(mat, 1, mean, na.rm = TRUE)
mat.m <- sweep(mat, 1, mat.m, "-")
mat.sd <- apply(mat, 1, sd, na.rm = TRUE)
anom <- sweep(mat.m, 1, mat.sd, "/")
anom.m <- apply(anom, 2, mean, na.rm = TRUE)

# ---- Visualisation anomalies ----
matplot(depths, anom, type = "l", lty = 1, las = 1, xlab = "Depth (m)", ylab = "Normalized anomaly", col = "darkgrey")
lines(depths, anom.m, lwd = 2)
abline(h = 0, lty = 2)

# ---- FPCA ----
mat <- t(anom)  # colonnes = profondeurs, lignes = profils
depthmat <- matrix(depths, ndepth, nsites)

neval <- 100 # output size
# Construire les listes
Ly <- lapply(1:ncol(mat), function(i) mat[, i])          # valeurs anomalies par profil
Lt <- lapply(1:ncol(mat), function(i) depths)            # profondeurs pour chaque profil

# Lancer FPCA
res <- FPCA(Ly = Ly, Lt = Lt,
            optns = list(
              useBinnedData = 'OFF',
              methodMuCovEst = 'smooth', #The method to estimate the mean and covariance in the case of dense functional data
              userBwCov = 5, #The bandwidth value for the smoothed covariance function
              userBwMu = 7, #The bandwidth value for the smoothed mean function 
              kernel = "epan", #Smoothing kernel choice
              nRegGrid = neval, #The number of support points in each direction of covariance surface
              dataType = 'Sparse', #The type of design we have
              methodXi = 'CE' #The method to estimate the PC scores
            ))
plot(res)

library(fields)

nHalf = sum(res$fittedCorr < 0)  
npix = neval*neval
Min = -1
Max = 1
Thresh = 0

## Couleurs
rc1 = colorRampPalette(colors = c("red", "white"), space = "Lab")(nHalf)    
nc2 = npix - nHalf
rc2 = colorRampPalette(colors = c("white", "blue"), space = "Lab")(nc2)
rampcols = c(rc1, rc2)

rb1 = seq(Min, Thresh, length.out = nHalf + 1)
rb2 = seq(Thresh, Max, length.out = npix - nHalf + 1)[-1]
rampbreaks = c(rb1, rb2)

# Affichage
fields::image.plot(
  res$workGrid, res$workGrid, res$fittedCorr[,neval:1],
  col = rampcols, breaks = rampbreaks, las=1,
  xlab = "Depth (m)", ylab = "Depth (m)",
  main = "Correlation matrix (Pearson)",
  cex.main=1,axes=FALSE,frame.plot=TRUE )

contour(res$workGrid,res$workGrid,res$fittedCorr[,neval:1],add=TRUE,labcex=1)
lab <- depths
pos <- seq(min(lab), max(lab), length.out = 40)
pos <- floor(pos)
axis(1,at=pos,labels=pos,las=2)
axis(2,at=rev(pos),labels=pos,las=2)
box()
grid()


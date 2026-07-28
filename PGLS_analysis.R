require(ape)
require(MCMCglmm)

#### Load FMR database ####
dir <- "C:/..."
data <- read.csv(paste0(dir,"/FMR_database.csv"), sep=",")

#### Load Pantheria ####
pantheria <- read.table(paste0(dir,"/PanTHERIA_1-0_WR05_Aug2008.txt"), header=T, row.names=NULL)
pantheria[which(pantheria==-999, arr.ind = T)]<-NA
plot(log(X18.1_BasalMetRate_mLO2hr) ~ log(X5.2_BasalMetRateMass_g), pantheria)
mod_BMR <- lm(log(X18.1_BasalMetRate_mLO2hr) ~ log(X5.2_BasalMetRateMass_g), pantheria, na.action = "na.exclude")
pantheria$BMR_M <- residuals(mod_BMR)
pant_traits <- data.frame(
  Species = paste0(pantheria$MSW05_Species,"_", pantheria$MSW05_Genus),
  ActivityCycle = pantheria$X1.1_ActivityCycle, #  (1) nocturnal only, (2) nocturnal/crepuscular, cathemeral, crepuscular or diurnal/crepuscular and (3) diurnal only. 
  BMR_M = pantheria$BMR_M,
  HomeRange = pantheria$X22.2_HomeRange_Indiv_km2,
  TrohpicLev = pantheria$X6.2_TrophicLevel #  (1) herbivore, (2) omnivore and (3) carnivore
)
data$ActivityCycle <- as.factor(pant_traits$ActivityCycle[match(data$Species, pant_traits$Species)])
data$BMR_M <- pant_traits$BMR_M[match(data$Species, pant_traits$Species)]
data$HomeRange <- pant_traits$HomeRange[match(data$Species, pant_traits$Species)]
data$TrohpicLev <- as.factor(pant_traits$TrohpicLev[match(data$Species, pant_traits$Species)])

plot(FMR_M ~ ActivityCycle, data) # (1) nocturnal only, (2) nocturnal/crepuscular, cathemeral, crepuscular or diurnal/crepuscular and (3) diurnal only. 
plot(FMR_M ~ BMR_M, data); abline(0,1)
plot(FMR_M ~ log(HomeRange), data)
plot(FMR_M ~ TrohpicLev, data) #  (1) herbivore, (2) omnivore and (3) carnivore

#### Non-phylo model selection ####
require(lme4)
require(lmerTest)

data_mammals_pgls <- data.frame(Species=data$Species, 
                                FMR_Watt=data$FMR_Watt, 
                                FMR_M=data$FMR_M, 
                                Mass=data$Mass,
                                TALOCmean=data$TALOCmean,
                                TMEAN=data$TMEAN, 
                                Tanomalies=data$Tanomalies 
                                #HomeRange=data$HomeRange, # un-comment these variables (from PanTHERIA) for the analysis with reduced database (n=211)
                                #ActivityCycle=data$ActivityCycle,
                                #BMR_M=data$BMR_M,
                                #TrophLev=data$TrohpicLev
)
data_mammals_pgls <- data_mammals_pgls[complete.cases(data_mammals_pgls),] # Remove NAs

# Including variable "Torpor" to identify individuals that displayed torpor during EWL measurements (total of 19 individuals of Microcebus murinus) 
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 61.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 59.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 43.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 55.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 38.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 54.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 47.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 86.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 97.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 84.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 79.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2000). Daily energy expenditure of the grey mouse lemur (Microcebus murinus): A small primate that uses torpor. Journal of Comparative Physiology B, 170(8), 633-641. https://doi.org/10.1007/s003600000146" & data_mammals_pgls$Mass == 75.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 63.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 51.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 52.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 60.5 & round(data_mammals_pgls$FMR_Watt) == 1] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 65.5] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 64.0] <- "YES"
data_mammals_pgls$Torpor[data_mammals_pgls$Article == "Schmid, J., & Speakman, J. R. (2009). Torpor and energetic consequences in free-ranging grey mouse lemurs (Microcebus murinus): A comparison of dry and wet forests. Naturwissenschaften, 96(5), 609-620. https://doi.org/10.1007/s00114-009-0515-z" & data_mammals_pgls$Mass == 46.8] <- "YES"

# remove individuals with torpor=yes
data_mammals_pgls <- data_mammals_pgls[-which(data_mammals_pgls$Torpor == "YES"),]

length(data_mammals_pgls$Species) # Obs
length(unique(data_mammals_pgls$Species)) # Spp

phylo <- read.tree(paste0(dir,"/FritzTree_mammals_consensus.tre")) # Consensus tree for mammals (Fritz et al. 2009)
data_mammals_pgls$Species[which(data_mammals_pgls$Species=="Eremitalpa_granti")] <- "Chrysochloris_stuhlmanni"

data_mammals_pgls$animal <- "not in tree"
for(i in 1:nrow(data_mammals_pgls)){
  species <- data_mammals_pgls$Species[i]
  if(rlang::is_empty(phylo$tip.label[grepl(species,phylo$tip.label)])){ # if sp is not in tree leave "not in tree"
  } else {
    data_mammals_pgls$animal[i]<-phylo$tip.label[grepl(species,phylo$tip.label)]} # else put the sp name from the tree
}

length(unique(data_mammals_pgls$animal[which(!data_mammals_pgls$animal=="not in tree")])) # species in tree
length(unique(data_mammals_pgls$Species[which(data_mammals_pgls$animal=="not in tree")])) # not in tree

nameslist <- phylo$tip.label
treenameslist <- as.data.frame(table(data_mammals_pgls$animal))
Speciestoretain <- intersect(treenameslist$Var1, nameslist)
pruned.tree <- drop.tip(phylo,phylo$tip.label[-match(Speciestoretain,phylo$tip.label)])
plot(pruned.tree, cex=0.5)
tree <- pruned.tree
tree <- compute.brlen(tree)
tree$node.label<-NULL
is.ultrametric(tree) # checks
is.rooted(tree)
any(duplicated(tree$node.label))

# Re-order dataframe to match species in data and tree labels
data_mammals_pgls_sub <- data_mammals_pgls[(data_mammals_pgls$animal %in% tree$tip.label),]
data_mammals_pgls_sub <- data_mammals_pgls[which(!data_mammals_pgls$animal=="not in tree"),]
length(unique(data_mammals_pgls_sub$animal))
nrow(data_mammals_pgls_sub)

###### Temperature and anomaly ####
require(lme4)
require(lmerTest)
require(MuMIn)
require(performance)

# models
mod_null <- MCMCglmm(FMR_M ~ 1, 
                     random=~Species, data=data_mammals_pgls_sub, 
                     pedigree=tree, nitt=10000)
mod_Tmean <- MCMCglmm(FMR_M ~ TMEAN, 
                      random=~Species, data=data_mammals_pgls_sub, 
                      pedigree=tree, nitt=10000)
mod_Tmean2 <- MCMCglmm(FMR_M ~ TMEAN + I(TMEAN^2), 
                       random=~Species, data=data_mammals_pgls_sub, 
                       pedigree=tree, nitt=10000)
mod_Tanom <- MCMCglmm(FMR_M ~ Tanomalies, 
                      random=~Species, data=data_mammals_pgls_sub, 
                      pedigree=tree, nitt=10000)
mod_Tanom2 <- MCMCglmm(FMR_M ~ Tanomalies + I(Tanomalies^2), 
                       random=~Species, data=data_mammals_pgls_sub, 
                       pedigree=tree, nitt=10000)
mod_Tmean_Tanom <- MCMCglmm(FMR_M ~ TMEAN + Tanomalies, 
                            random=~Species, data=data_mammals_pgls_sub, 
                            pedigree=tree, nitt=10000)
mod_Tmean2_Tanom <- MCMCglmm(FMR_M ~ TMEAN + Tanomalies + I(TMEAN^2), 
                             random=~Species, data=data_mammals_pgls_sub, 
                             pedigree=tree, nitt=10000)
mod_Tmean_Tanom2 <- MCMCglmm(FMR_M ~ TMEAN + Tanomalies + I(Tanomalies^2), 
                             random=~Species, data=data_mammals_pgls_sub, 
                             pedigree=tree, nitt=10000)
mod_Tmean2_Tanom2 <- MCMCglmm(FMR_M ~ TMEAN + Tanomalies + I(TMEAN^2) + I(Tanomalies^2), 
                              random=~Species, data=data_mammals_pgls_sub, 
                              pedigree=tree, nitt=10000)
mod_TmeanXTanom <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies, 
                            random=~Species, data=data_mammals_pgls_sub, 
                            pedigree=tree, nitt=10000)
mod_TmeanXTanom2 <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2), 
                             random=~Species, data=data_mammals_pgls_sub, 
                             pedigree=tree, nitt=10000)
mod_Tmean2XTanom <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(TMEAN^2), 
                             random=~Species, data=data_mammals_pgls_sub, 
                             pedigree=tree, nitt=10000)
mod_Tmean2XTanom2 <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2), 
                              random=~Species, data=data_mammals_pgls_sub, 
                              pedigree=tree, nitt=10000)

# DIC comparison
mod_null$DIC
mod_Tmean$DIC
mod_Tmean2$DIC
mod_Tanom$DIC
mod_Tanom2$DIC
mod_Tmean_Tanom$DIC
mod_Tmean2_Tanom$DIC
mod_Tmean_Tanom2$DIC
mod_Tmean2_Tanom2$DIC
mod_TmeanXTanom$DIC 
mod_TmeanXTanom2$DIC 
mod_Tmean2XTanom$DIC 
mod_Tmean2XTanom2$DIC

# Check model MCMC
plot(mod_Tmean2XTanom2)
summary(mod_TmeanXTanom2)

###### Including covariates ####

mod_Tmean2XTanom2 <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2), 
                             random=~Species, data=data_mammals_pgls_sub, 
                             pedigree=tree, nitt=10000)
mod_Tmean2XTanom2_TL <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2) +
                              TrophLev, 
                              random=~Species, data=data_mammals_pgls_sub, 
                              pedigree=tree, nitt=10000)
mod_Tmean2XTanom2_HR <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2) +
                              log(HomeRange)+I(log(HomeRange)^2), 
                              random=~Species, data=data_mammals_pgls_sub, 
                              pedigree=tree, nitt=10000)
mod_Tmean2XTanom2_AC <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2) +
                              ActivityCycle, 
                              random=~Species, data=data_mammals_pgls_sub, 
                              pedigree=tree, nitt=10000)
mod_Tmean2XTanom2_AC_HR <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2) +
                                 ActivityCycle + log(HomeRange)+I(log(HomeRange)^2), 
                                 random=~Species, data=data_mammals_pgls_sub, 
                                 pedigree=tree, nitt=10000)

# DIC comparison
mod_Tmean2XTanom2$DIC
mod_Tmean2XTanom2_TL$DIC
mod_Tmean2XTanom2_HR$DIC 
mod_Tmean2XTanom2_AC$DIC
mod_Tmean2XTanom2_AC_HR$DIC

# Check model MCMC
plot(mod_Tmean2XTanom2_AC)
summary(mod_Tmean2XTanom2_AC)
summary(mod_Tmean2XTanom2_TL)

# Parameters of the best model (without covariates)
best_model <- function(Tmean, Tanom){
  0.6065489-0.0647042*Tmean+0.1043187*Tanom+0.0202398*Tanom^2+0.0009998*Tmean^2-0.0098475*Tanom*Tmean
}

Tcold_anom_0 <- best_model(Tmean = 11.10, Tanom = 0)
Tcold_anom_minus1 <- best_model(Tmean = 11.10, Tanom = -1)
Tcold_anom_plus1 <- best_model(Tmean = 11.10, Tanom = 1)

Tcold_anom_minus1/Tcold_anom_0
Tcold_anom_plus1/Tcold_anom_0

Twarm_anom_0 <- best_model(Tmean = 23.83, Tanom = 0)
Twarm_anom_minus1 <- best_model(Tmean = 23.83, Tanom = -1)
Twarm_anom_plus1 <- best_model(Tmean = 23.83, Tanom = 1)

(Twarm_anom_minus1+1)/(Twarm_anom_0+1)
(Twarm_anom_plus1+1)/(Twarm_anom_0+1)

#### PLOTS BEST MODEL ####

require(ggplot2)
require(RColorBrewer)
summary(data_mammals_pgls_sub$TMEAN) # Observed ranges of mean T

# Separate data below and above the mean T
data_plot_lowT <- data_mammals_pgls_sub[data_mammals_pgls_sub$TMEAN<17.04,]
data_plot_highT <- data_mammals_pgls_sub[data_mammals_pgls_sub$TMEAN>17.04,]

# Relación FMR vs mean T
ggplot() + theme_classic() + #500 x 400
  theme(axis.text = element_text(colour = "black", size=15),
        axis.title = element_text(size=18)) +
  ylab("Mass-specific Field Metabolic Rate") + xlab("Mean Temperature (ºC)") +
  geom_point(data_plot_lowT, mapping=aes(y=FMR_M, x=TMEAN), col="#4393C3", size=2) + 
  geom_point(data_plot_highT, mapping=aes(y=FMR_M, x=TMEAN), col="#B2182B", size=2) + 
  geom_function(fun=function(x) best_model(Tanom=0,Tmean=x), lwd=1, col="black") 

# FMR vs anomaly
ggplot() + theme_classic() +  #500 x 400
  theme(axis.text = element_text(colour = "black", size=15),
        axis.title = element_text(size=18)) +
  ylab("Mass-specific Field Metabolic Rate") + xlab("Temperature anomaly (ºC)") +
  geom_point(data_plot_lowT, mapping=aes(y=FMR_M, x=Tanomalies), col="#4393C3", size=2) + 
  geom_point(data_plot_highT, mapping=aes(y=FMR_M, x=Tanomalies), col="#B2182B", size=2) + 
  geom_vline(xintercept = 0, lty="dashed")+
  geom_function(fun=function(x) best_model(Tanom=x,Tmean=11.10), lwd=1, col="#4393C3") +
  geom_function(fun=function(x) best_model(Tanom=x,Tmean=23.83), lwd=1, col="#B2182B") 

FMR_heatmap <- array(NA, dim=c(100,100))
Tmean <- seq(-10,30,length.out=100)
Tanom <- seq(-6,6, length.out=100)
for(i in 1:100){
  for(j in 1:100){
    FMR_heatmap[i,j] <- best_model(Tmean = Tmean[i], Tanom = Tanom[j])
  }
}
require(reshape2)
require(RColorBrewer)
require(metR)
FMR_plot <- melt(FMR_heatmap)
FMR_plot$Var1 <- rep(Tmean, 100)
FMR_plot$Var2 <- sort(rep(Tanom, 100))
ggplot() + 
  theme_bw() + theme(axis.text = element_text(colour = "black", size=10),
                     axis.title = element_text(size=12),
                     legend.title = element_text(angle=90),
                     legend.title.position = "right",
                     panel.grid = element_blank()) +
  ylab("Temperature anomaly (ºC)") + xlab("Mean Temperature (ºC)") + labs(fill="Mass-specific FMR") +
  metR::geom_contour_fill(FMR_plot, mapping=aes(x=Var1, y=Var2, z=value), bins=10) +
  metR::geom_contour2(FMR_plot, mapping=aes(x=Var1, y=Var2, z=value), skip=1, bins=10) +
  # metR::geom_contour2(aes(z = value, label = stat(level)), skip=1, bins=10) +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_rug(data_mammals_pgls_sub, mapping=aes(y=Tanomalies, x=TMEAN), size=1) +
  # geom_point(datos_mammals_pgls_sub, mapping=aes(y=Tanomalies, x=TMEAN)) +
  scale_fill_gradientn(colors = brewer.pal(8, "Reds"))

# Additional pannel (R1)
summary(data_mammals_pgls_sub$TMEAN) # Observed ranges of mean T

data_plot_lowT <- data_mammals_pgls_sub[data_mammals_pgls_sub$TMEAN<17.04,]
data_plot_highT <- data_mammals_pgls_sub[data_mammals_pgls_sub$TMEAN>17.04,]

ggplot() + theme_classic() +  
  theme(axis.text = element_text(colour = "black", size=15),
        axis.title = element_text(size=18)) +
  ylab("Mass-specific Field Metabolic Rate") + xlab("Temperature anomaly (ºC)") +
  geom_point(data_plot_lowT, mapping=aes(y=FMR_M, x=Tanomalies), col="blue", size=2) + 
  geom_point(data_plot_highT, mapping=aes(y=FMR_M, x=Tanomalies), col="red", size=2) + 
  
  geom_vline(xintercept = 0, lty="dashed")+
  geom_function(fun=function(x) best_model(Tanom=x,Tmean=11.10), lwd=1, col="blue") +
  geom_function(fun=function(x) best_model(Tanom=x,Tmean=23.83), lwd=1, col="red") 

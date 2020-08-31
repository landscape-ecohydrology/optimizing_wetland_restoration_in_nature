# ===============================================================================================
# Decoupled Wetland Removal 
# DATE: 06 June 2019
#
# This function is used to generate: 
# a) the N removal potential of the wetlands in a HUC region
# b) the N mass removed of the wetlands in a HUC region
# within a Monte Carlo framework and outputs the relevant summary statistics of the MCA 

# ===============================================================================================
HUC_ProcessWetlandSurplus.fun <- function(HUC_NUM = '07') {
  DataWD <- "D:/OneDrive - University of Waterloo/Research/Wetland Loss in the US/Data Large/"

# HUC Boundary Spatial -------------------------------------------------------------------------
setwd(paste('D:/GIS DATA FILES/Watershed Boundary Dataset/Reprojected/HUC8_',HUC_NUM, sep = ""))
HUC_Boundary <-  st_read(paste('HUC8_', HUC_NUM, '.shp', sep =""))  # read HUC data
HUC_Boundary <- dplyr::select(HUC_Boundary, c(HUC8))                # keep only two columns
HUC_Boundary$HUC_area_m2 <- as.vector(st_area(HUC_Boundary))        # calculate areas of each HUC in m2

#N_SURPLUS --------------------------------------------------------------------------------------
N_SURPLUS_DATA  <- readr::read_csv("D:/OneDrive - University of Waterloo/Research/Wetland Loss in the US/Data Large/N Surplus/huc_n_surplus.csv")
colnames(N_SURPLUS_DATA) <- c('HUC8', 'NSUR_kgha')

HUC_Boundary <- HUC_Boundary %>%
  left_join(N_SURPLUS_DATA) %>%
  mutate(N_SUR_kg = NSUR_kgha * HUC_area_m2 / 10000)

# Wetland Shapes Spatial -------------------------------------------------------------------------
wetlands <- readRDS(paste(DataWD,'HUC_Wetlands/HUC8_',HUC_NUM, '.rds', sep = "")) %>%
  mutate(BIN = as.factor(floor(log10(wet_area_m2)))) %>%
  filter(WETLAND_TY != 'Lake') %>%
  filter(WETLAND_TY != 'Estuarine and Marine Deepwater') 

# Perform Monte Carlo Simulations for each HUC ----------------------------------------------------
# Generate a set of n values for tau, k, CA, R for each wetland -----------------------------------


HUC8_Index <-as.character(unique(wetlands$HUC8))

RemovalsPerHUC8 <- function(index, wetlands_df, HUC_Boundary_df){
  num_sim = 250
  wetlands_HUC <- filter(wetlands, HUC8 == index) %>% # Test for one HUC8 watershed 
    na.omit()
  
  TAU_MC <- replicate(num_sim, apply(wetlands_HUC, 1, function(x) 10^runif(1,0.09097,0.2909 ) 
                                                         * (as.numeric(x[3]) ^ runif(1,0.2083,0.2455) )))
                  
  K_MC   <- TAU_MC ^ runif(1,-0.9486, -0.8491) * 10 ^ runif(1,-0.535,-0.3831)

  R_MC <- 1-exp(-TAU_MC * K_MC)
  remove(TAU_MC, K_MC)
  
  CA_MC  <- replicate(num_sim, apply(wetlands_HUC, 1, function(x) runif(1, 3.5,20.16) * (as.numeric(x[3]) )))
    
    
  CA_RATIO <- colSums(CA_MC) / HUC_Boundary$HUC_area_m2[HUC_Boundary$HUC8 == index]
  CA_RATIO[CA_RATIO<1] <- 1  
  CA_MC <- CA_MC / CA_RATIO
  
  WeightedRp <- colSums(CA_MC*R_MC) / HUC_Boundary$HUC_area_m2[HUC_Boundary$HUC8 == index]
  
  NSUR_to_wet_kg_MC <-  CA_MC * runif(1, 0.4,0.6) * HUC_Boundary$NSUR_kgha[HUC_Boundary$HUC8 == index]/10000

  Ra_HUC <- colSums(R_MC * NSUR_to_wet_kg_MC)
  
  OUTPUT <- data.frame(HUC8 = index, NumWet = nrow(wetlands_HUC), 
                         sum_wetSA = sum(wetlands_HUC$wet_area_m2), sum_CA = median(colSums(CA_MC)),
              Rp5  = quantile(WeightedRp, 0.05),
              Rp25 = quantile(WeightedRp, 0.25),
              Rp50 = quantile(WeightedRp, 0.50),
              Rp75 = quantile(WeightedRp, 0.75),
              Rp95 = quantile(WeightedRp, 0.95),
              MeanRp = mean(WeightedRp), 
              sdRp = sd(WeightedRp),
              Ra5  = quantile(Ra_HUC, 0.05),
              Ra25 = quantile(Ra_HUC, 0.25),
              Ra50 = quantile(Ra_HUC, 0.50),
              Ra75 = quantile(Ra_HUC, 0.75),
              Ra95 = quantile(Ra_HUC, 0.95),
              MeanRa = mean(Ra_HUC), 
              sdRa = sd(Ra_HUC))

 
}

start_time <- Sys.time()
HUC8_Removals <- do.call('rbind', lapply(HUC8_Index, FUN =  RemovalsPerHUC8, wetlands_df = wetlands, 
                                          HUC_Boundary_df = HUC_Boundary))
end_time <- Sys.time()
end_time - start_time


HUC_OUTPUT <- HUC8_Removals %>%
  left_join(HUC_Boundary) %>%
  mutate(CA_RATIO = sum_CA / HUC_area_m2, HUC2 = HUC_NUM)

saveRDS(HUC_OUTPUT, paste('D:/OneDrive - University of Waterloo/Research/Wetland Loss in the US/WetlandLoss/Outputs/Wetlands in Nature/MonteCarlo/OUTPUT_RpRa_HUC8_', HUC_NUM, '.rds', sep=""))


}



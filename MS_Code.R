## Header and Preamble ------------------------------------------------------------
##
## Script name: MS Analyses
##
## Purpose of script: Clean code to reproduce MS model results and analyses
##
##

## Notes:
##
##
## # FIX BEFORE GITHUB
##

## Load up the packages we will need:

library(gratia) # For looking at model directionality
library(tidyverse)

##
######## End


# Load in and prep data ------------------------------------------------------------

MegaDS <- 
  read.csv("data/exported_data/summary_Chamber.csv")

NEP <- 
  MegaDS %>% 
    filter(Treatment == "Light") %>% 
    filter(!Species == "ASU") %>% 
    filter(!Species == "Sediment") %>% 
    mutate(Site = as.factor(Site),
           Date = as.factor(Date)) 

ER <- 
  MegaDS %>% 
  filter(Treatment == "Dark") %>% 
  filter(!Species == "ASU") %>% 
  filter(!Species == "Sediment") %>% 
  mutate(Site = as.factor(Site),
         Date = as.factor(Date),
         rate_neg = rate*-1) 


# Metabolism Models -------------------------------------------------------

# Code for model comparisons 

#### NEP ---------------

Base_NEP <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(u_Light_Middle, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(ABR_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Depth_End, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = NEP
)

summary(Base_NEP)

# remove ABR

Mod1_NEP <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(u_Light_Middle, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Depth_End, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = NEP
)

summary(Mod1_NEP)

# remove Depth

Mod2_NEP <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(u_Light_Middle, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = NEP
)

summary(Mod2_NEP)

# remove sal

Mod3_NEP <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(u_Light_Middle, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = NEP
)

summary(Mod3_NEP)

Mod4_NEP <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(u_Light_Middle, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = NEP
)

summary(Mod4_NEP)


AIC(
  Base_NEP,
  Mod1_NEP,
  Mod2_NEP,
  Mod3_NEP
)

# Mod2 has lowest AIC

# Directional 

draw(Mod2_NEP, residuals = FALSE) 



### ER -----------------


Base_ER <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(ABR_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Depth_End, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = ER
)

summary(Base_ER)

# Remove ABR

Mod1_ER <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Depth_End, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = ER
)

summary(Mod1_ER)

# remove Depth

Mod2_ER <- bam(
  rate ~ s(Temp_u, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = ER
)

summary(Mod2_ER)


AIC(
  Base_ER,
  Mod1_ER,
  Mod2_ER
)

# Model 2

# Direction
draw(Mod3_ER, residuals = FALSE)



# Scale Up ending ---------------------------------------------------------

# SAV area taken from Hensel 2023's paper (analysis modeleled after Alvaro 2025)

# Create a data frame with SAV area by scenario
sav_area <- tribble(
  ~Scenario,     ~Species, ~SAV_Area_ha,
  "2022",        "Rm",     9778.34,
  "2022",        "Zm",     8104.41,
  "NR_2060",     "Rm",     13717.48,
  "NR_2060",     "Zm",     2525.20,
  "NO_NR_2060",  "Rm",     4851.73,
  "NO_NR_2060",  "Zm",     2094.89
) %>%
  mutate(SAV_Area_m2 = SAV_Area_ha * 10000)

# Get average GPP and R per species across spatial data set (June site's only)

Scale_Meto <-
  Site_Avg_Area %>% 
  left_join(
    Environment %>% 
      select(Date, Site)
  ) %>% 
  # subset to spatial dataset
  mutate(Date = as_date(Date)) %>% 
  filter(!(Site == "Dameron Marsh")) %>% 
  filter(!(Date < "2022-06-05" | Date > "2022-07-07")) %>% 
  filter(!as_date(Date) == "2022-07-05") %>% 
  filter(!as_date(Date) == "2022-06-20") %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  mutate(
    GPP = as.numeric(GPP),
    R   = as.numeric(R)
  ) %>% 
  reframe(
    GPP_u = mean(GPP_day, na.rm = TRUE),
    GPP_sd = sd(GPP_day, na.rm = TRUE),
    R_u = mean(R_day, na.rm = TRUE),
    R_sd = sd(R_day, na.rm = TRUE)
  )


# Now join and multiply

Scale_results <- 
  Scale_Meto %>% 
  left_join(sav_area, by = "Species") %>% 
  mutate(
    GPP_area = GPP_u*SAV_Area_m2,
    GPP_area_sd = GPP_sd*SAV_Area_m2,
    R_area = R_u*SAV_Area_m2,
    R_area_sd = R_sd*SAV_Area_m2) %>% 
  dplyr::select(-c(SAV_Area_m2)) %>% 
  # Getting back to hectares from m2
  mutate(across(GPP_area:R_area_sd, ~ .x / 10000)) %>% 
  # Getting from mmol O2 (or C) to molC
  mutate(across(GPP_area:R_area_sd, ~ .x / 1000))
# Now values area in mol C per hectare per day 


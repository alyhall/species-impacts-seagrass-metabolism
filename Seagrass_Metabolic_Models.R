## Header and Preamble ---------------------------
##
## Script name: Seagrass Metabolic Models
##

## Load up the packages needed: 
library(readxl)
library(ggpubr)
library(rstatix)
library(mgcv)
library(tidyverse)

##
## When using or reading this code, we recommended turning on code soft-wrap.
########End

## Data Load In ----------------

# Read in the compiled metabolism csv (Seagrass_Metabolism_Dataframe.csv)

# Here we name that dataframe as "MegaDS" in subsequent analyses

MegaDS <- 
  read.csv("Seagrass_Metabolism_Dataframe.csv")

# Seperate by metabolic process for each GAMM model and only retain
# seagrass ecosystem measurements:

# NEP

NEP <-
  MegaDS %>% 
  filter(Treatment == "Light") %>% 
  filter(!Species == "ASU") %>% 
  filter(!Species == "Sediment") %>% 
  mutate(
    Julian = yday(Date),
  ) %>% 
  arrange(Species, Julian) %>%
  group_by(Species) %>%
  mutate(
    Site = as.factor(Site),
    Date = as.factor(Date)
  ) %>%
  ungroup()

# merged_data is your combined dataset
NEP <- 
  NEP[
    , c(
      "rate", "Temp_u", "u_Light_Middle", 
      "pH", "Sal", "Turb",
      "Julian", "ABR_Bio", "A_Above_Bio", 
      "U_Leaf_Area", "Depth_End", "Species")
  ]


# ER

ER <-
  MegaDS %>% 
  filter(Treatment == "Dark") %>% 
  filter(!Species == "ASU") %>% 
  filter(!Species == "Sediment") %>% 
  mutate(
    Julian = yday(Date),
  ) %>% 
  arrange(Species, Julian) %>%
  group_by(Species) %>%
  mutate(
    Site = as.factor(Site),
    Date = as.factor(Date),
    rate_neg = rate*-1
  ) %>%
  ungroup()

ER <- ER[
  , c(
    "rate", "Temp_u", "pH", 
    "Sal", "Julian", "ABR_Bio", 
    "A_Above_Bio", "U_Leaf_Area", 
    "Depth_End", "Species")
]


# Model (GAMMS) Analysis ------------------------------------------------------------

# NEP

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


# ER


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

Mod3_ER <- bam(
  rate_neg ~ s(Temp_u, by = Species,  k = 8) + 
    s(pH, by = Species,  k = 8) + 
    s(Sal, by = Species,  k = 8) + 
    s(A_Above_Bio, by = Species,  k = 8) +
    s(U_Leaf_Area, by = Species,  k = 8) +
    s(Date, bs='re') +
    s(Site, bs='re'),
  data = ER
)

summary(Mod3_ER)


AIC(
  Base_ER,
  Mod1_ER,
  Mod2_ER,
  Mod3_ER
)
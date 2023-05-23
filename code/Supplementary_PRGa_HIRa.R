#This script generates the supplememntary Dose-response curves and EC50s for the PRGamides and HIRamides of the Thiel et al paper on Nematostella GPCR deorphanisation



rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# load packages -----------------------------------------------------------
{
library(tidyverse)
library(cowplot)
library(png)
library(patchwork)
library(drc)
library(usethis)
library(gitcreds)
library(devtools)
library(knitr)
library(vroom)
library(credentials)
library(plotly)
library(networkD3)
library(webshot2)
}

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
# read data ---------------------------------------------------------------
AllGPCR <- vroom("data/Nvec_PRGa_HIRa_full_subset.csv")

# BE CAREFUL OF THE HEADER OF YOUR CONCENTRATIONS, depending on how the values were saved this may or may not work------
#I solve it by copying the names on the headers of my file.
AllGPCRtoplot <- AllGPCR %>%
  pivot_longer(c(
    "0", "1.00E-13", "1.00E-12",	"1.00E-11", "1.00E-10",	
    "1.00E-09",	"1.00E-08",	"1.00E-07",	"1.00E-06",	
    "1.00E-05",	"1.00E-04"
    ),
    names_to = "concentration", values_to = "luminescence"
)


#convert conc values to double
AllGPCRtoplot <- AllGPCRtoplot %>% 
  mutate(concentration = as.double(concentration))

# delete the NA data from the row because otherwise normalization  --------

AllGPCRtoplot <- AllGPCRtoplot[!is.na(AllGPCRtoplot$luminescence), ]                 # Omit NA by column via is.na


#function to normalise to reference (zero ligand control)
normalize_to_ctr <- function(x) {
  return (100*(x - x[1]) / (max(x) - x[1]))
}

grids <- function(axis = c("xy", "x", "y"), color = "azure2", size = NULL, linetype = NULL)
{
  axis <- match.arg(axis)
  grid.major <- element_line(color = color, size = size,
                             linetype = linetype)
  grid.minor <- element_line(color = color, size = 0.25,
                             linetype = linetype)
  
  switch(axis,
         xy = theme(panel.grid.major = grid.major, panel.grid.minor = grid.minor),
         x = theme(panel.grid.major.x = grid.major, panel.grid.minor.x = grid.minor),
         y = theme(panel.grid.major.y = grid.major, panel.grid.minor.y = grid.minor)
  )
}

#This is normalization with all the different receptors
#accomplished by grouping by replicate, receptor and peptide will perform the individual
#normalization!
AllGPCRtoplot <- AllGPCRtoplot %>%
  group_by(Replicate, Receptor, Peptide)%>%
  mutate('norm_luminescence'=normalize_to_ctr(luminescence))


# plotting function -------------------------------------------------------

DRC_plot <- function(Receptor_name){
  AllGPCRtoplot %>% 
    filter(Receptor == Receptor_name) %>%
    ggplot(aes(
      x = concentration, 
      y = norm_luminescence, 
      group = Peptide, 
      colour = Peptide)
      ) +
    geom_boxplot(aes(
      group = paste(concentration, Peptide)), 
      outlier.shape=NA, size=0.5, width=0.2
      ) + 
    geom_smooth(
      method = drm, method.args = list(fct = L.4()), 
      se = FALSE, linewidth = 0.8
      ) +
    scale_x_log10(
      breaks = c(1e-13,1e-12,1e-10,1e-8,1e-6,1e-4), 
      limits = c(5e-14,1e-3)
      ) +
    scale_y_continuous(limits = c(-10, 110), breaks = c(0, 50, 100)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 12, angle = 90), 
          axis.text.y = element_text(size = 12), 
          legend.text = element_text(size=15), 
          legend.title = element_text(size=12),
          legend.key.size = unit(2.5, "mm"),
          legend.position = "bottom",
          legend.margin = margin(unit(-30,"mm")),
          axis.title=element_text(size=18), 
          axis.title.x=element_text(margin = margin(t = 8)),
          panel.background = element_blank(),
          plot.title = element_text(size = 16, hjust = 0.3)) +
    stat_summary(fun.y = mean, geom = "point", shape = 20, size = 1.5) +
    labs(
      x = "", y = "", 
      colour = "", title = Receptor_name
      ) +
    scale_color_manual(
      values = c(
        "grey40", "#E69F00", "#56B4E9", 
        "#009E73", "#F0E442", "#0072B2", 
        "#D55E00", "#CC79A7")
      )
  
  ggsave(
    paste("pictures/", Receptor_name, ".png", sep = ""), 
    width = 1000, height = 1000, limitsize = TRUE, 
    units = c("px"), bg='white'
  )
}

# plot and save all receptor plots ----------------------------------------
# Plot with fitted curve and LL4 analysis DRC
#plots out the dataset with the corresponding 4-parameter log-logit dose response curves

AllGPCRtoplot %>% 
  ungroup() %>%
  dplyr::select(Receptor) %>%
  unique() %>%
  pull(Receptor) %>%
  sapply(function(Receptor) DRC_plot(Receptor))

# EC50 and slope calculation --------------------------------------------------------

ListREceptors <- unique(AllGPCRtoplot$Receptor)
ListPeptides <- unique(AllGPCRtoplot$Peptide)
df <- ""
for (i in ListREceptors) {
  RecSubset <- filter(AllGPCRtoplot, Receptor == i)
  for (k in ListPeptides) {
    SubsetRecPep <- filter(RecSubset, Peptide == k)
    if (nrow(SubsetRecPep) > 0) {
      #      write.csv(SubsetRecPep, file = paste("Subs_", i, k), sep = ",", row.names = FALSE);
      model <- drm(
        norm_luminescence~concentration, 
        data = SubsetRecPep, 
        fct=LL.4(names =c("Slope", "Lower Limit", "Upper Limit", "ED50" ))
      )
      Slope <- model$coefficients[1]
      EC50 <- model$coefficients[4]
      dfrow <- c(Receptor = i, Peptide = k, EC50 = EC50, Slope = Slope)
      df <- rbind(df, dfrow)
      write.csv(df, file = "supplements/slope_and_EC50table_HIRaPRGa.csv", sep = ",", row.names = FALSE);
    }
  }
}


# save the table as supplementary file ------------------------------------

readr::write_csv(AllGPCRtoplot, file="supplements/Supplementary_table_HIRaPRGa.csv", na="", quote="none")

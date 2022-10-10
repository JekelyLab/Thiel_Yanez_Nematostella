rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# load packages -----------------------------------------------------------
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

#check working directory- should be the project dir, all directories will be defined relative to this dir
getwd()


# read data ---------------------------------------------------------------
AllGPCR <- vroom("/home/dt380/Dropbox/R_hub/Nvec_GPCRs/data/00_Nvec_Dose_response_assays_cleaned.csv")

# BE CAREFUL OF THE HEADER OF YOUR CONCENTRATIONS, depending on how the values were saved this may or may not work------
#I solve it by copying the names on the headers of my file.
AllGPCRtoplot <- AllGPCR %>%
  pivot_longer(c("0",	"1.00E-13", "1.00E-12",	"1.00E-11", "1.00E-10",	"1.00E-09",	"1.00E-08",	"1.00E-07",	"1.00E-06",	"1.00E-05",	"1.00E-04"), 
               names_to = "concentration", values_to = "luminescence")

#AllGPCRtoplot <- AllGPCR %>%
 # pivot_longer(starts_with(c("0", "1e")),
  #             names_to = "concentration", values_to = "luminescence")

#convert conc valus to double
AllGPCRtoplot$concentration <- as.double(AllGPCRtoplot$concentration)

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

# Plot with fitted curve and LL4 analysis DRC -----------------------------
#plots out the dataset with the corresponding 4-parameter log-logit dose response curves
AllGPCRtoplot %>% 
  ggplot(aes(x=concentration, y=norm_luminescence, colour=Peptide, group=Peptide)) +
  #Remove this comment and the next # to plot the boxplot. However, there is not
  #enough space in the case of multiple peptides in the same graph for all the boxes
  #geom_boxplot(aes(group=concentration, colour=Peptide, group=Peptide), outlier.shape=NA, size=0.15, width=0.2) +
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE) +
  scale_x_log10(breaks = c(1e-14,1e-13,1e-12,1e-10,1e-8,1e-6,1e-4),limits = c(1e-13,1e-3)) +
  theme_gray(base_size = 11) + 
  grids(linetype = "longdash") + #to change the colour if grids, change it in the function above +
  theme(axis.text=element_text(size = 14), 
        legend.text = element_text(size=16), 
        legend.title=element_text(size=20),
        axis.title=element_text(size=26), 
        axis.title.x=element_text(margin = margin(t = 10)),
        panel.background = element_blank()) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4) +
  stat_summary(fun.data = mean_se, geom="errorbar", size=0.3) +
  facet_wrap(vars(Receptor))




# EC50 calculation --------------------------------------------------------

ListREceptors <- unique(AllGPCRtoplot$Receptor)
ListPeptides <- unique(AllGPCRtoplot$Peptide)
df <- ""
for (i in ListREceptors) {
  RecSubset <- filter(AllGPCRtoplot, Receptor==i)
  for (k in ListPeptides) {
    SubsetRecPep <- filter(RecSubset, Peptide==k)
    if (nrow(SubsetRecPep)>0) {
      write.csv(SubsetRecPep, file = paste("Subs_", i, k), sep = ",", row.names = FALSE);
      model <- drm(norm_luminescence~concentration, data = SubsetRecPep, fct=LL.4(names =c("Slope", "Lower Limit", "Upper Limit", "ED50" )))
      EC50 <- model$coefficients[4]
      dfrow <- c(Receptor=i, Peptide=k, EC50=EC50)
      df <- rbind(df, dfrow)
      write.csv(df, file = "Ec50table.csv", sep = ",", row.names = FALSE);
    }
  }}

# Save the plot as pdf and png --------------------------------------------

ggsave("pictures/GPCR_curve.pdf", 
       width = 2000, 
       height = 1600, limitsize = TRUE, 
       units = c("px"))
ggsave("pictures/GPCR_curve.png", 
       width = 1700, 
       height = 1400, limitsize = TRUE, 
       units = c("px"), bg='white')




# save the table as supplementary file ------------------------------------

readr::write_csv(GPCR, file="supplements/Supplementary_table1.csv", na="", quote="none")



# assemble figure ---------------------------------------------------------

#define layout for patchwork to assemble figure panels
layout <-"
AABB"

Fig1 <- panelA + panelB + 
  plot_layout(design = layout, heights = c(1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

#save
ggsave("figures/Figure1.pdf", limitsize = FALSE, 
         units = c("px"), Fig1, width = 1600, height = 800)
ggsave("figures/Figure1.png", limitsize = FALSE, 
         units = c("px"), Fig1, width = 1600, height = 800, bg='white')




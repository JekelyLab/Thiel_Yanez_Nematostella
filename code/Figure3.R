{
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
library(plotly)
library(networkD3)
library(webshot2)
}

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
# read data ---------------------------------------------------------------
AllGPCR <- vroom("data/00_Nvec_Dose_response_assays_cleaned.csv")

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

DRC_plot("GLWLp.R018b")
# plot and save all receptor plots ----------------------------------------
# Plot with fitted curve and LL4 analysis DRC
#plots out the dataset with the corresponding 4-parameter log-logit dose response curves

AllGPCRtoplot %>% 
  ungroup() %>%
  dplyr::select(Receptor) %>%
  unique() %>%
  pull(Receptor) %>%
  sapply(function(Receptor) DRC_plot(Receptor))


# EC50 calculation --------------------------------------------------------

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
      EC50 <- model$coefficients[4]
      dfrow <- c(Receptor = i, Peptide = k, EC50 = EC50)
      df <- rbind(df, dfrow)
      write.csv(df, file = "supplements/EC50table.csv", sep = ",", row.names = FALSE);
    }
  }
}

# draw a table of EC50 values ---------------------------------------------

{
table <- plot_ly(
  type = 'table',
  columnwidth = c(
    3.8, 2.2, 2.5, 0.1,
    3.8, 2.2, 2.5, 0.1,
    3.8, 2.2, 2.5, 0.1,
    3.8, 2.2, 2.5
    ),
  columnorder = c(
    0, 1, 2, 3, 
    4, 5, 6, 7,
    8, 9, 10, 11,
    12, 13, 14
    ),
  header = list(
    values = c(
      "Receptor", "Peptide", "EC50", "", 
      "Receptor", "Peptide", "EC50", "", 
      "Receptor", "Peptide", "EC50", "", 
      "Receptor", "Peptide", "EC50"
      ),
    align = c("center"),
    line = list(width = 1, color = 'black'),
    fill = list(color = c(
      "#E69F00", "#F0E442", "#cccccc", "#cccccc",
      "#E69F00", "#F0E442", "#cccccc", "#cccccc",
      "#E69F00", "#F0E442", "#cccccc", "#cccccc",
      "#E69F00", "#F0E442", "#cccccc")),
    font = list(
      family = "Arial", size = 14, color = "black")
  ),
  cells = list(
    values = rbind(
      as_tibble(df[2:11, ]) %>%
        select(Receptor) %>%
        pull(), 
      as_tibble(df[2:11, ]) %>%
        select(Peptide) %>%
        pull(), 
      formatC(as.double(df[2:11, 3]), 
              format = "e", digits = 2),
      "",
      as_tibble(df[12:21, ]) %>%
        select(Receptor) %>%
        pull(), 
      as_tibble(df[12:21, ]) %>%
       select(Peptide) %>%
        pull(), 
      formatC(as.double(df[12:21, 3]), 
              format = "e", digits = 2),
      "",
      as_tibble(df[22:31, ]) %>%
        select(Receptor) %>%
        pull(), 
      as_tibble(df[22:31, ]) %>%
        select(Peptide) %>%
        pull(), 
      formatC(as.double(df[22:31, 3]), 
              format = "e", digits = 2),
      "",
      c(as_tibble(df[32:40, ]) %>%
          select(Receptor) %>%
          pull(), ""),
      c(as_tibble(df[32:40, ]) %>%
          select(Peptide) %>%
          pull(), ""),
      c(formatC(as.double(df[32:40, 3]), 
                format = "e", digits = 2), "")
    ),
    align = c("center"),
    line = list(color = "black", width = 0.3),
    font = list(family = "Arial", size = 12, 
                color = c("black"))
  )
)

table

saveNetwork(table, "pictures/EC50_table.html")
webshot2::webshot(url="pictures/EC50_table.html",
                  file="pictures/EC50_table.png",
                  vwidth=850, vheight=500, #define the size of the browser window
                  cliprect = c(58, 23, 784, 233), zoom=2)
}

#histogram of min EC50 values (for each receptor only the lowest value)
as_tibble(df[2:40, ]) %>%                      # Specify data frame
  group_by(Receptor) %>%
  rename(EC50 = starts_with("EC50")) %>%
  mutate(EC50 = as.double(EC50)) %>%
  summarise_at(vars(EC50),  # Specify column
               list(minEC50 = min)) %>%     #select only the min value for each receptor
  ggplot(aes(minEC50)) +
  geom_histogram(color = "grey50", fill = Okabe_Ito[1]) +
  scale_x_log10(breaks = c(1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)) +
  theme_minimal() +
  labs(x = bquote(EC[50])) +
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
        plot.title = element_text(size = 16, hjust = 0.3)) 

ggsave(
  paste("pictures/EC50_histogram.png"), 
  width = 1000, height = 1000, limitsize = TRUE, 
  units = c("px"), bg='white'
)

# save the table as supplementary file ------------------------------------

readr::write_csv(AllGPCRtoplot, file="supplements/Supplementary_table1.csv", na="", quote="none")


# assemble figure ---------------------------------------------------------

#read panels

{
GLWLp.R018a <- ggdraw() + draw_image(readPNG("pictures/GLWLp.R018b.png")) + 
  draw_label("norm. luminescence", x = 0.05, y = 0.5, size = 10, angle = 90)
GLWLp.R018b <- ggdraw() + draw_image(readPNG("pictures/GLWLp.R018b.png"))
HIRa.R021 <- ggdraw() + draw_image(readPNG("pictures/HIRa.R021.png"))
HIRa.R029 <- ggdraw() + draw_image(readPNG("pictures/HIRa.R029.png")) + 
  draw_label("norm. luminescence", x = 0.05, y = 0.5, size = 10, angle = 90)
FLRNa.R026 <- ggdraw() + draw_image(readPNG("pictures/FLRNa.R026.png"))
FLRNa.R197 <- ggdraw() + draw_image(readPNG("pictures/FLRNa.R197.png"))
FLRNa.R230 <- ggdraw() + draw_image(readPNG("pictures/FLRNa.R230.png"))
PFHa.R036 <- ggdraw() + draw_image(readPNG("pictures/PFHa.R036.png"))
QWa.R069 <- ggdraw() + draw_image(readPNG("pictures/QWa.R069.png"))
QGRFa.R070 <- ggdraw() + draw_image(readPNG("pictures/QGRFa.R070.png")) + 
  draw_label("norm. luminescence", x = 0.05, y = 0.5, size = 10, angle = 90)
QGRFa.R234 <- ggdraw() + draw_image(readPNG("pictures/QGRFa.R234.png"))
QITRFa.R196 <- ggdraw() + draw_image(readPNG("pictures/QITRFa.R196.png"))
LRWa1.R019 <- ggdraw() + draw_image(readPNG("pictures/LRWa1.R019.png"))
LRWa.R193 <- ggdraw() + draw_image(readPNG("pictures/LRWa.R193.png"))
LRWa3.R204 <- ggdraw() + draw_image(readPNG("pictures/LRWa3.R204.png"))
LRWa2.R213  <- ggdraw() + draw_image(readPNG("pictures/LRWa2.R213.png")) + 
  draw_label("norm. luminescence", x = 0.05, y = 0.5, size = 10, angle = 90)
PRGa.R028 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R028.png"))
PRGa.R032 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R032.png"))
PRGa.R198 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R198.png"))
PRGa.R199 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R199.png"))
PRGa.R200 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R200.png"))
PRGa.R202 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R202.png")) + 
  draw_label("norm. luminescence", x = 0.05, y = 0.5, size = 10, angle = 90)
PRGa.R210 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R210.png"))
PRGa.R211 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R211.png"))
PRGa.R219 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R219.png"))
PRGa.R220 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R220.png"))
PRGa.R221 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R221.png"))
PRGa.R222 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R222.png")) + 
  draw_label("norm. luminescence", x = 0.05, y = 0.5, size = 10, angle = 90)
PRGa.R223 <- ggdraw() + draw_image(readPNG("pictures/PRGa.R223.png"))

EC50_table <- ggdraw() + draw_image(readPNG("pictures/EC50_table.png"))
EC50_hist <- ggdraw() + draw_image(readPNG("pictures/EC50_histogram.png"))


#define layout for patchwork to assemble figure panels
layout <- "
##AbBc
######
CdDeEf
######
FgGhHi
######
IjJkKl
######
LmMnNo
######
OppppP
"

Fig3 <- GLWLp.R018a + GLWLp.R018b + HIRa.R021 + 
  HIRa.R029 + FLRNa.R026 + FLRNa.R197 + FLRNa.R230 + PFHa.R036 + QWa.R069 + 
  QGRFa.R070 + QGRFa.R234 + QITRFa.R196 + LRWa1.R019  + 
  LRWa.R193 + LRWa3.R204 + LRWa2.R213 + PRGa.R028 + 
  PRGa.R032 + PRGa.R198 + PRGa.R199 + PRGa.R200 + 
  PRGa.R202 + PRGa.R210 + PRGa.R211 + PRGa.R219 + 
  PRGa.R220 + PRGa.R221 + PRGa.R222 + PRGa.R223 +
  EC50_table + EC50_hist +
  plot_layout(design = layout, heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1)) +
  plot_annotation(tag_levels = "i") & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("figures/Figure3.pdf", limitsize = FALSE, 
         units = c("px"), Fig3, width = 3400, height = 3700)

ggsave("figures/Figure3.png", limitsize = FALSE, 
         units = c("px"), Fig3, width = 3400, height = 3700, bg='white')


}

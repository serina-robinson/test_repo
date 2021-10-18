# Install packages
pacman::p_load("tidyverse", "readxl", "RColorBrewer")

# Read in the dataset
file.list <- list.files("data/05122020_bioactivity_assay/", pattern='*.xlsx', full.names = T)
file.list <- setNames(file.list, file.list) 
file.list
# Combine
df <- map_df(file.list, read_excel, .id = "id") %>%
  janitor::clean_names() %>%
  dplyr::select(-temperature_c) %>%
  dplyr::mutate(time_pt = as.numeric(gsub("\\h.xlsx", "", word(id, sep = "_", -1)))) %>%
  dplyr::mutate(organism = word(word(id, sep = "//", 2), sep = "_", 1)) %>%
  dplyr::mutate(row_lett = rep(toupper(letters[1:8]), nrow(.)/8)) %>%
  reshape2::melt(., id.vars = c("id", "organism", "time_pt", "row_lett")) %>%
  dplyr::mutate(well_position = paste0(row_lett, word(variable, sep = "x", 2))) %>%
  dplyr::filter(well_position %in% c("B2", "B3", "C4", "C5", "D6", "D7",
                                     "E8", "E9", "F10", "F11", "G2", "G3")) %>%
  dplyr::mutate(rep_number = case_when(well_position %in% c("B2", "C4", "D6", "E8", "F10") ~ "rep1",
                                       well_position %in% c("B3", "C5", "D7", "E9", "F11") ~ "rep2",
                                       well_position %in% c("G2") ~ "rep3",
                                       well_position %in% c("G3") ~ "rep4")) %>%
  dplyr::mutate(treatment = case_when(grepl("B|G", well_position) ~ "solvent_control",
                                      grepl("C", well_position) ~ "50 µM ampicillin",
                                      grepl("D", well_position) ~ "50 µM peptide",
                                      grepl("E", well_position) ~ "25 µM peptide",
                                      grepl("F", well_position) ~ "50 µM chloramphenicol"))
  

# Set random number seed for random color palette
set.seed(1234)
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c(pal[c(1,2, 4, 3,5, 8)], "dodgerblue", "goldenrod", "chartreuse",  "blue", "gold1", "black")

pdf("output/05122020_bioactivity_results.pdf")
ggplot(df, aes(x = time_pt, y = value, group = rep_number, color = rep_number)) +
  geom_point() + 
  geom_line() +
  labs(y = "Optical density (600nm)", x = "Time (hrs)") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 10),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="top") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) +
  facet_grid(rows = vars(organism),
             cols = vars(treatment),
             scales = "free")
dev.off()


pdf("output/05122020_log10_bioactivity_results.pdf")
ggplot(df, aes(x = time_pt, y = value, group = rep_number, color = rep_number)) +
  geom_point() + 
  geom_line() +
  scale_y_log10(breaks = scales::pretty_breaks(n = 3)) +
  labs(y = "Optical density (600nm)", x = "Time (hrs)") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 10),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="top") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) +
  facet_grid(rows = vars(organism),
             cols = vars(treatment),
             scales = "free")
dev.off()


### Now filter out to only keep the peptides and the controls
df_filt$treatment
df_filt <- df %>%
  dplyr::filter(grepl("peptide|solvent", treatment)) #%>%
 # dplyr::filter(time_pt != 48)
df_filt

pal3 <- c(pal[c(1,2)], "gray60")
pdf("output/05122020_overlay_results_scales_fixed.pdf", width = 4, height = 10)
ggplot(df_filt, aes(x = time_pt, y = value, group = treatment, color = treatment)) +
  geom_point(alpha = 0.2) + 
  scale_y_log10(breaks = scales::pretty_breaks(n = 3)) +
  labs(y = "Optical density (600nm)", x = "Time (hrs)") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 10),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="right") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal3) +
  facet_wrap(vars(organism), nrow = 7, ncol = 1)
dev.off()

pdf("output/05122020_overlay_results_scales_free.pdf", width = 4, height = 10)
ggplot(df_filt, aes(x = time_pt, y = value, group = treatment, color = treatment)) +
  geom_point(alpha = 0.2) + 
  scale_y_log10(breaks = scales::pretty_breaks(n = 3)) +
  labs(y = "Optical density (600nm)", x = "Time (hrs)") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 10),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="right") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal3) +
  facet_wrap(vars(organism), nrow = 7, ncol = 1, scales = "free")
dev.off()

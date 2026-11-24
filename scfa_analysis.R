library("tidyverse")
library('dplyr')
load("/Users/mattigo/Documents/CU_Folders/Projects_Ongoing/NFLPA_Microbiome/Test_folder_2/NFLPA_LIBRA_US_VMP_Microbiome.RData")

scfa_list = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Eubacterium rectale", "Eubacterium hallii",
              "Subdoligranulum variabile", "Roseburia hominis", "Anaerobutyricum soehngenii", "Coprococcus comes",
              "Coprococcus catus", "Butyricicoccus pullicaecorum", "Agathobacter rectalis", "Bacteroides thetaiotaomicron",
              "Bacteroides vulgatus", "Bacteroides uniformis", "Bacteroides fragilis", "Parabacteroides distasonis",
              "Megasphaera elsdenii", "Veillonella parvula", "Veillonella dispar", "Ruminococcus obeum", "Blautia obeum",
              "Bifidobacterium longum", "Bifidobacterium adolescentis", "Bifidobacterium breve", "Lactobacillus rhamnosus",
              "Lactobacillus plantarum", "Akkermansia muciniphila", "Ruminococcus gnavus")



relative_taxa_adbund_filtered = as.data.frame(input_data_rel)
relative_taxa_adbund_filtered = taxa_species_relative #%>% dplyrfilter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")
relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]
relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")
relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance")
meta$SampleID = rownames(meta)
all_taxa_meta_data = plyr::join(relative_taxa_adbund_filtered_update3_melt, meta, by="SampleID") %>%
  mutate(Species = str_remove(Species, "s__")) %>%
  mutate(Species = str_replace(Species, "_", " "))
all_taxa_meta_data2 = all_taxa_meta_data %>% mutate(Abundance_log = log10(Abundance + 1E-10))
all_taxa_meta_data2$LIBRA.score.quartile = factor(all_taxa_meta_data2$LIBRA.score.quartile, levels = c("4", "3", "2", "1"))
all_taxa_meta_data3 = all_taxa_meta_data2 %>% arrange(Species, (LIBRA.score)) %>%
  mutate(Libra.score.half = case_when(LIBRA.score.quartile == "1" | LIBRA.score.quartile == "2" ~ "Lower",
                                      T~"Upper"))
all_taxa_meta_data_scfa = all_taxa_meta_data3 %>% dplyr::filter(Species %in% scfa_list)
all_taxa_meta_data_scfa = data.frame(all_taxa_meta_data_scfa)
all_taxa_meta_data_update_scfa = all_taxa_meta_data_scfa %>% group_by(id) %>% mutate(scfa_abund = sum(Abundance)) %>% ungroup()
all_taxa_meta_data_update_scfa$scfa_abund
all_taxa_meta_data_update_scfa_2 = all_taxa_meta_data_update_scfa %>% dplyr::select(SampleID, scfa_abund, LIBRA.score, LIBRA.score.quartile, Libra.score.half, Age, Homeless.ever,BMI, Race, Ethnicity,Sex) %>% distinct()

range(all_taxa_meta_data_update_scfa_2$scfa_abund, na.rm = TRUE)

mean(all_taxa_meta_data_update_scfa_2$scfa_abund)
sd(all_taxa_meta_data_update_scfa_2$scfa_abund)

p1a_scfa <- ggplot(data = all_taxa_meta_data_update_scfa_2, aes(x = LIBRA.score, y = scfa_abund)) +
  geom_point(size=3, aes(color=Sex)) +
  geom_smooth(se = F, method = "lm", linewidth=1.5) +
  #stat_cor(label.y = 0.65, label.x = 6, size = 6, fontface="bold")+
  theme_bw() +
  labs(tag = "A", x="LIBRA Score", y="Relative Abundance of SCFA Bacteria") +
  theme(text = element_text(size = 12, face = "bold"))#,
        #legend.position = "none")
p1a_scfa

all_taxa_meta_data_update_scfa_2$LIBRA.score.quartile = factor(all_taxa_meta_data_update_scfa_2$LIBRA.score.quartile, levels = c("1","2","3","4"))
p1b_scfa <- ggplot(data = all_taxa_meta_data_update_scfa_2, aes(x = LIBRA.score.quartile, y = scfa_abund)) +
  geom_boxplot(aes(fill = LIBRA.score.quartile)) + 
  geom_jitter(aes(color = LIBRA.score.quartile)) + 
  scale_color_manual(values=(colvec)) + 
  scale_fill_manual(values=(colvec))+
  theme_bw() + 
  labs(tag = "B", x="LIBRA Score Quartile", y="Relative Abundance of SCFA Bacteria") +
  stat_compare_means(method = "kruskal.test", label.y.npc="top", label.x = 2.75, size = 5) + # Add Kruskal-Wallis p-value
  theme(plot.tag = element_text(size = 12, face = "bold"),
        text = element_text(size = 12, face = "bold"),
        legend.position = "none") 
p1b_scfa

summary(lm(scfa_abund~LIBRA.score, data = all_taxa_meta_data_update_scfa_2))

summary(lm(scfa_abund~LIBRA.score*Sex+Homeless.ever +Age + BMI, data = all_taxa_meta_data_update_scfa_2))
car::Anova(lm(scfa_abund~ LIBRA.score*Sex+Homeless.ever + Age + BMI, data = all_taxa_meta_data_update_scfa_2))
car::vif(lm(scfa_abund~LIBRA.score.quartile*Sex + Age + BMI + Homeless.ever + Race + Ethnicity, data = all_taxa_meta_data_update_scfa_2))
all_taxa_meta_data_update_scfa_2


p_scfa = ggarrange(p1a_scfa,p1b_scfa, ncol=2,nrow=1)
ggsave(paste0("/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/NFLPA_Microbiome/Analysis/Figures/MERIT_SCFA_by_LIBRA_grant_fig_", Sys.Date(),".jpeg"),p_scfa, height = 5, width = 14, units="in" )


##########PART2

load("/Users/mattigo/Documents/CU_Folders/Projects_Ongoing/NFLPA_Microbiome/Test_folder_2/NFLPA_LIBRA_US_VMP_Microbiome.RData")

####update SCFA bacteria
scfa_list = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Eubacterium rectale", "Eubacterium hallii",
              "Subdoligranulum variabile", "Roseburia hominis", "Anaerobutyricum soehngenii", "Coprococcus comes",
              "Coprococcus catus", "Butyricicoccus pullicaecorum", "Agathobacter rectalis", "Bacteroides thetaiotaomicron",
              "Bacteroides vulgatus", "Bacteroides uniformis", "Bacteroides fragilis", "Parabacteroides distasonis",
              "Megasphaera elsdenii", "Veillonella parvula", "Veillonella dispar", "Ruminococcus obeum", "Blautia obeum",
              "Bifidobacterium longum", "Bifidobacterium adolescentis", "Bifidobacterium breve", "Lactobacillus rhamnosus",
              "Lactobacillus plantarum", "Akkermansia muciniphila", "Ruminococcus gnavus")

meta

#all_taxa_meta_data_update_scfa = all_taxa_meta_data_update %>% filter(Species %in% scfa_list) %>% group_by(Participant, Time_num) %>% mutate(scfa_abund = sum(total_abd))

#all_taxa_meta_data_update_scfa_2 = all_taxa_meta_data_update_scfa %>% select(Participant, age, sex, Ever_Homeless, Time_num, scfa_abund, LIBRA_Half,Total_LIBRA_Score, BMI_Current) %>% distinct()



relative_taxa_adbund_filtered = taxa_species_relative #%>% dplyr::filter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")

relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]

relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")


relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance") 

meta$SampleID = rownames(meta)
all_taxa_meta_data = plyr::join(relative_taxa_adbund_filtered_update3_melt, meta, by="SampleID") %>% 
  mutate(Species = str_remove(Species, "s__")) %>% 
  mutate(Species = str_replace(Species, "_", " "))

all_taxa_meta_data2 = all_taxa_meta_data %>% mutate(Abundance_log = log10(Abundance + 1E-10))
all_taxa_meta_data2$LIBRA.score.quartile = factor(all_taxa_meta_data2$LIBRA.score.quartile, levels = c("4", "3", "2", "1"))
all_taxa_meta_data3 = all_taxa_meta_data2 %>% arrange(Species, (LIBRA.score)) %>%
  mutate(Libra.score.half = case_when(LIBRA.score.quartile == "1" | LIBRA.score.quartile == "2" ~ "Lower",
                                      T~"Upper"))


all_taxa_meta_data_scfa = all_taxa_meta_data3 %>% dplyr::filter(Species %in% scfa_list)

all_taxa_meta_data_update_scfa = all_taxa_meta_data_scfa  %>% group_by(SampleID) %>% mutate(scfa_abund = sum(Abundance))
all_taxa_meta_data_update_scfa_2 = all_taxa_meta_data_update_scfa %>% select(SampleID, scfa_abund, LIBRA.score, LIBRA.score.quartile, Libra.score.half, Age, Homeless.ever,BMI) %>% distinct()

colvec <- rev(c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A"))

p1scfa <- ggplot(data = all_taxa_meta_data_update_scfa_2 , aes(x = factor(LIBRA.score.quartile, levels=c("1","2","3","4")), y = scfa_abund)) +
  geom_boxplot(aes(fill = LIBRA.score.quartile)) + 
  geom_jitter(aes(color = LIBRA.score.quartile)) + 
  scale_color_manual(values=(colvec)) + 
  scale_fill_manual(values=(colvec))+
  theme_bw() + 
  labs(tag = "B", x="LIBRA Score Quartile", y="Relative Abundance of SCFA Bacteria") +
  stat_compare_means(method = "kruskal.test", label.y.npc="top", label.x = 2.75, size = 3) + # Add Kruskal-Wallis p-value
  theme(plot.tag = element_text(size = 12, face = "bold"),
        legend.position = "none") 

p1scfa <- ggplot(data = all_taxa_meta_data_update_scfa_2 , aes(x = (Libra.score.half), y = scfa_abund)) +
  geom_boxplot(aes(fill = Libra.score.half)) + 
  geom_jitter(aes(color = Libra.score.half)) + 
  #scale_color_manual(values=(colvec)) + 
  #scale_fill_manual(values=(colvec))+
  theme_bw() + 
  labs(tag = "B", x="LIBRA Score Quartile", y="Relative Abundance of SCFA Bacteria") +
  stat_compare_means(method = "kruskal.test", label.y.npc="top", label.x = 1, size = 3) + # Add Kruskal-Wallis p-value
  theme(plot.tag = element_text(size = 12, face = "bold"),
        legend.position = "none") 

all_taxa_meta_data_update_scfa_2$Homeless.ever = 
model = (lm(scfa_abund ~ LIBRA.score + Age + BMI + Homeless.ever, data = all_taxa_meta_data_update_scfa_2))
summary(model)
TukeyHSD(model)

p1a_scfa <- ggplot(data = all_taxa_meta_data_update_scfa_2, aes(x = LIBRA.score, y = scfa_abund)) +
  geom_point() + 
  geom_smooth(se = F, method = "lm") +
  stat_cor(label.y.npc="top", label.x = 5.25, size = 3)+
  theme_bw() + 
  labs(tag = "A", x="LIBRA Score", y="Relative Abundance of SCFA Bacteria") +
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave(paste0("Figures/Figure_1b_Alpha_Boxplots", Sys.Date(),".jpeg"),p1b, height = 10, width = 14, units="in" )


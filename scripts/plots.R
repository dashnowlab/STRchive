library('plotly')
library('ggplot2')
library('cowplot')
library('jsonlite')
theme_set(theme_cowplot())
options(stringsAsFactors = FALSE)

disease.loci = fromJSON("data/STRchive-loci.json")

# Switch between units
disease.loci$benign_min_bp = disease.loci$benign_min * disease.loci$motif_len
disease.loci$benign_max_bp = disease.loci$benign_max * disease.loci$motif_len
disease.loci$int_min_bp = disease.loci$intermediate_min * disease.loci$motif_len
disease.loci$int_max_bp = disease.loci$intermediate_max * disease.loci$motif_len
disease.loci$path_min_bp = disease.loci$pathogenic_min * disease.loci$motif_len
disease.loci$path_max_bp = disease.loci$pathogenic_max * disease.loci$motif_len

disease.loci$benign_min_bp[disease.loci$benign_min_bp == 0] = NA
disease.loci$benign_max_bp[disease.loci$benign_max_bp == 0] = NA

disease.loci$disease_id = reorder(disease.loci$disease_id, -disease.loci$path_min_bp)

disease.loci = subset(disease.loci, !grepl("conflicting evidence", disease.loci$details))

# Allele size
p_size = ggplot(disease.loci, aes(x = disease_id)) +
  geom_linerange(aes(ymin = benign_max_bp, ymax = path_max_bp), 
                 linewidth = 0.5, linetype = "dotted", alpha = 0.3, color = "gray") +
  geom_linerange(aes(ymin = path_min_bp, ymax = path_max_bp, color = 'Pathogenic'
  ), linewidth = 2) +
  geom_linerange(aes(ymin = int_min_bp, ymax = int_max_bp, color = 'Intermediate'
  ), linewidth = 2) +
  geom_linerange(aes(ymin = benign_min_bp, ymax = benign_max_bp, color = 'Benign'
  ), linewidth = 2) +
  geom_point(data = subset(disease.loci, path_min_bp == path_max_bp), aes(y = path_min_bp, color = 'Pathogenic'), shape = 16, size = 2) +
  geom_point(data = subset(disease.loci, int_min_bp == int_max_bp), aes(y = int_min_bp, color= 'Intermediate'), shape = 16, size = 2) +
  geom_point(data = subset(disease.loci, benign_min_bp == benign_max_bp), aes(y = benign_min_bp, color= 'Benign'), shape = 16, size = 2) +
  scale_y_continuous(name = 'Allele size in base pairs', trans = 'log10') +
  scale_x_discrete(name = 'Disease') +
  scale_colour_manual(values = c('Benign' = '#00AFBB', 'Intermediate' = '#E7B800', 'Pathogenic' = '#FC4E07'), 
                      breaks = c('Benign', 'Intermediate', 'Pathogenic'),
                      name = 'Allele size') +
  theme(panel.grid.major.x = element_line(color = 'lightgrey', linetype = 'longdash')) +
  coord_flip()

htmltools::save_html(ggplotly(p_size, height = 1200), "data/plots/plotly_path_size.html")

disease.loci$inheritance = factor(disease.loci$inheritance, levels = c("AD/AR", "AD", "AR", "XD", "XR"))

# Age of onset
p_age = ggplot(subset(disease.loci, !is.na(disease.loci$age_onset_min) & disease.loci$inheritance != '') , 
               aes(x = reorder(disease_id, -age_onset_min), color = inheritance)) +
  geom_linerange(aes(ymin = age_onset_min, ymax = age_onset_max), linewidth = 1) +
  geom_point(aes(y = age_onset_min), size = 0.5) +
  geom_point(aes(y = age_onset_max), size = 0.5) +
  geom_linerange(aes(ymin = typ_age_onset_min, ymax = typ_age_onset_max,
  ), linewidth = 2) +
  scale_y_continuous(name = 'Age of onset (years)', 
                     breaks = c(seq(0, 100, 25), 18)) +
  scale_x_discrete(name = 'Disease') +
  scale_color_brewer(palette = 'Paired', direction = -1) +
  geom_segment(aes(x = 1, y = 18, xend = 66, yend = 18), linetype = 'longdash', color = 'lightgrey') +
  coord_flip()

htmltools::save_html(ggplotly(p_age, height = 1200), "data/plots/plotly_age_onset.html")

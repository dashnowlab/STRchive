library('plotly')
library('ggplot2')
library('cowplot')
theme_set(theme_cowplot())
options(stringsAsFactors = FALSE)

# Reference genome: HG38

# Pre-R data clean-up
# ```{bash, eval=FALSE}
# cut -f 5 str-disease-loci-2022.tsv > repeats.txt
# python normalise_str.py > repeats_norm.txt

disease.loci = read.csv('STR-disease-loci.csv')

# Switch between units
disease.loci$norm_min_bp = disease.loci$normal_min * disease.loci$repeatunitlen
disease.loci$norm_max_bp = disease.loci$normal_max * disease.loci$repeatunitlen
disease.loci$int_min_bp = disease.loci$intermediate_min * disease.loci$repeatunitlen
disease.loci$int_max_bp = disease.loci$intermediate_max * disease.loci$repeatunitlen
disease.loci$path_min_bp = disease.loci$pathogenic_min * disease.loci$repeatunitlen
disease.loci$path_max_bp = disease.loci$pathogenic_max * disease.loci$repeatunitlen

# Allele size
p_size = ggplot(disease.loci, aes(x = disease_id)) +
  geom_linerange(aes(ymin = norm_min_bp, ymax = norm_max_bp + 1, color = 'Normal'
                     ), linewidth = 1.5) +
  geom_linerange(aes(ymin = int_min_bp, ymax = int_max_bp + 1, color = 'Intermediate*'
                     ), linewidth = 1.5) +
  geom_linerange(aes(ymin = path_min_bp, ymax = path_max_bp + 1, color = 'Pathogenic'
                     ), linewidth = 1.5) +
  scale_y_continuous(name = 'Allele size in base pairs', trans = 'log10') +
  scale_x_discrete(name = 'Disease') +
  scale_colour_brewer(palette = 'Accent', breaks = c('Normal', 'Intermediate*', 'Pathogenic'), name = 'Allele size') + 
  theme(panel.grid.major.x = element_line(color = 'lightgrey', linetype = 'longdash')) +
  coord_flip()

htmltools::save_html(ggplotly(p_size, height = 1000), "images/plotly_path_size.html")

# Age of onset
p_age = ggplot(subset(disease.loci, !is.na(disease.loci$age_onset_min) & disease.loci$Inheritance != '') , aes(x = reorder(disease_id, -age_onset_min), color = Inheritance)) +
  geom_linerange(aes(ymin = age_onset_min, ymax = age_onset_max + 1,
                     ), size = 1.5) +
  scale_y_continuous(name = 'Age of onset (years)') +
  scale_x_discrete(name = 'Disease') +
  geom_segment(aes(x = 1, y = 18, xend = 55, yend = 18), linetype = 'longdash', color = 'lightgrey') +
  coord_flip()

htmltools::save_html(ggplotly(p_age, height = 1000), "images/plotly_age_onset.html")

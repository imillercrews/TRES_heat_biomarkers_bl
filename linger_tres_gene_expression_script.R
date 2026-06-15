#### Linger TRES project 
### Analysis of gene expression
# R 4.3.1 
 
## setwd
setwd("/N/project/Multibrain_IMC/linger_TRES/")

#### load libraries ####
## install packages
# BiocManager::install("DESeq2")
# BiocManager::install("vsn")
# BiocManager::install("ComplexHeatmap")
# devtools::install_github("RRHO2/RRHO2", build_opts = c("--no-resave-data", "--no-manual"))

## load libraries
library(tidyverse)
library("DESeq2")

# GO annotation
library(AnnotationForge)

# DEG
library(pvclust)
library(vsn)

# variance
library(variancePartition)

# interaction
library(ComplexHeatmap)
library(RRHO2)

# WGCNA
library(WGCNA)

# GO term analysis
library(clusterProfiler)
library(org.Hrustica.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(GO.db)


#### create GO database ####
### create GO database for barn swallow genes
## use Annotation forge
makeOrgPackageFromNCBI(
  version = "0.1",
  author = "Isaac Miller-Crews <imillerc@iu.edu>",
  maintainer = "Isaac Miller-Crews <imillerc@iu.edu>",
  outputDir = "/N/project/Snseq_IMC/Bird_hypo_snseq/data/references/updated_ref/tree_swallow/", 
  tax_id = "43150",               # Replace with your NCBI Taxonomy ID
  genus = "Hirundo",         # Replace with your Genus
  species = "rustica"           # Replace with your Species
)

# install database 
install.packages("/N/project/Snseq_IMC/Bird_hypo_snseq/data/references/updated_ref/tree_swallow/org.Hrustica.eg.db/", 
                 repos = NULL,
                 type = "source")

# load database
library(org.Hrustica.eg.db)

#### QC of mapping ####
### linger
linger.map.stats = read_tsv('data/Linger_rnaseq/multiqc_data_1/star_summary_table.txt') %>% 
  filter(Sample != '.')

## get summary stats
linger.map.stats %>% 
  filter(Sample != 'GSF3634-2881-90843_BL_S13') %>% 
  dplyr::select(`Total reads`,
                `Uniq aligned`) %>% 
  psych::describe()

### during
during.map.stats = read_tsv('data/During_rnaseq/multiqc_data/star_summary_table.txt') %>% 
  filter(Sample != '.')

## get summary stats
during.map.stats %>% 
  dplyr::select(`Total reads`,
                `Uniq aligned`) %>% 
  psych::describe()


#### process gene expression files ####
### load gene counts
## get list of gene count files
gene.count.files = list.files(path = 'data/Linger_rnaseq/bams/',
                              pattern = "ReadsPerGene.out.tab", 
                              full.names = TRUE)

# get sample names
sample.names = data.frame(Sample = gsub("_ReadsPerGene.out.tab",
                    "", 
                    basename(gene.count.files))) %>% 
  separate_wider_delim('Sample',
                       delim = '_BL',
                       names = c('Sample',
                                 NA)) %>% 
  separate_wider_delim('Sample',
                       delim = 'GSF3634-',
                       names = c(NA,
                                 'Sample')) %>% 
  mutate(Sample = str_replace_all(Sample,
                                  '-',
                                  '_')) %>% 
  mutate(Sample = paste0('Sample_',
                         Sample,
                         '_BL')) 
 

## get gene count data
# create dummy 
data.linger.gene.wide = data.frame(Gene = as.character())

# loop through all files
for (i in 1:length(gene.count.files)) {
  # get sample name
  tmp.sample = sample.names$Sample[i]
  
  # load tmp file
  # just keep gene and reverse strand counts
  tmp = read.table(gene.count.files[i], 
                         skip = 4) %>% 
    transmute(Gene = V1,
              V4 = V4) %>% 
    dplyr::rename(!!tmp.sample := 'V4')
  
  # combine with all data
  data.linger.gene.wide = data.linger.gene.wide %>% 
    full_join(tmp)
}

# remove temp files
rm(tmp)


# fix gene name
data.linger.gene.wide = data.linger.gene.wide %>% 
  separate_wider_delim('Gene',
                       delim = 'gene-',
                       names = c(NA,
                                 'Gene')) 


# save file
write.csv(data.linger.gene.wide,
          'data/data.linger.gene.wide.csv',
          row.names = FALSE)




#### During: process gene expression files ####
### load gene counts
## get list of gene count files
gene.count.files.during = list.files(path = 'data/During_rnaseq/bams/',
                                     pattern = "ReadsPerGene.out.tab", 
                                     full.names = TRUE)

# get sample names
sample.names.during = data.frame(Sample = gsub("_ReadsPerGene.out.tab",
                                               "", 
                                               basename(gene.count.files.during))) %>% 
  separate_wider_delim('Sample',
                       delim = '-BL',
                       names = c('Sample',
                                 NA)) %>% 
  separate_wider_delim('Sample',
                       delim = '-GSF2996-',
                       names = c(NA,
                                 'Sample')) %>% 
  separate_wider_delim('Sample',
                       delim = '21',
                       names = c(NA,
                                 'Sample')) %>% 
  mutate(Sample = paste0(Sample,
                         '_21_BL'))


## get gene count data
# create dummy 
data.during.gene.wide = data.frame(Gene = as.character())

# loop through all files
for (i in 1:length(gene.count.files.during)) {
  # get sample name
  tmp.sample = sample.names.during$Sample[i]
  
  # load tmp file
  # just keep gene and reverse strand counts
  tmp = read.table(gene.count.files.during[i], 
                   skip = 4) %>% 
    transmute(Gene = V1,
              V4 = V4) %>% 
    dplyr::rename(!!tmp.sample := 'V4')
  
  # combine with all data
  data.during.gene.wide = data.during.gene.wide %>% 
    full_join(tmp)
}

# remove temp files
rm(tmp)


# fix gene name
data.during.gene.wide = data.during.gene.wide %>% 
  separate_wider_delim('Gene',
                       delim = 'gene-',
                       names = c(NA,
                                 'Gene')) 


# save file
write.csv(data.during.gene.wide,
          'data/data.during.gene.wide.csv',
          row.names = FALSE)




#### Linger: load data ####
### gene expression data
data.linger.gene.wide = read.csv('data/data.linger.gene.wide.csv')


# pivot longer
data.linger.gene = data.linger.gene.wide %>% 
  pivot_longer(cols = -c(Gene),
               names_to = 'Sample',
               values_to = 'Counts')

# create wide format for deseq2 
data.linger.gene.mat = data.linger.gene.wide %>% 
  column_to_rownames('Gene')

### load sample data
data.linger.sample = read.csv('data/treeSwallowCounts.samples.linger.csv') %>% 
  dplyr::rename('Treatment' = 'X') %>% 
  dplyr::rename('Sample' = 'Name.BL')


## fix sample names
data.linger.sample = data.linger.sample %>% 
  mutate(Sample = paste0('Sample_',
                         Sample,
                         '_BL')) %>% 
  column_to_rownames('Sample')

### load sample meta data
data.linger.sample.meta = read.csv('data/treeSwallowCounts.samples.meta.linger.clean.csv')

## fix sample names
data.linger.sample.meta = data.linger.sample.meta %>% 
  mutate(Sample = str_replace_all(Sample,
                                  '-',
                                  '_')) %>% 
  mutate(Sample = paste0('Sample_',
                         Sample,
                         '_BL'))

# check matching names
data.linger.sample.meta %>% 
  filter(Sample %in% rownames(data.linger.sample)) %>% 
  nrow()

### get annotation data
## get chromosome data
data.anno = read_tsv('data/barn_swallow_ncbi_gene.tsv')

# create gene categories 
data.anno = data.anno %>% 
  mutate(Gene_category = case_when(str_detect(Name,
                                              'globin') ~ 'globin',
                                   str_detect(Name,
                                              'heat shock') ~ 'HSP',
                                   Chromosome == 'Z' ~ 'Z',
                                   Chromosome == 'W' ~ 'W',
                                   is.na(Chromosome) ~ 'Unplaced',
                                   TRUE ~ 'Auto'))

# check number of genes
data.anno %>% 
  count(Gene_category)

### add sex to sample data
## get list of samples by sex
# use Z:A ratio (chromosome 5) with threshold of 0.8
data.linger.sample.sex = data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto ,
         Sex = ifelse(Z_A_ratio > 0.8,
                      'Male',
                      'Female')) %>% 
  dplyr::select(Sample,
                Sex)
  
## add to sample data
data.linger.sample = data.linger.sample %>% 
  rownames_to_column('Sample') %>% 
  full_join(data.linger.sample.sex) %>% 
  column_to_rownames('Sample')

### prepare data
# reorder column names
data.linger.gene.mat = data.linger.gene.mat[, rownames(data.linger.sample)]

# check that columns and rows match
all(rownames(data.linger.sample) %in% colnames(data.linger.gene.mat))
all(rownames(data.linger.sample) == colnames(data.linger.gene.mat))

## filter out globin genes
data.linger.gene.mat = data.linger.gene.mat %>% 
  rownames_to_column('Gene') %>% 
  filter(!(Gene %in% c(data.anno %>% 
                         filter(Gene_category == 'globin') %>% 
                         pull(Symbol)))) %>% 
  column_to_rownames('Gene')


#### During: load data ####
### gene expression data
data.during.gene.wide = read.csv('data/data.during.gene.wide.csv')


# pivot longer
data.during.gene = data.during.gene.wide %>% 
  pivot_longer(cols = -c(Gene),
               names_to = 'Sample',
               values_to = 'Counts')

# create wide format for deseq2 
data.during.gene.mat = data.during.gene.wide %>% 
  column_to_rownames('Gene')

### load sample data
## create dataframe
data.during.sample = data.frame(Sample = colnames(data.during.gene.wide)[-1],
                                Treatment = c('Hot',
                                              'Hot',
                                              'Hot',
                                              'Con',
                                              'Con',
                                              'Con')) %>% 
  column_to_rownames('Sample')

### get annotation data
## get chromosome data
data.anno = read_tsv('data/barn_swallow_ncbi_gene.tsv')

# create gene categories 
data.anno = data.anno %>% 
  mutate(Gene_category = case_when(str_detect(Name,
                                              'globin') ~ 'globin',
                                   Chromosome == 'Z' ~ 'Z',
                                   Chromosome == 'W' ~ 'W',
                                   is.na(Chromosome) ~ 'Unplaced',
                                   TRUE ~ 'Auto'))

# # check number of genes
# data.anno %>% 
#   count(Gene_category)

### add sex to sample data
## get list of samples by sex
# use Z:A ratio (chromosome 5) with threshold of 0.8
data.during.sample.sex = data.during.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto ,
         Sex = ifelse(Z_A_ratio > 0.8,
                      'Male',
                      'Female')) %>% 
  dplyr::select(Sample,
                Sex)

## add to sample data
data.during.sample = data.during.sample %>% 
  rownames_to_column('Sample') %>% 
  full_join(data.during.sample.sex) %>% 
  column_to_rownames('Sample')

### prepare data
# reorder column names
data.during.gene.mat = data.during.gene.mat[, rownames(data.during.sample)]

# check that columns and rows match
all(rownames(data.during.sample) %in% colnames(data.during.gene.mat))
all(rownames(data.during.sample) == colnames(data.during.gene.mat))

## filter out globin genes
data.during.gene.mat = data.during.gene.mat %>% 
  rownames_to_column('Gene') %>% 
  filter(!(Gene %in% c(data.anno %>% 
                         filter(Gene_category == 'globin') %>% 
                         pull(Symbol)))) %>% 
  column_to_rownames('Gene')


#### Sample meta data ####
### Heat and mass stats
## descriptive stats
# by treatment
psych::describeBy(data.linger.sample.meta %>% 
                        dplyr::select(
                                      Avg.Trial.Temp,
                                      Avg.Trial.Temp.1,
                                      Mass_D1,
                                      Mass_D2) %>% 
                    mutate(Mass_diff = Mass_D1-Mass_D2),
                      na.rm = T,
                    group = data.linger.sample.meta$Treatment)

# by sex
psych::describeBy(data.linger.sample.meta %>% 
                    left_join(data.linger.sample.sex) %>% 
                    dplyr::select(
                      Avg.Trial.Temp,
                      Avg.Trial.Temp.1,
                      Mass_D1,
                      Mass_D2) %>% 
                    mutate(Mass_diff = Mass_D1-Mass_D2),
                  na.rm = T,
                  group = data.linger.sample.meta %>% 
                    left_join(data.linger.sample.sex) %>% 
                    pull(Sex))

# by sex and treatment
psych::describeBy(data.linger.sample.meta %>% 
                    left_join(data.linger.sample.sex) %>% 
                    dplyr::select(
                      Avg.Trial.Temp,
                      Avg.Trial.Temp.1,
                      Mass_D1,
                      Mass_D2) %>% 
                    mutate(Mass_diff = Mass_D1-Mass_D2),
                  na.rm = T,
                  group = data.linger.sample.meta %>% 
                    left_join(data.linger.sample.sex) %>% 
                    mutate(ID=paste0(Sex,Treatment)) %>% 
                    pull(ID))

# overall
psych::describe(data.linger.sample.meta %>% 
                    left_join(data.linger.sample.sex) %>% 
                    dplyr::select(
                      Avg.Trial.Temp,
                      Avg.Trial.Temp.1,
                      Mass_D1,
                      Mass_D2) %>% 
                    mutate(Mass_diff = Mass_D1-Mass_D2),
                  na.rm = T)


## stats
# heat trial day
t.test(Avg.Trial.Temp ~ Treatment,
       data = data.linger.sample.meta)

# heat next day
t.test(Avg.Trial.Temp.1 ~ Treatment,
       data = data.linger.sample.meta)

# mass
t.test(Mass_D1 ~ Treatment,
       data = data.linger.sample.meta)

# mass, sex
aov(Mass_D1 ~ Sex*Treatment,
    data = data.linger.sample.meta %>% 
      left_join(data.linger.sample.sex)) %>% 
  summary()

# mass diff
t.test(Mass_diff ~ Treatment,
       data = data.linger.sample.meta %>% 
         mutate(Mass_diff = Mass_D1-Mass_D2))

# mass diff, sex
aov(Mass_diff ~ Sex*Treatment,
       data = data.linger.sample.meta %>% 
         left_join(data.linger.sample.sex) %>% 
         mutate(Mass_diff = Mass_D1-Mass_D2)) %>% 
  summary()


### graph nest temperature
# scatter plot
data.linger.sample.meta %>% 
  left_join(data.linger.sample.sex) %>% 
  ggplot(aes(x = Avg.Trial.Temp,
             y = Avg.Trial.Temp.1,
             color = Treatment)) +
  geom_point(size = 3) +
  theme_classic() +
  ylab('Next day temp') +
  xlab('Trial temp') +
  scale_color_manual(values = c('blue',
                                'red')) +
  xlim(min(data.linger.sample.meta$Avg.Trial.Temp.1)-1,
       max(data.linger.sample.meta$Avg.Trial.Temp)+1)+
  ylim(min(data.linger.sample.meta$Avg.Trial.Temp.1)-1,
       max(data.linger.sample.meta$Avg.Trial.Temp)+1) +
  coord_fixed()
ggsave('figures/DEG_trial_temp/Nest Temperatures.png')

# boxplot
data.linger.sample.meta %>% 
  left_join(data.linger.sample.sex) %>% 
  pivot_longer(cols = c('Avg.Trial.Temp',
                        'Avg.Trial.Temp.1'),
               names_to = 'Temp_day',
               values_to = 'Nest_temp') %>% 
  mutate(Temp_day = ifelse(Temp_day == 'Avg.Trial.Temp',
                           'D12',
                           'D13')) %>% 
  ggplot(aes(x = Treatment,
             y = Nest_temp,
             color = Treatment)) +
  geom_boxplot(aes(group = Treatment),
               outlier.shape = NA) +
  geom_point(size = 3) +
  theme_classic() +
  ylab('Nest temp') +
  xlab('') +
  scale_color_manual(values = c('blue',
                                'red')) +
  facet_grid(.~Temp_day)
ggsave('figures/DEG_trial_temp/Nest Temperatures boxplot.png')
#### graph to check sex ####
## QC
# check proportion of total genes
# paper
data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  group_by(Gene_category,
           Sample) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(
    Percent = Sum.counts/sum(Sum.counts)) %>% 
  mutate(Sex = ifelse(Gene_category %in% c('Z','W'),
                      'Sex',
                      'Auto'),
         Sex.sample = case_when(Gene_category == 'Z' ~ Percent,
                                TRUE ~ 0),
         Sex.sample = max(Sex.sample)) %>%
  filter(Gene_category %in% c('globin',
                              'Auto')) %>% 
  ggplot(aes(x = Gene_category,
             y = Percent,
             fill = Gene_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  xlab('') +
  ylab('Percent of total reads') +
  scale_fill_manual(values = c('white',
                               'grey')) +
  ylim(0,1) +
  theme(legend.position = 'none')
ggsave('figures/paper/QC_linger_a.pdf',
       width = 3,
       height = 3)

## use total gene counts
# compare Z and W counts with chromosome
data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5,
                           23)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Chromosome) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'Chromosome',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/`5`,
         W_A_ratio = W/`23`)  %>% 
  left_join(data.linger.sample.sex) %>% 
  dplyr::select(Sex,
                Z_A_ratio,
                W_A_ratio) %>% 
  psych::describe.by(group = 'Sex',
                     na.rm = T)

### graph sex
# check proportion of total genes 
data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  group_by(Gene_category,
        Sample) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(
         Percent = Sum.counts/sum(Sum.counts)) %>% 
  mutate(Sex = ifelse(Gene_category %in% c('Z','W'),
                      'Sex',
                      'Auto'),
         Sex.sample = case_when(Gene_category == 'Z' ~ Percent,
                                TRUE ~ 0),
         Sex.sample = max(Sex.sample)) %>%
  ggplot(aes(x = reorder(Sample,
                         Sex.sample),
             y = Percent,
             fill = Gene_category)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1)) +
  facet_grid(Sex~.,
             scales = 'free') +
  xlab('') +
  ylab('Percent of total reads')
ggsave('figures/QC/Gene category percent of reads.png',
       width = 20,
       height = 10)

# compare Z counts with chromosome
data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto,
         W_A_ratio = W/Auto) %>% 
  ggplot(aes(x = Z_A_ratio,
             y = W_A_ratio)) +
  geom_vline(xintercept = 1,
             linetype = 'dashed') +
  geom_vline(xintercept = 0.5,
             linetype = 'dashed') +
  geom_point() +
  theme_classic() +
  xlab('Z to Chr 5 ratio') +
  ylab('W to Chr 5 ratio') 
ggsave('figures/QC/Sex Z_A ratio vs W reads.png')

# compare Z counts with chromosome
# paper
data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto,
         W_A_ratio = W/Auto) %>% 
  ggplot(aes(x = Z_A_ratio,
             y = W_A_ratio)) +
  geom_vline(xintercept = 1,
             linetype = 'dashed') +
  geom_vline(xintercept = 0.5,
             linetype = 'dashed') +
  geom_point(size = 3) +
  theme_classic() +
  xlab('Z to Chr 5 reads ratio') +
  ylab('W to Chr 5 reads ratio') 
ggsave('figures/paper/QC_linger_b.pdf',
       height = 3,
       width = 3)

# add label
data.linger.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto,
         W_A_ratio = W/Auto) %>% 
  ggplot(aes(x = Z_A_ratio,
             y = W_A_ratio)) +
  geom_vline(xintercept = 1,
             linetype = 'dashed') +
  geom_vline(xintercept = 0.5,
             linetype = 'dashed') +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = Sample),
                            max.overlaps = 20) +
  theme_classic() +
  xlab('Z to Chr 5 ratio') +
  ylab('W to Chr 5 ratio') 
ggsave('figures/QC/Sex Z_A ratio vs W reads label.png',
       height = 10,
       width = 10)

# #### get GO terms ####
# library(mygene)
# 
# 
# # 
# # Your starting data frame matching your counts
# gene_counts <- data.frame(
#   gene_id = c("tree_swallow000729", "tree_swallow000XXX"),
#   refseq_id = c("XP_018772500.1", "XP_014724698.1")
# )
# 
# # 1. Strip the version decimals (important!)
# gene_counts$refseq_clean <- gsub("\\..*$", "", gene_counts$refseq_id)
# 
# # 2. Query MyGene for functional GO annotations across multiple species
# go_results <- getGenes(gene_counts$refseq_clean, fields = c("go", "symbol", "taxid"))
# 
# # 3. Clean and isolate GO Molecular Function / Biological Process IDs
# processed_go <- as.data.frame(go_results) %>%
#   rowwise() %>%
#   mutate(
#     GO_IDs = paste(unique(c(go.BP[[1]]$id, go.MF[[1]]$id, go.CC[[1]]$id)), collapse = "; ")
#   ) %>%
#   select(query, symbol, taxid, GO_IDs)
# 
# # Merge back with your original table
# final_go_table <- merge(gene_counts, processed_go, by.x = "refseq_clean", by.y = "query")
# 
# 
# 
# 
# # BiocManager::install('mygene')
# library(mygene)
# 
# # 1. Your starting tree swallow count-table keys and mixed IDs
# data.anno.chr <- data.anno %>% 
#   dplyr::select(TranscriptID,
#                 AvesBlastID) 
# 
# # Strip version decimals so the database recognizes them
# data.anno.chr$refseq_clean <- gsub("\\..*$", "", data.anno.chr$AvesBlastID)
# 
# 
# # 1. Verify your data vector has no missing values
# refseq_vector <- na.omit(data.anno.chr$refseq_clean)
# 
# 
# # 2. Use queryManyDb (NOT getGenes) and specify scopes="refseq"
# # We query universally first to map the mixed protein IDs to their core gene info
# lookup_results <- queryMany(
#   refseq_vector,
#   scopes = "refseq",                   # Tells the API these are XP_ / NP_ / NM_ identifiers
#   fields = c("symbol",
#              "go") # Asks for the symbol and native genomic data
#   # species = 43150                      # Force the return architecture of the Barn Swallow
# )
# 
# 
# 
# lookup_df2 <- lookup_df %>%
#   group_by(query) %>%
#   mutate(
#     go.BP = paste(unique(c(go.BP[[1]]$id)), collapse = "; "),
#     go.MF = paste(unique(c(go.MF[[1]]$id)), collapse = "; "),
#     go.CC = paste(unique(c(go.CC[[1]]$id)), collapse = "; ")
#   ) %>%
#   dplyr::select(query, symbol, go.BP, go.MF, go.CC)
# 
# 
# # Merge the symbols back to your original dataframe
# data.anno.chr$gene_symbol <- lookup_df2$symbol
# 
# # Remove any genes that couldn't be converted to a symbol to avoid API errors
# clean_symbols <- lookup_df2$symbol[!is.na(lookup_df2$symbol)]
# 
# 
# # --- STEP 2: Species-Specific Query via Gene Symbols ---
# # Now use standard symbols to query the bird genome of your choice.
# # We use queryManyDb() instead of getGenes() because it handles symbol text searches safely.
# 
# swallow_orthologs <- getGenes(
#   clean_symbols,
#   scopes = "symbol",
#   fields = c("symbol", "genomic_pos", "genomic_pos_hg19"),
#   species = 43150 # Your target bird TaxID (e.g., 72873 for Hirundo rustica)
# )
# 
# 
# 
# # Convert to a standard data frame safely without hitting the rename bug
# swallow_df <- as.data.frame(swallow_orthologs)
# 
# head(swallow_df)
# 
# # --- STEP 3: Extract Chromosome Z safely ---
# # Parse the nested results to isolate your Z-linked genes
# final_z_map <- swallow_df %>%
#   rowwise() %>%
#   mutate(
#     # Pull the chromosome name out of the nested lists safely
#     Chromosome = ifelse(!is.null(genomic_pos$chr), genomic_pos$chr, 
#                         ifelse(!is.null(genomic_pos_hg19$chr), genomic_pos_hg19$chr, NA))
#   ) %>%
#   dplyr::select(query, symbol, Chromosome) %>%
#   filter(tolower(Chromosome) == "z") # Captures both "Z" and "z"
# 
# 
# 
# 
# 
# 
# 
# 
# # 2. Query MyGene to get the Barn Swallow coordinates via orthology
# # We pass "species=147047" to force the structural metadata to come from Hirundo rustica
# swallow_orthologs <- getGenes(
#   data.anno.chr$refseq_clean, 
#   fields = c("symbol"))
# 
# as.data.frame(swallow_orthologs) %>%
#   head()
# 
# # 3. Process the results into a clean dataframe
# # MyGene returns the chromosome in the 'genomic_pos.chr' or 'genomic_pos_hg19.chr' nested field
# processed_coords <- as.data.frame(swallow_orthologs) %>%
#   rowwise() %>%
#   mutate(
#     # Handle variations in how the API returns chromosome fields for non-model birds
#     Barn_Swallow_Chr = ifelse(!is.null(genomic_pos$chr), genomic_pos$chr, 
#                               ifelse(!is.null(genomic_pos_hg19$chr), genomic_pos_hg19$chr, NA)),
#     Start = ifelse(!is.null(genomic_pos$start), genomic_pos$start, NA)
#   ) %>%
#   dplyr::select(query, symbol, Barn_Swallow_Chr, Start)
# 
# # 4. Merge back to your original Tree Swallow Gene Table
# final_table <- merge(my_data, processed_coords, by.x = "refseq_clean", by.y = "query", all.x = TRUE)
# 
# # 5. Filter specifically for Z-linked genes
# z_linked_genes <- final_table %>% 
#   filter(Barn_Swallow_Chr == "Z" | Barn_Swallow_Chr == "z")

#### Use DESQ2 ####
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### create DESEQ dataset
dds <- DESeqDataSetFromMatrix(countData = data.linger.gene.mat,
                              colData = data.linger.sample,
                              design = ~ Treatment + Sex)
# check data
dds

# remove rows with low counts
# need to have a count of at least 5 and be present in at least half of samples
keep = rowSums(counts(dds) >= 5) >= nrow(data.linger.sample)/2

# remove low expressed genes
dds = dds[keep,]

#check data
dds
# 9304 genes

### run DEG
## run DESeq
dds = DESeq(dds)

## graph mean variance
# save graph
png('figures/QC/Deseq2 mean variance.png')
plotDispEsts(dds)
dev.off()


## get results
# heat
res.treat = results(dds,
              contrast = c('Treatment',
                           'Hot',
                           'Con'))

summary(res.treat)

# check direction
png('figures/DEG/Heat direction check.png')
plotCounts(dds,
           gene=which.min(res.treat$padj),
           intgroup="Treatment")
dev.off()


# create dataframe
res.treat.df = res.treat %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.treat.df,
          'data/res.treat.df.csv',
          row.names = F)

# sex
res.sex = results(dds,
                    contrast = c('Sex',
                                 'Male',
                                 'Female'))

summary(res.sex)

# check direction
png('figures/DEG/Sex direction check.png')
plotCounts(dds,
           gene=which.min(res.sex$padj),
           intgroup="Sex")
dev.off()

# create dataframe
res.sex.df = res.sex %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.sex.df,
          'data/res.sex.df.csv',
          row.names = F)


#### graph DEG results ####
### create volcano plot
## heat 
res.treat.df %>% 
  mutate(Sig = ifelse(padj <= 0.153,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Sig == 'Sig',
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.15),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Heat DEG (padj < 0.153)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG/Heat volcano plot.png',
       height = 10,
       width = 10)  

# get chromosome counts
res.treat.df %>% 
  mutate(Sig = case_when(padj <= 0.155 & log2FoldChange > 0 ~ 'Hot',
                         padj <= 0.155 & log2FoldChange < 0 ~ 'Con',
                         TRUE ~ 'none')) %>%
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Sig,
        Gene_category)

## heat boxplots
## get DEGs
res.treat.df.deg = res.treat.df %>% 
  filter(padj < 0.16) %>% 
  pull(Gene)

# normalize and scale counts
deg.mat = assay(vst(dds,
                    blind = F))[res.treat.df.deg,] %>% 
  t() %>% 
  scale() %>%
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(-c(Sample),
               names_to = 'Gene',
               values_to = 'z.score') %>% 
  left_join(data.linger.sample.meta) %>% 
  left_join(data.linger.sample.sex)


## graph each DEG boxplots
# treatment temp
# paper
deg.mat %>% 
  mutate(ID = paste0(Treatment,
                     '.',
                     Sex)) %>% 
  ggplot(aes(x = ID,
             y = z.score,
             color = Treatment)) +
  geom_boxplot(aes(group = ID),
               outlier.shape = NA) +
  geom_point(size = 2.5,
             aes(shape = Sex)) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab('Gene expression (z-score)') +
  xlab('') +
  scale_color_manual(values = c('#464b9f',
                                '#f26622')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 3) 
ggsave('figures/paper/Figure 1.pdf',
       height = 3,
       width = 6.5,
       dpi = 720)

## Sex 
res.sex.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Gene %in% c(res.sex.df %>% 
                                      mutate(direction = ifelse(log2FoldChange > 0,
                                                                'up',
                                                                'down')) %>% 
                                      group_by(direction) %>% 
                                      slice_min(order_by = padj,
                                                n = 5) %>% 
                                      pull(Gene)),
                        Gene,
                        NA),
         Label = ifelse(abs(log2FoldChange) > 5.9,
                        Gene,
                        Label)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Sex DEG (padj < 0.05)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG/Sex volcano plot.png',
       height = 10,
       width = 10)  

# get sex counts
res.sex.df %>% 
  mutate(Sig = case_when(padj <= 0.01 & log2FoldChange > 0 ~ 'Male',
                         padj <= 0.01 & log2FoldChange < 0 ~ 'Female',
                         TRUE ~ 'none')) %>%
  count(Sig)

# check chromosome
res.sex.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig')) %>%
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Gene_category,
                 shape = Sig),
             size = 5) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Sex DEG (padj < 0.05)') 
ggsave('figures/DEG/Sex volcano plot chrom.png',
       height = 10,
       width = 10) 

# paper
res.sex.df %>% 
  mutate(Sig = ifelse(padj <= 0.01,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Gene %in% c(res.sex.df %>% 
                                      mutate(direction = ifelse(log2FoldChange > 0,
                                                                'up',
                                                                'down')) %>% 
                                      group_by(direction) %>% 
                                      slice_min(order_by = padj,
                                                n = 5) %>% 
                                      pull(Gene)),
                        Gene,
                        NA),
         Label = ifelse(abs(log2FoldChange) > 5.9,
                        Gene,
                        Label)) %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  mutate(Chromosome_type = case_when(Chromosome == 'Z' ~ 'Z',
                                     Chromosome == 'W' ~ 'W',
                                     TRUE ~ 'Auto'),
         Sig = ifelse(is.na(Sig),
                      'Not sig',
                      Sig)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.01),
             linetype = 'dashed') +
  geom_point(aes(color = Sig,
                 shape = Chromosome_type),
             size = 3) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Sex DEG (padj < 0.01)') +
  scale_color_manual(values = c('grey',
                                'black')) +
  scale_shape_manual(values = c(16,
                                17,
                                15)) +
  labs(shape = 'Chromosome',
       color = 'Significance') +
  theme(legend.position = 'inside' ,
        legend.position.inside = c(0.9,
                                   0.75))
ggsave('figures/paper/DEG_linger_sex_1.pdf',
       height = 6.5,
       width = 6.5)  

# get chromosome counts
res.sex.df %>% 
  mutate(Sig = case_when(padj <= 0.01 & log2FoldChange > 0 ~ 'Male',
                         padj <= 0.01 & log2FoldChange < 0 ~ 'Female',
                         TRUE ~ 'none')) %>%
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Sig,
        Gene_category)
  
  
#### normalize data with vst and graph PCA ####
#VST
vsd = vst(dds,
          blind = T)


## convert to matrix
vsd.mat = assay(vsd)

# get top 1000 most variable genes
rv = rowVars(vsd.mat)
select = order(rv,
               decreasing = TRUE)[seq_len(min(1000, length(rv)))]

# Run PCA on transposed matrix
pca <- prcomp(t(vsd.mat[select, ]))

# Access PC values for each sample
pca_results <- as.data.frame(pca$x)

# the contribution to the total variance for each component
percentVar.tmp <- pca$sdev^2 / sum(pca$sdev^2 )

### compare loadings to chromosome
## check loadings
# top and bottom 10
pca$rotation %>% 
  as.data.frame() %>% 
  slice_max(order_by = PC1,
            n = 10) %>% 
  bind_rows(pca$rotation %>% 
              as.data.frame() %>% 
          slice_min(order_by = PC1,
                    n = 10)) %>% 
  rownames_to_column('Symbol') %>% 
  left_join(data.anno) %>% 
  mutate(sign = sign(PC1)) %>% 
    count(sign,
          Gene_category)

# top 10 overall 
pca$rotation %>% 
  as.data.frame() %>% 
  mutate(abs.pc1 = abs(PC1)) %>% 
  slice_max(order_by = abs.pc1,
            n = 10) %>% 
  rownames_to_column('Symbol') %>% 
  left_join(data.anno) %>% 
  mutate(sign = sign(PC1)) %>% 
  count(sign,
        Gene_category)
  

#### graph PCA
### graph PCA var
data.frame(percentVar = percentVar.tmp,
           PC = seq(1:length(percentVar.tmp))) %>%
  mutate(PC = paste('PC',
                    PC,
                    sep ='')) %>%
  filter(percentVar >= 0.01) %>%
  ggplot(aes(x = reorder(PC,
                         -percentVar),
             y = percentVar)) +
  geom_point() +
  geom_segment(aes(x=reorder(PC,
                             -percentVar),
                   xend=reorder(PC,
                                -percentVar),
                   y=0,
                   yend=percentVar)) +
  theme_classic() +
  ggtitle('PC variance') +
  xlab('PCs') +
  ylab('Percent variance')
ggsave('figures/QC/PCA/PC variance.png')

### plot pca
## PC 1 vs PC 2
# treatment
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar.tmp[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar.tmp[2]),"% variance")) +
  coord_fixed() +
  theme_classic() 
ggsave('figures/QC/PCA/PCA PC1 vs PC2.png')

# paper
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar.tmp[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar.tmp[2]),"% variance")) +
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('#464b9f',
                                '#f26622')) +
  theme(legend.position = 'none')
ggsave('figures/paper/QC_linger_2b_nolegend.pdf',
       height = 3,
       width = 3.25)

## PC 3 vs PC 4
# treatment
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC3,
             y = PC4,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC3: ",round(100*percentVar.tmp[3]),"% variance")) +
  ylab(paste0("PC4: ",round(100*percentVar.tmp[4]),"% variance")) +
  coord_fixed() +
  theme_classic() 
ggsave('figures/QC/PCA/PCA PC3 vs PC4.png')


### calculate percent variation due to sex and treatment
## use variance partition on normalized 
varPart = fitExtractVarPartModel(assay(
  vst(
    dds,
    blind = F)), 
  formula = ~ Treatment + Sex, 
  data.linger.sample)

# graph global distribution
plotVarPart(sortCols(varPart))
ggsave('figures/DEG/Variance partition.png')

# paper
plotVarPart(sortCols(varPart)) +
  ggtitle('~ Treatment + Sex')
ggsave('figures/paper/DEG_linger_QC_1a.pdf',
       height = 3,
       width = 3.25)

# average percentage of variation due to treatment across the genome
mean(varPart$Treatment) * 100
# 4.694733

# average percentage of variation due to treatment across the genome
mean(varPart$Sex) * 100
# 11.67194

#### graph VST normalization ####
## plot row sd vs row mean
# normalized
png('figures/QC/Vst normalization.png')
meanSdPlot(assay(vsd))
dev.off()


## graph expression per sample
# raw
data.frame(assay(dds)) %>%
  mutate(
    Gene_id = row.names(dds)
  ) %>%
  pivot_longer(-Gene_id) %>%
  ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Expression Linger",
    x = "sample",
    y = "expression"
  )
ggsave('figures/QC/Expression per sample linger.png')

# normalized
data.frame(assay(vsd)) %>%
  mutate(
    Gene_id = row.names(vsd)
  ) %>%
  pivot_longer(-Gene_id) %>%
  ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized Expression Linger",
    x = "sample",
    y = "normalized expression"
  )
ggsave('figures/QC/VST Normalized Expression per sample linger.png')

### check data for outliers
### create pvclust matrix for each tissue
## raw expression
# create tmp dataframe
tmp.data = assay(dds)
# only use selected genes
tmp.data = tmp.data[select,]

## run pvclust
tmp.p = pvclust::pvclust(tmp.data,
                         method.dist="cor",
                         method.hclust="average",
                         nboot=1000,
                         quiet = F)

# graph
png("figures/QC/Sample dendrogam counts.png",
    height = 10,
    width = 10,
    unit = 'in',
    res = 480)
plot(tmp.p)
dev.off()

## normalized expression
# create tmp dataframe
tmp.data = assay(vsd)
# only use selected genes
tmp.data = tmp.data[select,]

## run pvclust
tmp.p = pvclust::pvclust(tmp.data,
                         method.dist="cor",
                         method.hclust="average",
                         nboot=1000,
                         quiet = F)
# graph
png("figures/QC/Sample dendrogam counts normalized.png",
    height = 10,
    width = 10,
    unit = 'in',
    res = 480)
plot(tmp.p)
dev.off()

# paper
pdf("figures/paper/QC_linger_2a.pdf",
    height = 6,
    width = 6.5)
plot(tmp.p)
dev.off()


#### Remove sex chromsomes: normalize data with vst and graph PCA ####
### create DESEQ dataset
# remove sex chromsomes from dds
dds.no.sex <- DESeqDataSetFromMatrix(countData = data.linger.gene.mat[c(data.anno %>% 
                                                                          filter(Symbol %in% rownames(dds)) %>% 
                                                                          filter(!(Chromosome %in% c('Z',
                                                                                                     'W'))) %>% 
                                                                          pull(Symbol) %>% 
                                                                          unique()),],
                              colData = data.linger.sample,
                              design = ~ Treatment + Sex)
# check data
dds.no.sex

# remove rows with low counts
# need to have a count of at least 5 and be present in at least half of samples
keep = rowSums(counts(dds.no.sex) >= 5) >= nrow(data.linger.sample)/2

# remove low expressed genes
dds.no.sex = dds.no.sex[keep,]

#check data
dds.no.sex
# 8762 genes

#VST
vsd.no.sex = vst(dds.no.sex,
          blind = T)


## convert to matrix
vsd.no.sex.mat = assay(vsd.no.sex)

# get top 1000 most variable genes
rv = rowVars(vsd.no.sex.mat)
select = order(rv,
               decreasing = TRUE)[seq_len(min(1000, length(rv)))]

# Run PCA on transposed matrix
pca <- prcomp(t(vsd.no.sex.mat[select, ]))

# Access PC values for each sample
pca_results <- as.data.frame(pca$x)

# the contribution to the total variance for each component
percentVar.tmp <- pca$sdev^2 / sum(pca$sdev^2 )

### compare loadings to chromosome
## check loadings
# top and bottom 10
pca$rotation %>% 
  as.data.frame() %>% 
  slice_max(order_by = PC1,
            n = 10) %>% 
  bind_rows(pca$rotation %>% 
              as.data.frame() %>% 
              slice_min(order_by = PC1,
                        n = 10)) %>% 
  rownames_to_column('Symbol') %>% 
  left_join(data.anno) %>% 
  mutate(sign = sign(PC1)) %>% 
  count(sign,
        Gene_category)

# top 10 overall 
pca$rotation %>% 
  as.data.frame() %>% 
  mutate(abs.pc1 = abs(PC1)) %>% 
  slice_max(order_by = abs.pc1,
            n = 10) %>% 
  rownames_to_column('Symbol') %>% 
  left_join(data.anno) %>% 
  mutate(sign = sign(PC1)) %>% 
  count(sign,
        Gene_category)


#### graph PCA
### graph PCA var
data.frame(percentVar = percentVar.tmp,
           PC = seq(1:length(percentVar.tmp))) %>%
  mutate(PC = paste('PC',
                    PC,
                    sep ='')) %>%
  filter(percentVar >= 0.01) %>%
  ggplot(aes(x = reorder(PC,
                         -percentVar),
             y = percentVar)) +
  geom_point() +
  geom_segment(aes(x=reorder(PC,
                             -percentVar),
                   xend=reorder(PC,
                                -percentVar),
                   y=0,
                   yend=percentVar)) +
  theme_classic() +
  ggtitle('PC variance') +
  xlab('PCs') +
  ylab('Percent variance')
# ggsave('figures/QC/PCA/PC variance.png')

### plot pca
## PC 1 vs PC 2
# treatment
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar.tmp[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar.tmp[2]),"% variance")) +
  coord_fixed() +
  theme_classic() 
# ggsave('figures/QC/PCA/PCA PC1 vs PC2.png')

## PC 3 vs PC 4
# treatment
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC3,
             y = PC4,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC3: ",round(100*percentVar.tmp[3]),"% variance")) +
  ylab(paste0("PC4: ",round(100*percentVar.tmp[4]),"% variance")) +
  coord_fixed() +
  theme_classic() 
# ggsave('figures/QC/PCA/PCA PC3 vs PC4.png')


### calculate percent variation due to sex and treatment
## use variance partition on normalized 
varPart.no.sex = fitExtractVarPartModel(assay(
  vst(
    dds.no.sex,
    blind = F)), 
  formula = ~ Treatment + Sex, 
  data.linger.sample)

# graph global distribution
plotVarPart(sortCols(varPart.no.sex)) +
  ggtitle('No sex chr')
ggsave('figures/DEG/Variance partition no sex chr.png')

# paper
plotVarPart(sortCols(varPart.no.sex)) +
  ggtitle('~ Treatment + Sex (no Z/W)')
ggsave('figures/paper/DEG_linger_QC_1b.pdf',
       height = 3,
       width = 3.25)

# average percentage of variation due to treatment across the genome
mean(varPart.no.sex$Treatment) * 100
# 4.806704

# average percentage of variation due to treatment across the genome
mean(varPart.no.sex$Sex) * 100
# 9.299909

#### PCA remove sex effect ####
### use limma to remove sex batch effects during normalization
vsd.sex.mat = limma::removeBatchEffect(assay(vsd),
                               batch = vsd$Sex)


# get top 1000 most variable genes
rv.sex = rowVars(vsd.sex.mat)
select.sex = order(rv.sex,
               decreasing = TRUE)[seq_len(min(1000, length(rv.sex)))]

# Run PCA on transposed matrix
pca.sex <- prcomp(t(vsd.sex.mat[select.sex, ]))

# Access PC values for each sample
pca.sex_results <- as.data.frame(pca.sex$x)

# the contribution to the total variance for each component
percentVar.tmp <- pca.sex$sdev^2 / sum(pca.sex$sdev^2 )

#### graph PCA
### graph PCA var
data.frame(percentVar = percentVar.tmp,
           PC = seq(1:length(percentVar.tmp))) %>%
  mutate(PC = paste('PC',
                    PC,
                    sep ='')) %>%
  filter(percentVar >= 0.01) %>%
  ggplot(aes(x = reorder(PC,
                         -percentVar),
             y = percentVar)) +
  geom_point() +
  geom_segment(aes(x=reorder(PC,
                             -percentVar),
                   xend=reorder(PC,
                                -percentVar),
                   y=0,
                   yend=percentVar)) +
  theme_classic() +
  ggtitle('PC variance') +
  xlab('PCs') +
  ylab('Percent variance')
ggsave('figures/QC/PCA_sex/PC variance.png')

### plot pca
## PC 1 vs PC 2
# treatment
pca.sex$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar.tmp[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar.tmp[2]),"% variance")) +
  coord_fixed() +
  theme_classic() 
ggsave('figures/QC/PCA_sex/PCA PC1 vs PC2.png')

## PC 3 vs PC 4
# treatment
pca.sex$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.linger.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC3,
             y = PC4,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC3: ",round(100*percentVar.tmp[3]),"% variance")) +
  ylab(paste0("PC4: ",round(100*percentVar.tmp[4]),"% variance")) +
  coord_fixed() +
  theme_classic() 
ggsave('figures/QC/PCA_sex/PCA PC3 vs PC4.png')





#### Use DESQ2 for interaction ####
### create DESEQ dataset
dds.interact <- DESeqDataSetFromMatrix(countData = data.linger.gene.mat,
                                       colData = data.linger.sample,
                                       design = ~ Treatment * Sex)
# check data
dds.interact

# remove rows with low counts
# need to have a count of at least 5 and be present in at least half of samples
keep = rowSums(counts(dds.interact) >= 5) >= nrow(data.linger.sample)/2

# remove low expressed genes
dds.interact = dds.interact[keep,]

#check data
dds.interact
# 9304 genes

### run DEG
## run DESeq
dds.interact = DESeq(dds.interact)

## get interaction results
# use LRT 
dds.interact.lrt = DESeq(dds.interact, 
                         test="LRT", 
                         reduced=~ Treatment + Sex)

# get results
# genes with significant interaction
res.interact.lrt = results(dds.interact.lrt)

# create dataframe
res.interact.lrt.df = res.interact.lrt %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.interact.lrt.df,
          'data/res.interact.lrt.df.csv',
          row.names = F)

### sex specific heat response
## run sex wald
dds.interact.wald = DESeq(dds.interact,
                          test = 'Wald')
  

# get names of interactions
resultsNames(dds.interact.wald)

# heat response in females
res.females.heat = results(dds.interact.wald, 
                         name="Treatment_Hot_vs_Con")
summary(res.females.heat)

# create dataframe
res.females.heat.df = res.females.heat %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.females.heat.df,
          'data/res.females.heat.df.csv',
          row.names = F)

# heat response in males
res.males.heat = results(dds.interact.wald, 
                        contrast=list("Treatment_Hot_vs_Con", 
                                      "TreatmentHot.SexMale"))
summary(res.males.heat)

# create dataframe
res.males.heat.df = res.males.heat %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.males.heat.df,
          'data/res.males.heat.df.csv',
          row.names = F)

## create list of interaction genes
# signifcant 0.1 adjusted pvalue
res.interact.df = res.interact.lrt.df %>% 
  dplyr::select(Gene,
                padj) %>% 
  inner_join(res.males.heat.df %>% 
               dplyr::select(Gene,
                             log2FoldChange)) %>% 
  inner_join(res.females.heat.df %>% 
               dplyr::select(Gene,
                             log2FoldChange),
             by = 'Gene',
             suffix = c('_male',
                        '_female')) %>% 
  mutate(delta = log2FoldChange_male - log2FoldChange_female,
         sex_modulated = case_when(padj < 0.1 & delta > 1 ~ 'male.stronger',
                                   padj < 0.1 & delta < -1 ~ 'female.stronger',
                                   TRUE ~ NA),
         direction = case_when(log2FoldChange_male > 0 & log2FoldChange_female > 0 ~ "up.up",
                               log2FoldChange_male < 0 & log2FoldChange_female < 0 ~ "down.down",
                               TRUE ~ "discordant"),
         
         sex_modulated_direction = case_when(
           sex_modulated == 'male.stronger'   & direction == "up.up"   ~ 'male.stronger.up',
           sex_modulated == 'male.stronger'   & direction == "down.down" ~ 'male.stronger.down',
           sex_modulated == 'female.stronger' & direction == "up.up"   ~ 'female.stronger.up',
           sex_modulated == 'female.stronger' & direction == "down.down" ~ 'female.stronger.down',
           sex_modulated == 'male.stronger'   & direction == "discordant" ~ 'male.discordant',
           sex_modulated == 'female.stronger' & direction == "discordant" ~ 'female.discordant',
           TRUE ~ NA))


#### Compare interaction #### 
### calculate percent variation due to sex and treatment
## use variance partition on normalized 
varPart.interact = fitExtractVarPartModel(assay(
  vst(
    dds.interact,
    blind = F)), 
  formula = ~ Treatment * Sex, 
  data.linger.sample)

# graph global distribution
plotVarPart(sortCols(varPart.interact))
ggsave('figures/DEG/Variance partition interaction.png')

# average percentage of variation due to treatment across the genome
mean(varPart.interact$Treatment) * 100
# 11.32502

# average percentage of variation due to treatment across the genome
mean(varPart.interact$Sex) * 100
# 15.99451

# average percentage of variation due to interaction sex and treatment across the genome
mean(varPart.interact$`Treatment:Sex`) * 100
# 14.67263


## get DEG
res.treat.df.deg = res.treat.df %>% 
  filter(padj < 0.16) %>% 
  pull(Gene)

# # add sex specific DEG
# res.treat.df.deg = c(res.treat.df.deg,
#                      'RYBP',
#                      'LOC120765528'
#                      )

# normalize and scale counts
deg.mat = assay(vst(dds,
                    blind = F))[res.treat.df.deg,] %>% 
  t() %>% 
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(-c(Sample),
               names_to = 'Gene',
               values_to = 'z.score') %>% 
  left_join(data.linger.sample.meta) %>% 
  left_join(data.linger.sample.sex)


## graph each DEG boxplots
# treatment temp
deg.mat %>% 
  mutate(ID = paste0(Sex,
                     '.',
                     Treatment)) %>% 
  ggplot(aes(x = ID,
             y = z.score,
             color = Treatment)) +
  geom_boxplot(aes(group = ID),
               outlier.shape = NA) +
  geom_point(size = 3,
             aes(shape = Sex)) +
  theme_classic() +
  ylab('Normalized gene expression (z-score)') +
  xlab('') +
  scale_color_manual(values = c('blue',
                                'red')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 3)
ggsave('figures/DEG/Boxplot DEG expression.png',
       height = 10,
       width = 15)

## graph each DEG across temp
# treatment temp
deg.mat %>% 
  ggplot(aes(x = Avg.Trial.Temp,
             y = z.score,
             color = Treatment)) +
  geom_smooth(method = 'lm',
              se = FALSE) +
  geom_point(size = 3,
             aes(shape = Sex)) +
  theme_classic() +
  ylab('Normalized gene expression (z-score)') +
  xlab('Trial day temp') +
  scale_color_manual(values = c('blue',
                                'red')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 4)
ggsave('figures/DEG/Temp/Trial_day_temp vs DEG expression.png',
       height = 10,
       width = 15)

# paper
deg.mat %>% 
  ggplot(aes(x = Avg.Trial.Temp,
             y = z.score,
             color = Treatment)) +
  geom_smooth(method = 'lm',
              se = FALSE) +
  geom_point(size = 1.5,
             aes(shape = Sex)) +
  theme_classic() +
  ylab('Normalized gene expression (z-score)') +
  xlab('Next day temp') +
  scale_color_manual(values = c('#464b9f',
                                '#f26622')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 3) +
  theme(legend.position = 'none',
        strip.text = element_text(size = 6))
ggsave('figures/paper/DEG_linger_temp_1.pdf',
       height = 3.25,
       width = 3.25)

# next day temp
deg.mat %>% 
    ggplot(aes(x = Avg.Trial.Temp.1,
               y = z.score,
               color = Treatment)) +
    geom_smooth(method = 'lm',
                se = FALSE) +
    geom_point(size = 3,
               aes(shape = Sex)) +
    theme_classic() +
    ylab('Normalized gene expression (z-score)') +
    xlab('Next day temp') +
    scale_color_manual(values = c('blue',
                                  'red')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 4)
ggsave('figures/DEG/Temp/Next_day_temp vs DEG expression.png',
       height = 10,
       width = 15)



### compare across males and females
## create volcano plot
# get overlapping DEG
res.sex.both.heat.df.deg = res.males.heat.df %>% 
  mutate(Sig.male = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig')) %>% 
  dplyr::select(Sig.male,
                Gene) %>% 
  left_join(res.females.heat.df %>% 
              mutate(Sig.female = ifelse(padj <= 0.05,
                                       'Sig',
                                       'Not sig')) %>% 
              dplyr::select(Sig.female,
                            Gene)) %>% 
  mutate(Sig.both = ifelse(Sig.male == 'Sig' & Sig.female == 'Sig',
                           'Both',
                           NA)) %>% 
  filter(Sig.both == 'Both') %>% 
  pull(Gene)


# paper
res.males.heat.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(padj <= 0.001 | abs(log2FoldChange) > 7.5,
                        Gene,
                        NA),
         Both = ifelse(Gene %in% res.sex.both.heat.df.deg,
                        'Both',
                        'none'),
         Sig = ifelse(is.na(Sig),
                      'Not sig',
                      Sig),
         Sig.color = case_when(Sig == 'Sig' & log2FoldChange > 0 ~ 'Heat',
                               Sig == 'Sig' & log2FoldChange < 0 ~ 'Con',
                               TRUE ~ 'none')) %>%
  arrange(Sig) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig.color,
                 shape = Both),
             size = 3) +
  ggrepel::geom_label_repel(aes(label = Label),
                            max.overlaps = 50) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Male heat DEG') +
  scale_color_manual(values = c('#464b9f',
                                '#f26622',
                                'grey')) +
  theme(legend.position = 'none')
ggsave('figures/paper/Deg_linger_heat_sex_1a.pdf',
       height = 3.5,
       width = 3.25)  


# heat males
res.males.heat.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(padj <= 0.01 | abs(log2FoldChange) > 7.5,
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label),
                            max.overlaps = 50) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Male heat DEG (padj < 0.05)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG/Male heat volcano plot.png',
       height = 10,
       width = 10)  

# get heat male counts
res.males.heat.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Sig = case_when(Sig == 'Sig' & log2FoldChange > 0 ~ 'Hot',
                         Sig == 'Sig' & log2FoldChange < 0 ~ 'Con',
                         TRUE ~ 'none')) %>% 
  count(Sig)

# heat females
res.females.heat.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(padj <= 0.015 | abs(log2FoldChange) > 4.5,
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label),
                            max.overlaps = 50) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Female heat DEG (padj < 0.05)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG/Female heat volcano plot.png',
       height = 10,
       width = 10)  

# paper
res.females.heat.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(padj <= 0.01 | abs(log2FoldChange) > 4.5,
                        Gene,
                        NA),
         Both = ifelse(Gene %in% res.sex.both.heat.df.deg,
                       'Both',
                       'none'),
         Sig = ifelse(is.na(Sig),
                      'Not sig',
                      Sig),
         Sig.color = case_when(Sig == 'Sig' & log2FoldChange > 0 ~ 'Heat',
                               Sig == 'Sig' & log2FoldChange < 0 ~ 'Con',
                               TRUE ~ 'none')) %>%
  arrange(Sig) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig.color,
                 shape = Both),
             size = 3) +
  ggrepel::geom_label_repel(aes(label = Label),
                            max.overlaps = 50) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Female heat DEG')  +
  scale_color_manual(values = c('#464b9f',
                                '#f26622',
                                'grey')) +
  theme(legend.position = 'none')
ggsave('figures/paper/Deg_linger_heat_sex_1b.pdf',
       height = 3.5,
       width = 3.25) 

# get heat female counts
res.females.heat.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Sig = case_when(Sig == 'Sig' & log2FoldChange > 0 ~ 'Hot',
                         Sig == 'Sig' & log2FoldChange < 0 ~ 'Con',
                         TRUE ~ 'none')) %>% 
  count(Sig)



## graph logfoldchange differences
res.interact.df %>% 
  filter(is.na(sex_modulated)) %>% 
  ggplot(aes(x = log2FoldChange_male,
             y = log2FoldChange_female,
             color = sex_modulated)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 'dashed') +
  geom_point() +
  geom_point(data = res.interact.df %>% 
               filter(!is.na(sex_modulated)),
             aes(x = log2FoldChange_male,
                 y = log2FoldChange_female,
                 color = sex_modulated)) +
  annotate('text',
           hjust= 1,
           x = max(res.interact.df$log2FoldChange_male),
           y = min(res.interact.df$log2FoldChange_female),
           size = 3,
           label = paste0('R-squared = ',
                          cor(res.interact.df$log2FoldChange_male,
                              res.interact.df$log2FoldChange_female) %>% 
                            round(3)))+
  theme_classic()
ggsave('figures/DEG/Interaction scatterplot.png')


### create heatmap of interaction genes
## get top10 genes by interaction 
res.interact.lrt.df.top = res.interact.lrt.df %>% 
  mutate(direction = ifelse(log2FoldChange > 0,
                            'up',
                            'down')) %>% 
  filter(padj < 0.1) %>% 
  group_by(direction) %>% 
  slice_min(order_by = pvalue,
            n = 10) 

# normalize and scale counts
interaction.mat = assay(vst(dds.interact.lrt,
                            blind = F))[res.interact.lrt.df.top$Gene,] %>% 
  t() %>% 
  scale() %>% 
  t()

# create heatmap annotations
interaction.mat.anno = HeatmapAnnotation(df = data.linger.sample, 
                        col = list(Treatment = c("Con" = "blue",
                                                 "Hot" = "red"),
                                   Sex = c("Male" = "lightblue",
                                           "Female" = "lightpink")))

# create heatmap
png('figures/DEG/Top interaction heatmap.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 720)
Heatmap(interaction.mat, 
        name = "Z-score", 
        top_annotation = interaction.mat.anno,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        column_split = paste0(data.linger.sample$Sex,
                              "_", 
                              data.linger.sample$Treatment),
        row_split = 2,
        cluster_columns = TRUE)
dev.off()

## All genes by interaction 
res.interact.lrt.df.top = res.interact.lrt.df %>% 
  mutate(direction = ifelse(log2FoldChange > 0,
                            'up',
                            'down')) %>% 
  filter(padj < 0.1) %>% 
  group_by(direction) %>% 
  slice_min(order_by = pvalue,
            n = 150) 

# normalize and scale counts
interaction.mat = assay(vst(dds.interact.lrt,
                            blind = F))[res.interact.lrt.df.top$Gene,] %>% 
  t() %>% 
  scale() %>% 
  t()

# create heatmap annotations
interaction.mat.anno = HeatmapAnnotation(df = data.linger.sample, 
                                         col = list(Treatment = c("Con" = "blue",
                                                                  "Hot" = "red"),
                                                    Sex = c("Male" = "lightblue",
                                                            "Female" = "lightpink")))

# create heatmap
png('figures/DEG/All interaction heatmap.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 720)
Heatmap(interaction.mat, 
        name = "Z-score", 
        top_annotation = interaction.mat.anno,
        show_row_names = F, 
        show_column_names = TRUE,
        column_split = paste0(data.linger.sample$Sex,
                              "_", 
                              data.linger.sample$Treatment),
        row_split = 2,
        cluster_columns = TRUE)
dev.off()

## check count of interaction genes
res.interact.lrt.df %>% 
  filter(padj < 0.1) %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Gene_category)
# Gene_category   n
# 1          Auto 139
# 2           HSP   2
# 3      Unplaced   2
# 4             Z  11


### RRHO2
## compare male vs female heat response
rrho_output = RRHO2_initialize(res.males.heat.df %>% 
                                 mutate(metric = sign(log2FoldChange) * -log10(pvalue)) %>% 
                                 dplyr::select(Gene,
                                               metric) %>% 
                                 na.omit(),
                               res.females.heat.df %>% 
                                 mutate(metric = sign(log2FoldChange) * -log10(pvalue)) %>% 
                                 dplyr::select(Gene,
                                               metric) %>% 
                                 na.omit(), 
                                labels = c("Male Heat", 
                                           "Female Heat"),
                               boundary = 0.1)


# get data frame of quadrant genes 
rrho_output.df = data.frame(Gene = rrho_output$genelist_uu$gene_list_overlap_uu,
                            Direction = 'up.up') %>% 
  full_join(data.frame(Gene = rrho_output$genelist_dd$gene_list_overlap_dd,
                       Direction = 'down.down'))%>% 
  full_join(data.frame(Gene = c(rrho_output$genelist_ud$gene_list_overlap_ud,
                                NA),
                       Direction = 'up.down'))%>% 
  full_join(data.frame(Gene = c(rrho_output$genelist_du$gene_list_overlap_du,
                                NA),
                       Direction = 'down.up')) %>% 
  na.omit()

## graph
# RRHO plot
png('figures/DEG/RRHO heat by sex.png')
RRHO2_heatmap(rrho_output)
dev.off()

# paper
pdf('figures/paper/DEG_linger_heat_sex_3.pdf',
    height = 6.5,
    width = 6.5)
RRHO2_heatmap(rrho_output)
dev.off()


# gene categories
rrho_output.df %>% 
  count(Direction) %>% 
  ggplot(aes(x = Direction,
             y = n,
             label = n)) +
  geom_label() +
  theme_classic() +
  ylab('Number of genes') +
  xlab('RRHO category')
ggsave('figures/DEG/RRHO genes per category.png') 

### graph logfc across sexes
# get max and min
max = max(abs(res.males.heat.df$log2FoldChange),
          abs(res.females.heat.df$log2FoldChange))

## graph DEG for each sex
# all DEG
res.males.heat.df %>%
  dplyr::select(Gene,
                log2FoldChange) %>%
  dplyr::rename(male_fc = log2FoldChange) %>%
  inner_join(res.females.heat.df %>%
               dplyr::select(Gene, log2FoldChange) %>%
               dplyr::rename(female_fc = log2FoldChange),
             by = "Gene") %>%
  mutate(Concordance = case_when(sign(male_fc) != sign(female_fc) ~ 'Dis',
                                 TRUE ~ 'Con')) %>% 
  left_join(res.treat.df %>%
              mutate(Sig = ifelse(padj < 0.16,
                                   'Sig',
                                   'Not sig')) %>% 
              dplyr::select(Gene,
                            Sig)) %>%
  mutate(Sig = ifelse(is.na(Sig),
                      'Not sig',
                      Sig),
         label = ifelse(Sig == 'Sig',
                        Gene,
                        NA)) %>% 
  arrange(Sig) %>% 
  ggplot(aes(x = male_fc,
             y = female_fc,
             color = Sig)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = 'dashed') +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = label),
                            color = 'black',
                            force = 15,
                            force_pull = 0,
                            point.padding = 0.5,
                            box.padding = 1) +
  theme_classic() +
  xlim(-max,
       max) +
  ylim(-max,
       max) +
  scale_color_manual(values = c('black',
                                'red'))+
    ggtitle('All DEG')
ggsave('figures/DEG/Interaction scatterplot DEG.png')

# each sex
res.males.heat.df %>%
  mutate(Male.sig = ifelse(padj < 0.05,
                      'Sig',
                      'Not sig')) %>% 
  dplyr::select(Gene,
                log2FoldChange,
                Male.sig) %>%
  dplyr::rename(male_fc = log2FoldChange) %>%
  inner_join(res.females.heat.df %>%
               mutate(Female.sig = ifelse(padj < 0.05,
                                   'Sig',
                                   'Not sig')) %>%
               dplyr::select(Gene,
                             log2FoldChange,
                             Female.sig) %>%
               dplyr::rename(female_fc = log2FoldChange),
             by = "Gene") %>%
  mutate(Concordance = case_when(sign(male_fc) != sign(female_fc) ~ 'Dis',
                                 TRUE ~ 'Con'))  %>%
  mutate(Sig = case_when(Male.sig == 'Sig' & Female.sig != 'Sig' ~ 'Male',
                         Male.sig != 'Sig' & Female.sig == 'Sig' ~ 'Female',
                         Male.sig == 'Sig' & Female.sig == 'Sig' ~ 'Both',
                         TRUE ~ 'Not sig'),
         Sig.order = case_when(Sig == 'Both' ~ 1,
                               Sig == 'Male' ~ 2,
                               Sig == 'Female' ~ 3,
                               Sig == 'Not sig' ~ 0)) %>% 
  arrange(Sig.order) %>% 
  mutate(Label = ifelse(Gene %in% c('LOC120765528',
                                    'RYBP'),
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = male_fc,
             y = female_fc,
             color = Sig)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = 'dashed') +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = Label),
             color = 'black') +
  theme_classic() +
  xlim(-max,
       max) +
  ylim(-max,
       max) + 
  scale_color_manual(values = c('purple',
                                'darkred',
                                'darkblue',
                                'grey'))+
  labs(color = 'Sig (padj<0.05)') +
  ggtitle('Sex specific DEG')
ggsave('figures/DEG/Interaction scatterplot sex specific DEG.png')

# paper
res.males.heat.df =read.csv('data/res.males.heat.df.csv')

res.females.heat.df =read.csv('data/res.females.heat.df.csv')

res.males.heat.df %>%
  mutate(Male.sig = ifelse(padj < 0.05,
                           'Sig',
                           'Not sig')) %>% 
  dplyr::select(Gene,
                log2FoldChange,
                Male.sig) %>%
  dplyr::rename(male_fc = log2FoldChange) %>%
  inner_join(res.females.heat.df %>%
               mutate(Female.sig = ifelse(padj < 0.05,
                                          'Sig',
                                          'Not sig')) %>%
               dplyr::select(Gene,
                             log2FoldChange,
                             Female.sig) %>%
               dplyr::rename(female_fc = log2FoldChange),
             by = "Gene") %>%
  mutate(Concordance = case_when(sign(male_fc) != sign(female_fc) ~ 'Dis',
                                 TRUE ~ 'Con'))  %>%
  mutate(Sig = case_when(Male.sig == 'Sig' & Female.sig != 'Sig' ~ 'Male',
                         Male.sig != 'Sig' & Female.sig == 'Sig' ~ 'Female',
                         Male.sig == 'Sig' & Female.sig == 'Sig' ~ 'Both',
                         TRUE ~ 'Not sig'),
         Sig.order = case_when(Sig == 'Both' ~ 1,
                               Sig == 'Male' ~ 2,
                               Sig == 'Female' ~ 3,
                               Sig == 'Not sig' ~ 0),
         Sig.color = case_when(female_fc > 0 & Sig == 'Both' ~ 'Heat',
                               female_fc > 0 & Sig == 'Male' ~ 'Heat',
                               female_fc > 0 & Sig == 'Female' ~ 'Heat',
                               female_fc < 0 & Sig == 'Both' ~ 'Con',
                               female_fc < 0 & Sig == 'Male' ~ 'Con',
                               female_fc < 0 & Sig == 'Femal' ~ 'Con',
                               TRUE ~ 'none'),
         Sig.shape = case_when(Sig == 'Both' ~ 'Both',
                               Sig == 'Male' | Sig == 'Female' ~ 'none',
                               TRUE ~ 'none')) %>% 
  arrange(desc(Sig.color)) %>% 
  ggplot(aes(x = male_fc,
             y = female_fc,
             color = Sig.color,
             shape = Sig.shape)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_abline(slope = 1,
  #             intercept = 0,
  #             linetype = 'dashed') +
  geom_point(size = 3) +
  # ggrepel::geom_label_repel(aes(label = Label),
  #                           color = 'black') +
  theme_classic() +
  # xlim(-max,
  #      max) +
  # ylim(-max,
  #      max) + 
  scale_color_manual(values = c('#464b9f',
                                '#f26622',
                                'grey'))+
  labs(color = 'Sig (padj<0.05)') +
  ggtitle('Sex specific DEG') +
  coord_fixed(ratio = 2) +
  theme(legend.position = 'none') +
  # theme(legend.position = 'inside',
  #       legend.position.inside = c(0.9,
  #                                  0.2)) +
  xlab('Male log2FC')+
  ylab('Female log2FC')
ggsave('figures/paper/Deg_linger_heat_sex_2a.pdf',
       height = 3.25,
       width = 3.25)

# check overlap 
res.males.heat.df %>%
  mutate(Male.sig = ifelse(padj < 0.05,
                           'Sig',
                           'Not sig')) %>% 
  dplyr::select(Gene,
                log2FoldChange,
                Male.sig) %>%
  dplyr::rename(male_fc = log2FoldChange) %>%
  inner_join(res.females.heat.df %>%
               mutate(Female.sig = ifelse(padj < 0.05,
                                          'Sig',
                                          'Not sig')) %>%
               dplyr::select(Gene,
                             log2FoldChange,
                             Female.sig) %>%
               dplyr::rename(female_fc = log2FoldChange),
             by = "Gene") %>%
  mutate(Concordance = case_when(sign(male_fc) != sign(female_fc) ~ 'Dis',
                                 TRUE ~ 'Con')) %>% 
  mutate(Sig = case_when(Male.sig == 'Sig' & Female.sig != 'Sig' & male_fc > 0 ~ 'Male.up',
                         Male.sig == 'Sig' & Female.sig != 'Sig' & male_fc < 0 ~ 'Male.down',
                         Male.sig != 'Sig' & Female.sig == 'Sig' & female_fc > 0  ~ 'Female.up',
                         Male.sig != 'Sig' & Female.sig == 'Sig' & female_fc < 0  ~ 'Female.down',
                         Male.sig == 'Sig' & Female.sig == 'Sig' & male_fc > 0 ~ 'Both.up',
                         Male.sig == 'Sig' & Female.sig == 'Sig'& male_fc < 0 ~ 'Both.down',
                         TRUE ~ 'Not sig'),
         Sig.order = case_when(Sig == 'Both' ~ 1,
                               Sig == 'Male' ~ 2,
                               Sig == 'Female' ~ 3,
                               Sig == 'Not sig' ~ 0)) %>% 
  full_join(res.treat.df %>% 
              filter(padj < 0.155) %>% 
              mutate(DEG.heat = 'DEG.heat') %>% 
              dplyr::select(Gene,
                            DEG.heat)) %>% 
  count(Sig,
        DEG.heat,
        Concordance)

# Sig DEG.heat Concordance    n
# 1    Both.down DEG.heat         Con    1
# 2    Both.down     <NA>         Con   36
# 3      Both.up     <NA>         Con   14
# 4  Female.down     <NA>         Con    5
# 5    Female.up     <NA>         Con    2
# 6    Male.down     <NA>         Con   30
# 7      Male.up     <NA>         Con    8
# 8      Not sig DEG.heat         Con    4
# 9      Not sig DEG.heat         Dis    1
# 10     Not sig     <NA>         Con 8473
# 11     Not sig     <NA>         Dis  730


## heat boxplots both sexes
# get TOP DEG
# get overlapping DEG
res.sex.both.heat.df.deg.top =
  res.males.heat.df %>% 
  mutate(Sig.male = ifelse(padj <= 0.05,
                           'Sig',
                           'Not sig')) %>% 
  dplyr::select(Sig.male,
                Gene,
                padj) %>% 
  left_join(res.females.heat.df %>% 
              mutate(Sig.female = ifelse(padj <= 0.05,
                                         'Sig',
                                         'Not sig')) %>% 
              dplyr::select(Sig.female,
                            Gene,
                            padj),
            by = 'Gene',
            suffix = c('_male',
                       '_female')) %>% 
  mutate(Sig.both = ifelse(Sig.male == 'Sig' & Sig.female == 'Sig',
                           'Both',
                           NA)) %>% 
  filter(Sig.both == 'Both') %>% 
    dplyr::select(Gene,
                  padj_male,
                  padj_female) %>% 
    pivot_longer(cols = -c('Gene'),
                 names_to = 'Sex',
                 values_to = 'padj') %>% 
    group_by(Sex) %>% 
    slice_min(order_by = padj,
              n = 1,
              with_ties = FALSE) %>% 
    ungroup() %>% 
  pull(Gene) %>% 
  unique()
    
    
    
# normalize and scale counts
deg.sex.heat.mat = assay(vst(dds,
                    blind = F))[res.sex.both.heat.df.deg.top,] %>% 
  t() %>% 
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(-c(Sample),
               names_to = 'Gene',
               values_to = 'z.score') %>% 
  left_join(data.linger.sample.meta) %>% 
  left_join(data.linger.sample.sex)


## graph each DEG boxplots
# treatment temp
# paper
deg.sex.heat.mat %>% 
  mutate(ID = paste0(Treatment,
                     '.',
                     Sex)) %>% 
  ggplot(aes(x = ID,
             y = z.score,
             color = Treatment)) +
  geom_boxplot(aes(group = ID),
               outlier.shape = NA) +
  geom_point(size = 2.5,
             aes(shape = Sex)) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab('Gene expression (z-score)') +
  xlab('') +
  scale_color_manual(values = c('#464b9f',
                                '#f26622')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 2)  +
  theme(legend.position = 'none')
ggsave('figures/paper/DEG_linger_heat_sex_1b.pdf',
       height = 3.25,
       width = 3.25,
       dpi = 720)


#### WGCNA ####
options(stingsAsFactors = F)

### prepare data
# use full model to keep interaction effects
vsd.interact = vst(dds.interact, 
                   blind = FALSE)

# remove low variance genes
keep = rowVars(assay(vsd.interact)) > quantile(rowVars(assay(vsd.interact)), 0.25)

# pivot
datExpr = t(assay(vsd.interact))[, keep]

# check samples
gsg = goodSamplesGenes(datExpr,
                       verbose = 3)

# clean up
datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]

### check softthreshold
## pick soft threshold
# powers
powers = c(1:10, seq(12,20, 2))

# pick threshold
sft = pickSoftThreshold(datExpr,
                        powerVector = powers,
                        networkType = 'signed',
                        verbose = 3)

## graph 
png('figures/WGCNA/Softpower and connectivity threshold.png',
    width = 10,
    height = 5,
    units = 'in',
    res = 720)
par(mfrow = c(1,2))

# Scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.90, col="red") # Threshold line

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

# use softpower threshold of 10

### run WGCNA
net = blockwiseModules(
  datExpr,
  power = 10,
  maxBlockSize = 20000,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3)

## graph WGCNA
# dendrogram from WGCNA result
geneTree = net$dendrograms[[1]]

# module colors
moduleColors <- labels2colors(net$colors)

# plot
png('figures/WGCNA/WGCNA dendrogram.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 720)
plotDendroAndColors(
  geneTree,
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()


# plot genes per module
table(net$colors) %>% 
  as.data.frame() %>% 
  dplyr::rename('Modules' = 'Var1') %>% 
  mutate(Modules = paste0('ME',
                          Modules)) %>% 
  ggplot(aes(x = reorder(Modules,
                         -Freq),
             y = Freq,
             label = Freq)) +
  geom_point() +
  ggrepel::geom_label_repel() +
  theme_classic() +
  xlab('Modules')
ggsave('figures/WGCNA/Genes per module.png')


## extract module eigengenes
# set colors
# get module eigengenes
MEs = orderMEs(net$MEs)

## graph ME 
# create dendrogram
METree <- hclust(as.dist(1 - cor(MEs)), method = "average")

# graph correlation of modules
png('figures/WGCNA/Module correlations.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 720)
plot(METree, main = "Clustering of module eigengenes")

plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

### test ME vs treatment and sex
## create list of MEs with meta data
MEs.sample = MEs %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(cols = -c('Sample'),
               names_to = 'Module',
               values_to = 'MEs') %>% 
  left_join(data.linger.sample %>% 
              rownames_to_column('Sample'))  %>%
  mutate(ID = interaction(Sex, Treatment))

## statistics
# run linear model 
results_list <- lapply(colnames(MEs), function(me) {
  fit <- lm(MEs ~ Treatment * Sex, 
            data = MEs.sample %>% 
              filter(Module == me))
  summary(fit)$coefficients
})

names(results_list) <- colnames(MEs)

# create data frame
results_df <- imap(results_list, function(mat, module) {
  as.data.frame(mat) %>%
    tibble::rownames_to_column("term") %>%
    mutate(module = module)
}) %>%
  bind_rows() %>%
  filter(term != "(Intercept)") %>%
  select(module, 
         term, 
         Estimate, 
         `Pr(>|t|)`) %>%
  pivot_wider(names_from = term,
              values_from = c(Estimate, 
                              `Pr(>|t|)`))

# clean column name
colnames(results_df) = gsub("Pr\\(>\\|t\\|\\)", 
                             "p",
                             colnames(results_df))

colnames(results_df) = gsub(":", 
                            "_",
                            colnames(results_df))


# add fdr correction
results_df = results_df%>%
  mutate(
    padj_Treatment  = p.adjust(p_TreatmentHot, 
                               method = "fdr"),
    padj_Sex        = p.adjust(p_SexMale,
                               method = "fdr"),
    padj_Interaction = p.adjust(p_TreatmentHot_SexMale, 
                                method = "fdr")
  )

# select signficant modules
results_df_sig_modules = results_df %>% 
  filter(padj_Treatment < 0.1 | padj_Sex < 0.1 | padj_Interaction < 0.1 ) %>% 
  pull(module)

### graph boxplots of ME
## all modules 
MEs.sample %>%
  ggplot(aes(x = ID,
             y = MEs,
             group = ID)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ID)) + 
  facet_wrap(.~Module,
             scales = 'free') +
  theme_classic()
ggsave('figures/WGCNA/Boxplot MEs by sex and treatment.png',
       width = 20,
       height = 20,
       units = 'in')

## select modules to plot
# function to extract letters per module
get_letters <- function(df) {
  fit <- lm(MEs ~ ID, data = df)
  em <- emmeans::emmeans(fit, ~ ID)
  cld_res <- multcomp::cld(em, Letters = letters)
  
  cld_res %>%
    as.data.frame() %>%
    select(ID, .group) %>%
    rename(letter = .group)
}

# get letters 
letters_df <- MEs.sample %>% 
  filter(Module %in% results_df_sig_modules) %>% 
  group_by(Module) %>%
  group_modify(~ get_letters(.x))

# get y value for letters
y_pos <- MEs.sample %>%
  filter(Module %in% results_df_sig_modules) %>% 
  group_by(Module, ID) %>%
  summarise(y = max(MEs), .groups = "drop")

letters_df <- left_join(letters_df, y_pos,
                        by = c("Module","ID"))

# filter down to modules of interest
MEs.sample %>%
  filter(Module %in% results_df_sig_modules) %>% 
  ggplot(aes(x = ID,
             y = MEs,
             group = ID)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ID)) + 
  geom_text(data = letters_df,
            aes(x = ID, 
                y = y + 0.2, 
                label = letter),
            inherit.aes = FALSE,
            size = 5) +
  facet_wrap(.~Module,
             scales = 'free') +
  theme_classic()
ggsave('figures/WGCNA/Boxplot MEs by sex and treatment letters.png',
       width = 10,
       height = 10,
       units = 'in')

## check overlap of modules with Z and W genes
net$colors %>% 
  as.data.frame() %>% 
  dplyr::rename('Module' = '.') %>% 
  mutate(Module = paste0('ME',
                         Module)) %>% 
  rownames_to_column('Gene') %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Module,
        Gene_category) %>% 
  group_by(Module) %>% 
  mutate(total = sum(n)) %>% 
  ungroup() %>% 
  mutate(percent = 100*n/total) %>% 
  ggplot(aes(x = reorder(Module,
                         -total),
             y = percent,
             fill = Gene_category)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  xlab('Module') +
  ylab('Percent module membership') +
  theme(axis.text = element_text(angle = 45,
                                 hjust = 0.5,
                                 vjust = 0))
ggsave('figures/WGCNA/Module gene category percentage.png',
       width = 10,
       height = 5,
       units = 'in')






#### WGCNA batch corrected ####
### remove batch effects of sex

options(stingsAsFactors = F)

### prepare data
# use full model to keep interaction effects
vsd.interact = vst(dds.interact, 
                   blind = FALSE)

# remove batch effects
# use design to keep treatment/interaction
mat.resid = limma::removeBatchEffect(assay(vsd.interact),
                              batch = vsd.interact$Sex,
                              design = model.matrix(~ Treatment, data = colData(vsd.interact)))


# remove low variance genes
keep = rowVars(mat.resid) > quantile(rowVars(mat.resid), 0.25)
datExpr = t(mat.resid)[, keep]

# check samples
gsg = goodSamplesGenes(datExpr,
                       verbose = 3)

# clean up
datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]

### check softthreshold
## pick soft threshold
# powers
powers = c(1:10, seq(12,20, 2))

# pick threshold
sft = pickSoftThreshold(datExpr,
                        powerVector = powers,
                        networkType = 'signed',
                        verbose = 3)

## graph 
png('figures/WGCNA_batch/Softpower and connectivity threshold.png',
    width = 10,
    height = 5,
    units = 'in',
    res = 720)
par(mfrow = c(1,2))

# Scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.80, col="red") # Threshold line

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
abline(h=100, col="red") 
dev.off()

# use softpower threshold of 12

### run WGCNA
net = blockwiseModules(
  datExpr,
  power = 12,
  maxBlockSize = 20000,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3)

## graph WGCNA
# dendrogram from WGCNA result
geneTree = net$dendrograms[[1]]

# module colors
moduleColors <- labels2colors(net$colors)

# plot
png('figures/WGCNA_batch/WGCNA dendrogram.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 720)
plotDendroAndColors(
  geneTree,
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()


# plot genes per module
table(net$colors) %>% 
  as.data.frame() %>% 
  dplyr::rename('Modules' = 'Var1') %>% 
  mutate(Modules = paste0('ME',
                          Modules)) %>% 
  ggplot(aes(x = reorder(Modules,
                         -Freq),
             y = Freq,
             label = Freq)) +
  geom_point() +
  ggrepel::geom_label_repel() +
  theme_classic() +
  xlab('Modules')
ggsave('figures/WGCNA_batch//Genes per module.png')


## extract module eigengenes
# set colors
# get module eigengenes
MEs = orderMEs(net$MEs)

## graph ME 
# create dendrogram
METree <- hclust(as.dist(1 - cor(MEs)), method = "average")

# graph correlation of modules
png('figures/WGCNA_batch//Module correlations.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 720)
plot(METree, main = "Clustering of module eigengenes")

plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

### test ME vs treatment and sex
## create list of MEs with meta data
MEs.sample = MEs %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(cols = -c('Sample'),
               names_to = 'Module',
               values_to = 'MEs') %>% 
  left_join(data.linger.sample %>% 
              rownames_to_column('Sample'))  %>%
  mutate(ID = interaction(Sex, Treatment))

## statistics
# run linear model 
results_list <- lapply(colnames(MEs), function(me) {
  fit <- lm(MEs ~ Treatment * Sex, 
            data = MEs.sample %>% 
              filter(Module == me))
  summary(fit)$coefficients
})

names(results_list) <- colnames(MEs)

# create data frame
results_df <- imap(results_list, function(mat, module) {
  as.data.frame(mat) %>%
    tibble::rownames_to_column("term") %>%
    mutate(module = module)
}) %>%
  bind_rows() %>%
  filter(term != "(Intercept)") %>%
  select(module, 
         term, 
         Estimate, 
         `Pr(>|t|)`) %>%
  pivot_wider(names_from = term,
              values_from = c(Estimate, 
                              `Pr(>|t|)`))

# clean column name
colnames(results_df) = gsub("Pr\\(>\\|t\\|\\)", 
                            "p",
                            colnames(results_df))

colnames(results_df) = gsub(":", 
                            "_",
                            colnames(results_df))


# add fdr correction
results_df = results_df%>%
  mutate(
    padj_Treatment  = p.adjust(p_TreatmentHot, 
                               method = "fdr"),
    padj_Sex        = p.adjust(p_SexMale,
                               method = "fdr"),
    padj_Interaction = p.adjust(p_TreatmentHot_SexMale, 
                                method = "fdr")
  )

# select signficant modules
results_df_sig_modules = results_df %>% 
  filter(padj_Treatment < 0.1 | padj_Sex < 0.1 | padj_Interaction < 0.1 ) %>% 
  pull(module)

### graph boxplots of ME
## all modules 
MEs.sample %>%
  ggplot(aes(x = ID,
             y = MEs,
             group = ID)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ID)) + 
  facet_wrap(.~Module,
             scales = 'free') +
  theme_classic()+
  xlab('') +
  theme(axis.text = element_text(angle = 45,
                                 hjust = 1))
ggsave('figures/WGCNA_batch//Boxplot MEs by sex and treatment.png',
       width = 20,
       height = 20,
       units = 'in')

## select modules to plot
# function to extract letters per module
get_letters <- function(df) {
  fit <- lm(MEs ~ ID, data = df)
  em <- emmeans::emmeans(fit, ~ ID)
  cld_res <- multcomp::cld(em, Letters = letters)
  
  cld_res %>%
    as.data.frame() %>%
    select(ID, .group) %>%
    rename(letter = .group)
}

# get letters 
letters_df <- MEs.sample %>% 
  filter(Module %in% results_df_sig_modules) %>% 
  group_by(Module) %>%
  group_modify(~ get_letters(.x))

# get y value for letters
y_pos <- MEs.sample %>%
  filter(Module %in% results_df_sig_modules) %>% 
  group_by(Module, ID) %>%
  summarise(y = max(MEs), .groups = "drop")

letters_df <- left_join(letters_df, y_pos,
                        by = c("Module","ID"))

# filter down to modules of interest
MEs.sample %>%
  filter(Module %in% results_df_sig_modules) %>% 
  ggplot(aes(x = ID,
             y = MEs,
             group = ID)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ID)) + 
  geom_text(data = letters_df,
            aes(x = ID, 
                y = y + 0.2, 
                label = letter),
            inherit.aes = FALSE,
            size = 5) +
  facet_wrap(.~Module,
             scales = 'free') +
  theme_classic() +
  xlab('') +
  theme(axis.text = element_text(angle = 45,
                                 hjust = 1))
ggsave('figures/WGCNA_batch//Boxplot MEs by sex and treatment letters.png',
       width = 10,
       height = 10,
       units = 'in')

## check overlap of modules with Z and W genes
net$colors %>% 
  as.data.frame() %>% 
  dplyr::rename('Module' = '.') %>% 
  mutate(Module = paste0('ME',
                         Module)) %>% 
  rownames_to_column('Gene') %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Module,
        Gene_category) %>% 
  group_by(Module) %>% 
  mutate(total = sum(n)) %>% 
  ungroup() %>% 
  mutate(percent = 100*n/total) %>% 
  ggplot(aes(x = reorder(Module,
                         -total),
             y = percent,
             fill = Gene_category)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  xlab('Module') +
  ylab('Percent module membership') +
  theme(axis.text = element_text(angle = 45,
                                 hjust = 0.5,
                                 vjust = 0))
ggsave('figures/WGCNA_batch//Module gene category percentage.png',
       width = 10,
       height = 5,
       units = 'in')







#### During: graph to check sex ####
## QC
# check proportion of total genes
# paper
data.during.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  group_by(Gene_category,
           Sample) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(
    Percent = Sum.counts/sum(Sum.counts)) %>% 
  mutate(Sex = ifelse(Gene_category %in% c('Z','W'),
                      'Sex',
                      'Auto'),
         Sex.sample = case_when(Gene_category == 'Z' ~ Percent,
                                TRUE ~ 0),
         Sex.sample = max(Sex.sample)) %>%
  filter(Gene_category %in% c('globin',
                              'Auto')) %>% View()
  ggplot(aes(x = Gene_category,
             y = Percent,
             fill = Gene_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  xlab('') +
  ylab('Percent of total reads') +
  scale_fill_manual(values = c('white',
                               'grey')) +
  ylim(0,1) +
  theme(legend.position = 'none')
ggsave('figures/paper/QC_during_a.pdf',
       width = 3,
       height = 3)


### graph sex
## use total gene counts
# check proportion of total genes 
data.during.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  group_by(Gene_category,
           Sample) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(
    Percent = Sum.counts/sum(Sum.counts)) %>% 
  mutate(Sex = ifelse(Gene_category %in% c('Z','W'),
                      'Sex',
                      'Auto'),
         Sex.sample = case_when(Gene_category == 'Z' ~ Percent,
                                TRUE ~ 0),
         Sex.sample = max(Sex.sample)) %>%
  ggplot(aes(x = reorder(Sample,
                         Sex.sample),
             y = Percent,
             fill = Gene_category)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1)) +
  facet_grid(Sex~.,
             scales = 'free') +
  xlab('') +
  ylab('Percent of total reads')
ggsave('figures/QC_during/Gene category percent of reads.png',
       width = 20,
       height = 10)


# compare Z counts with chromosome
data.during.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto) %>% 
  ggplot(aes(x = Z_A_ratio,
             y = W)) +
  geom_vline(xintercept = 1,
             linetype = 'dashed') +
  geom_vline(xintercept = 0.5,
             linetype = 'dashed') +
  geom_point() +
  theme_classic() +
  xlab('Z to Chr 5 ratio') +
  ylab('Total W reads') 
ggsave('figures/QC_during/Sex Z_A ratio vs W reads.png')

# add label
data.during.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto) %>% 
  ggplot(aes(x = Z_A_ratio,
             y = W)) +
  geom_vline(xintercept = 1,
             linetype = 'dashed') +
  geom_vline(xintercept = 0.5,
             linetype = 'dashed') +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = Sample),
                            max.overlaps = 20) +
  theme_classic() +
  xlab('Z to Chr 5 ratio') +
  ylab('Total W reads') 
ggsave('figures/QC_during//Sex Z_A ratio vs W reads label.png',
       height = 10,
       width = 10)

# compare Z counts with chromosome
# paper
data.during.gene %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Chromosome %in% c('Z',
                           'W',
                           5)) %>% 
  filter(Gene_category %in% c('Auto',
                              'W',
                              'Z')) %>% 
  group_by(Sample,
           Gene_category) %>% 
  summarise(Sum.counts = sum(Counts)) %>% 
  pivot_wider(names_from = 'Gene_category',
              values_from = 'Sum.counts') %>% 
  mutate(Z_A_ratio = Z/Auto,
         W_A_ratio = W/Auto) %>% 
  ggplot(aes(x = Z_A_ratio,
             y = W_A_ratio)) +
  geom_vline(xintercept = 1,
             linetype = 'dashed') +
  geom_vline(xintercept = 0.5,
             linetype = 'dashed') +
  geom_point(size = 3) +
  theme_classic() +
  xlab('Z to Chr 5 reads ratio') +
  ylab('W to Chr 5 reads ratio') 
ggsave('figures/paper/QC_during_b.pdf',
       height = 3,
       width = 3)

#### During: Use DESQ2 ####
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### create DESEQ dataset
dds.during <- DESeqDataSetFromMatrix(countData = data.during.gene.mat,
                              colData = data.during.sample,
                              design = ~ Treatment + Sex)
# check data
dds.during

# remove rows with low counts
# need to have a count of at least 5 and be present in at least half of samples
keep = rowSums(counts(dds.during) >= 5) >= nrow(data.during.sample)/2

# remove low expressed genes
dds.during = dds.during[keep,]

#check data
dds.during
# 10240 genes

### run DEG
## run DESeq
dds.during = DESeq(dds.during)

## graph mean variance
# save graph
png('figures/QC_during/Deseq2 mean variance.png')
plotDispEsts(dds.during)
dev.off()


## get results
# heat
res.during.treat = results(dds.during,
                    contrast = c('Treatment',
                                 'Hot',
                                 'Con'))

summary(res.during.treat)

# check direction
png('figures/DEG_during/Heat direction check.png')
plotCounts(dds.during,
           gene=which.min(res.during.treat$padj),
           intgroup="Treatment")
dev.off()

# hsp90
png('figures/DEG_during/Heat direction check HSP90AA1.png')
plotCounts(dds.during,
           gene='HSP90AA1',
           intgroup="Treatment")
dev.off()


# create dataframe
res.during.treat.df = res.during.treat %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.during.treat.df,
          'data/res.during.treat.df.csv',
          row.names = F)

# sex
res.during.sex = results(dds.during,
                  contrast = c('Sex',
                               'Male',
                               'Female'))

summary(res.during.sex)

# check direction
png('figures/DEG_during/Sex direction check.png')
plotCounts(dds.during,
           gene=which.min(res.during.sex$padj),
           intgroup="Sex")
dev.off()

# create dataframe
res.during.sex.df = res.during.sex %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.during.sex.df,
          'data/res.during.sex.df.csv',
          row.names = F)


#### During: graph DEG results ####
### create volcano plot
## heat 
res.during.treat.df %>% 
  mutate(Sig = ifelse(padj <= 0.001,
                      'Sig',
                      'Not sig'),
         Label = ifelse(padj <= 1e-10,
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.001),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Heat DEG (padj < 0.001)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG_during/Heat volcano plot.png',
       height = 10,
       width = 10)  

# paper
res.during.treat.df %>% 
  mutate(Sig = ifelse(padj <= 0.001,
                      'Sig',
                      'Not sig'),
         Label = ifelse(padj <= 1e-10,
                        Gene,
                        NA),
         Sig.color = case_when(Sig == 'Sig' & log2FoldChange > 0 ~ 'Hot',
                               Sig == 'Sig' & log2FoldChange < 0 ~ 'Con',
                               TRUE ~ 'none')) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.001),
             linetype = 'dashed') +
  geom_point(aes(color = Sig.color),
             size = 2) +
  ggrepel::geom_label_repel(aes(label = Label),
                            size = 2) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Heat DEG') +
  scale_color_manual(values = c('#3774ba',
                                '#f9a01b',
                                'grey')) +
  theme(legend.position = 'none')
ggsave('figures/paper/DEG_during_heat_1a.pdf',
       height = 3.25,
       width = 3.25)  

## heat boxplots
## get DEGs
res.during.treat.df.deg = res.during.treat.df %>% 
  slice_min(order_by = padj,
            n = 6) %>% 
  pull(Gene)

# normalize and scale counts
deg.during.mat = assay(vst(dds.during,
                    blind = F))[res.during.treat.df.deg,] %>% 
  t() %>% 
  scale() %>%
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(-c(Sample),
               names_to = 'Gene',
               values_to = 'z.score') %>% 
  left_join(data.during.sample %>% 
              rownames_to_column('Sample'))


## graph each DEG boxplots
# treatment temp
# paper
deg.during.mat %>% 
  mutate(ID = paste0(Treatment,
                     '.',
                     Sex)) %>% 
  ggplot(aes(x = ID,
             y = z.score,
             color = Treatment)) +
  geom_boxplot(aes(group = ID),
               outlier.shape = NA) +
  geom_point(size = 2,
             aes(shape = Sex)) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab('Gene expression (z-score)') +
  xlab('') +
  scale_color_manual(values = c('#3774ba',
                                '#f9a01b')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 2) +
  theme(legend.position = 'none')
ggsave('figures/paper/DEG_during_heat_1b.pdf',
       height = 3.25,
       width = 3.25,
       dpi = 720)

# get chromosome counts
res.during.treat.df %>% 
  mutate(Sig = case_when(padj <= 0.001 & log2FoldChange > 0 ~ 'Hot',
                         padj <= 0.001 & log2FoldChange < 0 ~ 'Con',
                         TRUE ~ 'none')) %>%
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Sig,
        Gene_category)

## Sex 
res.during.sex.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Gene %in% c(res.during.sex.df %>% 
                                      mutate(direction = ifelse(log2FoldChange > 0,
                                                                'up',
                                                                'down')) %>% 
                                      group_by(direction) %>% 
                                      slice_min(order_by = padj,
                                                n = 5) %>% 
                                      pull(Gene)),
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Sex DEG (padj < 0.05)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG_during/Sex volcano plot.png',
       height = 10,
       width = 10)  

# get sex counts
res.during.sex.df %>% 
  mutate(Sig = case_when(padj <= 0.05 & log2FoldChange > 0 ~ 'Male',
                         padj <= 0.05 & log2FoldChange < 0 ~ 'Female',
                         TRUE ~ 'none')) %>%
  count(Sig)

# check chromosome
res.during.sex.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig')) %>%
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Gene_category,
                 shape = Sig),
             size = 5) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Sex DEG (padj < 0.05)') 
ggsave('figures/DEG_during/Sex volcano plot chrom.png',
       height = 10,
       width = 10)  

# paper
res.during.sex.df %>% 
  mutate(Sig = ifelse(padj <= 0.05,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Gene %in% c(res.during.sex.df %>% 
                                      mutate(direction = ifelse(log2FoldChange > 0,
                                                                'up',
                                                                'down')) %>% 
                                      group_by(direction) %>% 
                                      slice_min(order_by = padj,
                                                n = 5) %>% 
                                      pull(Gene)),
                        Gene,
                        NA),
         Label = ifelse(abs(log2FoldChange) > 5.9,
                        Gene,
                        Label)) %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  mutate(Chromosome_type = case_when(Chromosome == 'Z' ~ 'Z',
                                     Chromosome == 'W' ~ 'W',
                                     TRUE ~ 'Auto'),
         Sig = ifelse(is.na(Sig),
                      'Not sig',
                      Sig)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed') +
  geom_point(aes(color = Sig,
                 shape = Chromosome_type),
             size = 3) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Sex DEG (padj < 0.05)') +
  scale_color_manual(values = c('grey',
                                'black')) +
  scale_shape_manual(values = c(16,
                                17,
                                15)) +
  labs(shape = 'Chromosome',
       color = 'Significance') +
  theme(legend.position = 'inside' ,
        legend.position.inside = c(0.9,
                                   0.75))
ggsave('figures/paper/DEG_during_sex_1.pdf',
       height = 6.5,
       width = 6.5)  


# get chromosome counts
res.during.sex.df %>% 
  mutate(Sig = case_when(padj <= 0.05 & log2FoldChange > 0 ~ 'Male',
                         padj <= 0.05 & log2FoldChange < 0 ~ 'Female',
                         TRUE ~ 'none')) %>%
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Sig,
        Gene_category)


#### During: normalize data with vst and graph PCA ####
#VST
vsd.during = vst(dds.during,
          blind = T)


## convert to matrix
vsd.during.mat = assay(vsd.during)

# get top 1000 most variable genes
rv = rowVars(vsd.during.mat)
select = order(rv,
               decreasing = TRUE)[seq_len(min(1000, length(rv)))]

# Run PCA on transposed matrix
pca <- prcomp(t(vsd.during.mat[select, ]))

# Access PC values for each sample
pca_results <- as.data.frame(pca$x)

# the contribution to the total variance for each component
percentVar.tmp <- pca$sdev^2 / sum(pca$sdev^2 )

### compare loadings to chromosome
## check loadings
# top and bottom 10
pca$rotation %>% 
  as.data.frame() %>% 
  slice_max(order_by = PC1,
            n = 10) %>% 
  bind_rows(pca$rotation %>% 
              as.data.frame() %>% 
              slice_min(order_by = PC1,
                        n = 10)) %>% 
  rownames_to_column('Symbol') %>% 
  left_join(data.anno) %>% 
  mutate(sign = sign(PC1)) %>% 
  count(sign,
        Gene_category)

# top 10 overall 
pca$rotation %>% 
  as.data.frame() %>% 
  mutate(abs.pc1 = abs(PC1)) %>% 
  slice_max(order_by = abs.pc1,
            n = 10) %>% 
  rownames_to_column('Symbol') %>% 
  left_join(data.anno) %>% 
  mutate(sign = sign(PC1)) %>% 
  count(sign,
        Gene_category)

#### graph PCA
### graph PCA var
data.frame(percentVar = percentVar.tmp,
           PC = seq(1:length(percentVar.tmp))) %>%
  mutate(PC = paste('PC',
                    PC,
                    sep ='')) %>%
  filter(percentVar >= 0.01) %>%
  ggplot(aes(x = reorder(PC,
                         -percentVar),
             y = percentVar)) +
  geom_point() +
  geom_segment(aes(x=reorder(PC,
                             -percentVar),
                   xend=reorder(PC,
                                -percentVar),
                   y=0,
                   yend=percentVar)) +
  theme_classic() +
  ggtitle('PC variance') +
  xlab('PCs') +
  ylab('Percent variance')
ggsave('figures/QC_during/PCA/PC variance.png')

### plot pca
## PC 1 vs PC 2
# treatment
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.during.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar.tmp[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar.tmp[2]),"% variance")) +
  coord_fixed() +
  theme_classic() 
ggsave('figures/QC_during/PCA/PCA PC1 vs PC2.png')

# paper
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.during.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar.tmp[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar.tmp[2]),"% variance")) +
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('#3774ba',
                                '#f9a01b')) +
  theme(legend.position = 'none')
ggsave('figures/paper/QC_during_2b_nolegend.pdf',
       height = 3,
       width = 3.25)

## PC 3 vs PC 4
# treatment
pca$x %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(data.during.sample %>%
              rownames_to_column('Sample')) %>%
  ggplot(aes(x = PC3,
             y = PC4,
             color = Treatment,
             shape = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC3: ",round(100*percentVar.tmp[3]),"% variance")) +
  ylab(paste0("PC4: ",round(100*percentVar.tmp[4]),"% variance")) +
  coord_fixed() +
  theme_classic() 
ggsave('figures/QC_during/PCA/PCA PC3 vs PC4.png')

### calculate percent variation due to sex and treatment
## use variance partition on normalized 
varPart.during = fitExtractVarPartModel(assay(
  vst(
    dds.during,
    blind = F)), 
  formula = ~ Treatment + Sex, 
  data.during.sample)

# graph global distribution
plotVarPart(sortCols(varPart.during))
ggsave('figures/DEG_during//Variance partition.png')

# average percentage of variation due to treatment across the genome
mean(varPart.during$Treatment) * 100
# 26.83997

# average percentage of variation due to treatment across the genome
mean(varPart.during$Sex) * 100
# 29.80797



#### During: graph VST normalization ####
## plot row sd vs row mean
# normalized
png('figures/QC_during/Vst normalization.png')
meanSdPlot(assay(vsd.during))
dev.off()


## graph expression per sample
# raw
data.frame(assay(dds.during)) %>%
  mutate(
    Gene_id = row.names(dds.during)
  ) %>%
  pivot_longer(-Gene_id) %>%
  ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Expression during",
    x = "sample",
    y = "expression"
  )
ggsave('figures/QC_during/Expression per sample during.png')

# normalized
data.frame(assay(vsd.during)) %>%
  mutate(
    Gene_id = row.names(vsd.during)
  ) %>%
  pivot_longer(-Gene_id) %>%
  ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized Expression during",
    x = "sample",
    y = "normalized expression"
  )
ggsave('figures/QC_during/VST Normalized Expression per sample during.png')

### check data for outliers
### create pvclust matrix for each tissue
## raw expression
# create tmp dataframe
tmp.data = assay(dds.during)
# only use selected genes
tmp.data = tmp.data[select,]

## run pvclust
tmp.p = pvclust::pvclust(tmp.data,
                         method.dist="cor",
                         method.hclust="average",
                         nboot=1000,
                         quiet = F)

# graph
png("figures/QC_during/Sample dendrogam counts.png",
    height = 10,
    width = 10,
    unit = 'in',
    res = 480)
plot(tmp.p)
dev.off()

## normalized expression
# create tmp dataframe
tmp.data = assay(vsd.during)
# only use selected genes
tmp.data = tmp.data[select,]

## run pvclust
tmp.p = pvclust::pvclust(tmp.data,
                         method.dist="cor",
                         method.hclust="average",
                         nboot=1000,
                         quiet = F)
# graph
png("figures/QC_during/Sample dendrogam counts normalized.png",
    height = 10,
    width = 10,
    unit = 'in',
    res = 480)
plot(tmp.p)
dev.off()

# paper
pdf("figures/paper/QC_during_2a.pdf",
    height = 6,
    width = 6.5)
plot(tmp.p)
dev.off()

#### compare across experiments ####
## HSP 
# get list of HSP genes that are DEGs in linger
res.during.treat.df.deg.hsp = res.during.treat.df %>% 
  filter(log2FoldChange > 0) %>% 
  filter(padj < 0.001) %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(Gene_category == 'HSP') %>% 
  pull(Gene)
  
# boxplot
# normalize in each treatment
deg.both.mat.hsp = assay(vst(dds.during,
                           blind = F))[res.during.treat.df.deg.hsp,] %>% 
  t() %>% 
  # scale() %>%
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  pivot_longer(-c(Sample),
               names_to = 'Gene',
               values_to = 'Normalized') %>% 
  left_join(data.during.sample %>% 
              rownames_to_column('Sample')) %>% 
  mutate(Exp = 'during') %>% 
  rbind(assay(vst(dds,
                  blind = F))[res.during.treat.df.deg.hsp,] %>% 
          t() %>% 
          # scale() %>%
          as.data.frame() %>% 
          rownames_to_column('Sample') %>% 
          pivot_longer(-c(Sample),
                       names_to = 'Gene',
                       values_to = 'Normalized') %>% 
          left_join(data.linger.sample %>% 
                      rownames_to_column('Sample')) %>% 
          mutate(Exp = 'linger'))
  

# graph each DEG boxplots
deg.both.mat.hsp %>% 
  # filter(Gene == 'HSP90AA1') %>% 
  mutate(ID = paste0(Treatment,
                     Sex,
                     Exp),
         Treatment.exp = paste0(Exp,
                     Treatment)) %>% 
  ggplot(aes(x = Treatment.exp,
             y = Normalized,
             color = Treatment.exp)) +
  geom_boxplot(aes(group = Treatment.exp),
               outlier.shape = NA) +
  geom_point(size = 2,
             aes(shape = Sex),
             position = position_dodge(width = 0.1)) +
  theme_classic() +
  # theme(axis.text.x = element_blank()) +
  ylab('Normalized expression (log2-scaled)') +
  xlab('') +
  scale_color_manual(values = c('#3774ba',
                                '#f9a01b',
                                '#464b9f',
                                '#f26622')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 4)
ggsave('figures/Comparison/HSP DEGs across experiments.png')

# paper
deg.both.mat.hsp %>% 
  mutate(ID = paste0(Treatment,
                     Sex,
                     Exp),
         Treatment.exp = paste0(Exp,
                                Treatment)) %>% 
  ggplot(aes(x = Treatment.exp,
             y = Normalized,
             color = Treatment.exp)) +
  geom_boxplot(aes(group = Treatment.exp),
               outlier.shape = NA) +
  geom_point(size = 1.7,
             aes(shape = Sex),
             position = position_dodge(width = 0.3)) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab('Normalized expression (log2-scaled)') +
  xlab('') +
  scale_color_manual(values = c('#3774ba',
                                '#f9a01b',
                                '#464b9f',
                                '#f26622')) +
  facet_wrap(~Gene,
             scale = 'free',
             ncol = 4) +
  theme(legend.position = 'none')
ggsave('figures/paper/DEG_during_heat_2.pdf',
       width = 6.5,
       height = 4)

### RRHO2
## create gene lists 
comparison.results.treat.df = res.during.treat.df %>% 
  mutate(metric.during = sign(log2FoldChange) * -log10(pvalue),
         log2FoldChange.during = log2FoldChange,
         sig.during = ifelse(padj < 0.001,
                             'Sig',
                             'Not sig')) %>% 
  dplyr::select(Gene,
                metric.during,
                log2FoldChange.during,
                sig.during) %>% 
  full_join(res.treat.df %>% 
              mutate(metric.linger = sign(log2FoldChange) * -log10(pvalue),
                     log2FoldChange.linger = log2FoldChange,
                     sig.linger = ifelse(padj < 0.2,
                                         'Sig',
                                         'Not sig')) %>% 
              dplyr::select(Gene,
                            metric.linger,
                            log2FoldChange.linger,
                            sig.linger))

# create one without NA for RRHO
comparison.results.treat.df.na = comparison.results.treat.df %>% 
  dplyr::select(Gene,
                metric.during,
                metric.linger) %>% 
  na.omit()

## compare linger vs during response
rrho_output_comp = RRHO2_initialize(comparison.results.treat.df.na %>% 
                                      dplyr::select(Gene,
                                                    metric.during),
                                    comparison.results.treat.df.na %>% 
                                      dplyr::select(Gene,
                                                    metric.linger), 
                               labels = c("During heat", 
                                          "Linger Heat"),
                               boundary = 0.1)


# get data frame of quadrant genes 
rrho_output_comp.df = data.frame(Gene = c(rrho_output_comp$genelist_uu$gene_list_overlap_uu,
                                     NA),
                            Direction = 'up.up') %>% 
  full_join(data.frame(Gene = c(rrho_output_comp$genelist_dd$gene_list_overlap_dd,
                                NA),
                       Direction = 'down.down'))%>% 
  full_join(data.frame(Gene = c(rrho_output_comp$genelist_ud$gene_list_overlap_ud,
                                NA),
                       Direction = 'up.down'))%>% 
  full_join(data.frame(Gene = c(rrho_output_comp$genelist_du$gene_list_overlap_du,
                                NA),
                       Direction = 'down.up')) %>% 
  na.omit()

## graph
# RRHO plot
png('figures/Comparison/RRHO heat by experiment.png')
RRHO2_heatmap(rrho_output_comp)
dev.off()

# paper
pdf('figures/paper/Figure 2a.pdf',
    height = 10,
    width = 10)
RRHO2_heatmap(rrho_output_comp)
dev.off()

# gene categories
rrho_output_comp.df %>% 
  count(Direction) %>% 
  ggplot(aes(x = Direction,
             y = n,
             label = n)) +
  geom_label() +
  theme_classic() +
  ylab('Number of genes') +
  xlab('RRHO category')
ggsave('figures/Comparison/RRHO genes per category.png')

# gene categories
# heat shock
rrho_output_comp.df %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  count(Direction,
        Gene_category) %>% 
  filter(Gene_category == 'HSP') %>% 
  ggplot(aes(x = Direction,
             y = n,
             fill = Gene_category)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  ylab('Number of genes') +
  xlab('RRHO category')
ggsave('figures/Comparison/RRHO genes per category HSP.png')

# create scatterplot
comparison.results.treat.df %>% 
  mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                         sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                         sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                         is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                         sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                         Gene %in% c('LOC120765528','RYBP') ~ 'Linger',
                         TRUE ~ 'Not sig')) %>% 
  left_join(rrho_output_comp.df) %>% 
  mutate(Direction.alpha = ifelse(is.na(Direction),
                                  'none',
                                  'direction')) %>% 
  mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                       0,
                                       log2FoldChange.during),
         log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                        0,
                                        log2FoldChange.linger)) %>% 
  arrange(desc(Sig)) %>%  
  ggplot(aes(x = log2FoldChange.during,
             y = log2FoldChange.linger,
             color = Sig,
             shape = Direction.alpha)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('blue',
                                'red',
                                'grey')) +
  scale_shape_manual(values = c(17,
                                16)) +
  xlab('During log2FC') +
  ylab('Linger log2FC') +
  labs(shape = 'RRHO enriched',
       color = 'Sig DEG') 
ggsave('figures/Comparison/Scatterplot comparison DEG.png')

# get count
comparison.results.treat.df %>% 
  mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                         sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                         sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                         is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                         sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                         Gene %in% c('LOC120765528','RYBP') ~ 'Linger',
                         TRUE ~ 'Not sig')) %>% 
  left_join(rrho_output_comp.df) %>% 
  mutate(Direction.alpha = ifelse(is.na(Direction),
                                  'none',
                                  'direction')) %>% 
  mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                        0,
                                        log2FoldChange.during),
         log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                        0,
                                        log2FoldChange.linger)) %>% 
  count(Sig,
        Direction)

# label concordant genes
comparison.results.treat.df %>% 
  mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                         sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                         sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                         is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                         sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                         TRUE ~ 'Not sig')) %>% 
  left_join(rrho_output_comp.df) %>% 
  mutate(Direction.alpha = ifelse(is.na(Direction),
                                  'none',
                                  'direction')) %>% 
  mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                        0,
                                        log2FoldChange.during),
         log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                        0,
                                        log2FoldChange.linger),
         Label = case_when(Direction == 'up.up' ~ Gene,
                           Direction == 'down.down' ~ Gene,
                           TRUE ~ NA)) %>%
  arrange(desc(Sig)) %>%
  ggplot(aes(x = log2FoldChange.during,
             y = log2FoldChange.linger,
             color = Sig,
             shape = Direction.alpha)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  # ggrepel::geom_label_repel(aes(label = Label),
  #                           color = 'black',
  #                           hjust = 'outward') +
  ggrepel::geom_label_repel(aes(label = Label),
                            color = 'black',
                            position = ggpp::position_nudge_center(x = 2,
                                                                   y = 2,
                                                                   center_x = 0,
                                                                   center_y = 0),
                            box.padding = 1) +
  theme_classic() +
  scale_color_manual(values = c('blue',
                                'red',
                                'grey')) +
  scale_shape_manual(values = c(17,
                                16)) +
  xlab('During log2FC') +
  ylab('Linger log2FC') +
  labs(shape = 'RRHO enriched',
       color = 'Sig DEG') 
ggsave('figures/Comparison/Scatterplot comparison DEG label.png')

# paper
# label concordant genes
comparison.results.treat.df %>% 
  mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                         sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                         sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                         is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                         sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                         TRUE ~ 'Not sig'),
         Sig.color = case_when(Sig == 'During' & log2FoldChange.during > 0 ~ 'Heat, 4-hour',
                               Sig == 'During' & log2FoldChange.during < 0 ~ 'Control, 4-hour',
                               Sig == 'Linger' & log2FoldChange.linger > 0 ~ 'Heat, 24-hour',
                               Sig == 'Linger' & log2FoldChange.linger< 0 ~ 'Control, 24-hour',
                               TRUE ~ NA)) %>% 
  left_join(rrho_output_comp.df) %>% 
  mutate(Direction.alpha = ifelse(is.na(Direction),
                                  'none',
                                  'direction')) %>% 
  mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                        0,
                                        log2FoldChange.during),
         log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                        0,
                                        log2FoldChange.linger),
         Label = case_when(Direction == 'up.up' ~ Gene,
                           Direction == 'down.down' ~ Gene,
                           TRUE ~ NA),
         Sig.color = ifelse(is.na(Label),
                            Sig.color,
                            'Concordant'))  %>% 
  arrange(desc(Sig)) %>%
  ggplot(aes(x = log2FoldChange.during,
             y = log2FoldChange.linger)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color = 'grey',
             size = 2) +
  geom_point(data = comparison.results.treat.df %>% 
               mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                                      sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                                      sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                                      is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                                      sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                                      TRUE ~ 'Not sig'),
                      Sig.color = case_when(Sig == 'During' & log2FoldChange.during > 0 ~ 'Heat, 4-hour',
                                            Sig == 'During' & log2FoldChange.during < 0 ~ 'Control, 4-hour',
                                            Sig == 'Linger' & log2FoldChange.linger > 0 ~ 'Heat, 24-hour',
                                            Sig == 'Linger' & log2FoldChange.linger< 0 ~ 'Control, 24-hour',
                                            TRUE ~ NA)) %>% 
               left_join(rrho_output_comp.df) %>% 
               mutate(Direction.alpha = ifelse(is.na(Direction),
                                               'none',
                                               'direction')) %>% 
               mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                                     0,
                                                     log2FoldChange.during),
                      log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                                     0,
                                                     log2FoldChange.linger),
                      Label = case_when(Direction == 'up.up' ~ Gene,
                                        Direction == 'down.down' ~ Gene,
                                        TRUE ~ NA),
                      Sig.color = ifelse(is.na(Label),
                                         Sig.color,
                                         'Concordant'))  %>% 
               filter(!is.na(Sig.color)) %>% 
               arrange(desc(Sig)),
             aes(color = Sig.color),
             size = 3) +
  # ggrepel::geom_label_repel(aes(label = Label),
  #                           color = 'black',
  #                           hjust = 'outward') +
  ggrepel::geom_label_repel(aes(label = Label),
                            color = 'black',
                            position = ggpp::position_nudge_center(x = 2,
                                                                   y = 2,
                                                                   center_x = 0,
                                                                   center_y = 0),
                            box.padding = 1) +
  theme_classic() +
  scale_color_manual(values = c('black',
                                '#464b9f',
                                '#3774ba',
                                '#f26622',
                                '#f9a01b')) +
  xlab('4-hour') +
  ylab('24-hour') +
  labs(color = 'Sig DEG') +
  theme(legend.position = 'none')
ggsave('figures/paper/Figure 2b.pdf',
       height = 6,
       width = 6)

# label HSP genes
comparison.results.treat.df %>% 
  mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                         sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                         sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                         is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                         sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                         TRUE ~ 'Not sig')) %>% 
  left_join(rrho_output_comp.df) %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  mutate(Direction.alpha = ifelse(is.na(Direction),
                                  'none',
                                  'direction')) %>% 
  mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                        0,
                                        log2FoldChange.during),
         log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                        0,
                                        log2FoldChange.linger),
         Label = case_when(Gene_category == 'HSP' & Direction.alpha == 'direction' ~ Gene,
                           TRUE ~ NA)) %>% 
  arrange(desc(Sig)) %>% 
  ggplot(aes(x = log2FoldChange.during,
             y = log2FoldChange.linger,
             color = Sig,
             shape = Direction.alpha)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  # ggrepel::geom_label_repel(aes(label = Label),
  #                           color = 'black',
  #                           hjust = 'outward') +
  ggrepel::geom_label_repel(aes(label = Label,
                                color = Sig),
                            position = ggpp::position_nudge_center(x = 2,
                                                                   y = 2,
                                                                   center_x = 0,
                                                                   center_y = 0),
                            box.padding = 1,
                            max.overlaps = 100) +
  theme_classic() +
  scale_color_manual(values = c('blue',
                                'red',
                                'grey')) +
  scale_shape_manual(values = c(17,
                                16)) +
  xlab('During log2FC') +
  ylab('Linger log2FC') +
  labs(shape = 'RRHO enriched',
       color = 'Sig DEG') 
ggsave('figures/Comparison/Scatterplot comparison DEG label HSP.png')

# count HSP genes
comparison.results.treat.df %>% 
  mutate(Sig = case_when(sig.during == 'Sig' & sig.linger != 'Sig' ~ 'During',
                         sig.during == 'Sig' & is.na(sig.linger) ~ 'During',
                         sig.during != 'Sig' & sig.linger == 'Sig' ~ 'Linger',
                         is.na(sig.during) & sig.linger == 'Sig' ~ 'Linger',
                         sig.during == 'Sig' & sig.linger == 'Sig' ~ 'Both',
                         TRUE ~ 'Not sig')) %>% 
  left_join(rrho_output_comp.df) %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  mutate(Direction.alpha = ifelse(is.na(Direction),
                                  'none',
                                  'direction')) %>% 
  mutate(log2FoldChange.during = ifelse(is.na(log2FoldChange.during),
                                        0,
                                        log2FoldChange.during),
         log2FoldChange.linger = ifelse(is.na(log2FoldChange.linger),
                                        0,
                                        log2FoldChange.linger),
         Label = case_when(Gene_category == 'HSP' & Direction.alpha == 'direction' ~ Gene,
                           TRUE ~ NA)) %>% 
  arrange(desc(Sig)) %>% 
  filter(Gene_category == 'HSP') %>% 
  count(Gene_category,
        Sig,
        Direction)

## exploratory concordance
comparison.results.treat.df %>% 
  left_join(rrho_output_comp.df) %>% 
  filter(Direction == 'up.up' | Direction == 'down.dwon') %>% 
  mutate(pvalue.during = 10^metric.during,
         pvalue.linger = 10^metric.linger) %>% 
  View()

#### Use DESQ2 on nest temp ####
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#### trial temp

### create DESEQ dataset
dds.temp <- DESeqDataSetFromMatrix(countData = data.linger.gene.mat,
                              colData = data.linger.sample %>% 
                                rownames_to_column('Sample') %>% 
                                left_join(data.linger.sample.meta) %>% 
                                mutate(Avg.Trial.Temp.scale = scale(Avg.Trial.Temp)),
                              design = ~ Avg.Trial.Temp.scale + Sex)

# make sure it is numeric
colData(dds.temp)$Avg.Trial.Temp.scale <- as.numeric(colData(dds.temp)$Avg.Trial.Temp.scale)

# check data
dds.temp

# remove rows with low counts
# need to have a count of at least 5 and be present in at least half of samples
keep = rowSums(counts(dds.temp) >= 5) >= nrow(data.linger.sample)/2

# remove low expressed genes
dds.temp = dds.temp[keep,]

#check data
dds.temp
# 9304 genes

### run DEG
## run DESeq
dds.temp = DESeq(dds.temp)

## get results
# temperature
res.temp = results(dds.temp,
                    name = 'Avg.Trial.Temp.scale')

summary(res.temp)

# create dataframe
res.temp.df = res.temp %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.temp.df,
          'data/res.temp.df.csv',
          row.names = F)




### graph DEG results
### create volcano plot
## temp 
res.temp.df %>% 
  mutate(Sig = ifelse(padj <= 0.11,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Sig == 'Sig',
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.11),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Trial day DEG (padj < 0.11)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG_trial_temp/Trial day temp volcano plot.png',
       height = 10,
       width = 10)  

# temp 
res.temp.df %>% 
  mutate(Sig = ifelse(padj <= 0.11,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Sig == 'Sig',
                        Gene,
                        NA),
         Sig.color = case_when(log2FoldChange > 0 & Sig == 'Sig' ~ 'Heat',
                               log2FoldChange < 0 & Sig == 'Sig' ~ 'Con',
                               TRUE ~ 'none')) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.11),
             linetype = 'dashed') +
  geom_point(aes(color = Sig.color),
             size = 3) +
  ggrepel::geom_label_repel(aes(label = Label),
                            box.padding = 0.5,
                            point.padding = 0.5,
                            size = 2.5,
                            max.overlaps = Inf,
                            nudge_y = 0.25) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Nest temperature DEG') +
  scale_color_manual(values = c('#464b9f',
                                '#f26622',
                                'grey')) +
  theme(legend.position = 'none')
ggsave('figures/paper/DEG_linger_temp_1b.pdf',
       height = 3.25,
       width = 3.25)  

#### next day temp

### create DESEQ dataset
dds.next.temp <- DESeqDataSetFromMatrix(countData = data.linger.gene.mat,
                                   colData = data.linger.sample %>% 
                                     rownames_to_column('Sample') %>% 
                                     left_join(data.linger.sample.meta) %>% 
                                     mutate(Avg.Trial.next.temp.scale = scale(Avg.Trial.Temp.1)),
                                   design = ~ Avg.Trial.next.temp.scale + Sex)

# make sure it is numeric
colData(dds.next.temp)$Avg.Trial.next.temp.scale <- as.numeric(colData(dds.next.temp)$Avg.Trial.next.temp.scale)

# check data
dds.next.temp

# remove rows with low counts
# need to have a count of at least 5 and be present in at least half of samples
keep = rowSums(counts(dds.next.temp) >= 5) >= nrow(data.linger.sample)/2

# remove low expressed genes
dds.next.temp = dds.next.temp[keep,]

#check data
dds.next.temp
# 9304 genes

### run DEG
## run DESeq
dds.next.temp = DESeq(dds.next.temp)

## get results
# temperature
res.next.temp = results(dds.next.temp,
                   name = 'Avg.Trial.next.temp.scale')

summary(res.next.temp)

# create dataframe
res.next.temp.df = res.next.temp %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene')

# save DEG
write.csv(res.next.temp.df,
          'data/res.next.temp.df.csv',
          row.names = F)




### graph DEG results
### create volcano plot
## temp 
res.next.temp.df %>% 
  mutate(Sig = ifelse(padj <= 0.01,
                      'Sig',
                      'Not sig'),
         Label = ifelse(Sig == 'Sig',
                        Gene,
                        NA)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.01),
             linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('Adjusted pvalue (-log10)') +
  ggtitle('Next day temp DEG (padj < 0.01)') +
  scale_color_manual(values = c('grey',
                                'darkred'))
ggsave('figures/DEG_trial_temp/Next day temp volcano plot.png',
       height = 10,
       width = 10)  




#### GO analysis  ####
###  sex bias 
# load data
res.sex.df = read.csv('data/res.sex.df.csv')

## get annotation data
data.anno = read_tsv('data/barn_swallow_ncbi_gene.tsv') %>% 
  mutate(Gene_category = case_when(str_detect(Name,
                                              'globin') ~ 'globin',
                                   str_detect(Name,
                                              'heat shock') ~ 'HSP',
                                   Chromosome == 'Z' ~ 'Z',
                                   Chromosome == 'W' ~ 'W',
                                   is.na(Chromosome) ~ 'Unplaced',
                                   TRUE ~ 'Auto'))

## get list of sex-biased genes
## use 0.05 cutoff and assign bias on logfoldchange
res.sex.df.deg = res.sex.df %>% 
  left_join(data.anno,
            by = c('Gene' = 'Symbol')) %>% 
  filter(padj <= 0.01) %>% 
  filter(Gene_category %in% c('Auto',
                              'HSP',
                              'Unplaced')) %>% 
  mutate(sex.bias = ifelse(log2FoldChange > 0,
                           'Male.biased',
                           'Female.biased'))

# # get count
res.sex.df.deg %>%
  dplyr::count(sex.bias)
# sex.bias   n
# 1 Female.biased 101
# 2   Male.biased 59
  
### over representation analysis
## Males
ego.males = enrichGO(gene          = res.sex.df.deg %>% 
                         filter(sex.bias == 'Male.biased') %>% 
                         pull(Gene),
                       universe      = res.sex.df$Gene,
                       keyType = 'SYMBOL',
                       OrgDb         = org.Hrustica.eg.db,
                       ont           = "ALL",
                       pAdjustMethod = "fdr",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1,
                       minGSSize = 5,
                       maxGSSize = 500,
                       readable      = TRUE)

# # get pairwise terms
# ego.males = ego.males %>% 
#   pairwise_termsim()
# 
# # simplyify 
# ego.males.simple = simplify(ego.males, 
#                               cutoff=0.7, 
#                               by="p.adjust", 
#                               select_fun=min)
# 
# # graph network
# png('figures/DEG/Sex male-biased GO network.png')
# emapplot(ego.males.simple,
#          group_category = TRUE,    
#          node_label = "group")
# dev.off()
# 
# # graph dotplot
# # graph dotplot
# dotplot(ego.males.simple, 
#         showCategory=30, 
#         label_format=NULL) + 
#   ggtitle("Male-biased")+ 
#   facet_wrap(.~ONTOLOGY, 
#              scales = "free_y", 
#              space = "free_y")
# ggsave('figures/DEG/Sex male-biased GO dotplot.png')

## females
ego.females = enrichGO(gene          = res.sex.df.deg %>% 
                       filter(sex.bias == 'Female.biased') %>% 
                       pull(Gene),
                     universe      = res.sex.df$Gene,
                     keyType = 'SYMBOL',
                     OrgDb         = org.Hrustica.eg.db,
                     ont           = "ALL",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.1,
                     minGSSize = 5,
                     maxGSSize = 500,
                     readable      = TRUE)

# get pairwise terms
ego.females = ego.females %>% 
  pairwise_termsim()

# simplyify 
ego.females.simple = simplify(ego.females, 
                               cutoff=0.7, 
                               by="p.adjust", 
                               select_fun=min)

# graph network
png('figures/DEG/Sex female-biased GO network.png')
emapplot(ego.females.simple,
         group_category = TRUE,    
         node_label = "group")
dev.off()

# graph dotplot
dotplot(ego.females.simple, 
        showCategory=30, 
        label_format=NULL) + 
  ggtitle("Female-biased")+ 
  facet_wrap(.~ONTOLOGY, 
             scales = "free_y", 
             space = "free_y")
ggsave('figures/DEG/Sex female-biased GO dotplot.png')

# paper
dotplot(ego.females.simple, 
        showCategory=30, 
        label_format=NULL) + 
  ggtitle("Female-biased")+ 
  facet_wrap(.~ONTOLOGY, 
             scales = "free_y", 
             space = "free_y") +
  theme(axis.text.y = element_text(size = 8))
ggsave('figures/paper/DEG_linger_GO_b.pdf',
       height = 6.5,
       width = 6.5)

### gene set enrichment
# load data
res.treat.df = read.csv('data/res.treat.df.csv')

## heat response
# create gene list from DEG results
res.treat.df.list = res.treat.df %>% 
  mutate(metric = -log10(pvalue) * sign(log2FoldChange)) %>% 
  pull(metric)

# res.treat.df.list = res.treat.df %>% 
#   pull(stat)

# name list
names(res.treat.df.list) = res.treat.df$Gene

# remove NA and sort
res.treat.df.list = res.treat.df.list %>% 
  na.omit() %>% 
  sort(decreasing = T)


## run gene set enrichment
ego.heat = gseGO(geneList = res.treat.df.list,
                       keyType = 'SYMBOL',
                       OrgDb         = org.Hrustica.eg.db,
                       ont           = "ALL",
                       pAdjustMethod = "fdr",
                 pvalueCutoff = 0.1)

# get pairwise terms
ego.heat = ego.heat %>% 
  pairwise_termsim()

# simplyify 
ego.heat.simple = simplify(ego.heat, 
                              cutoff=0.7, 
                              by="p.adjust", 
                              select_fun=min)

# graph network
png('figures/DEG/Heat GO network.png')
emapplot(ego.heat.simple,
         group_category = TRUE,    
         node_label = "group")
dev.off()

# graph dotplot
dotplot(ego.heat.simple, 
        showCategory=30, 
        label_format=NULL) + 
  ggtitle("Heat")+ 
  facet_wrap(.~ONTOLOGY, 
             scales = "free_y", 
             space = "free_y")
ggsave('figures/DEG/Heat GO dotplot.png')

# paper
dotplot(ego.heat.simple, 
        showCategory=30, 
        label_format=NULL) + 
  ggtitle("Heat")+ 
  facet_wrap(.~ONTOLOGY, 
             scales = "free_y", 
             space = "free_y")
ggsave('figures/paper/DEG_linger_GO_a.pdf',
       height = 6.5,
       width = 6.5)

## check leading edge genes 
# get list of genes
# 6 genes
ego.heat.simple.genes = ego.heat.simple %>% 
  as.data.frame() %>% 
  pull(core_enrichment) %>% 
  strsplit(split = "/") %>% 
  unlist() %>% 
  unique()
  
# check gene descriptions
data.anno %>%
  filter(Symbol %in% ego.heat.simple.genes) %>% 
  dplyr::select(Name,
                Symbol)


# graph on volcano plot 
res.treat.df %>% 
  mutate(Sig = ifelse(padj <= 0.153,
                      'DEG',
                      'Not sig'),
         Label = ifelse(Sig == 'DEG',
                        Gene,
                        NA),
         Sig = ifelse(Gene %in% ego.heat.simple.genes,
                      'GSEA',
                      Sig),
         Label = ifelse(Gene %in% ego.heat.simple.genes,
                      Gene,
                      Label)) %>% 
  arrange(desc(Sig)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(pvalue))) +
  # geom_hline(yintercept = -log10(0.15),
  #            linetype = 'dashed') +
  geom_point(aes(color = Sig),
             size = 5) +
  ggrepel::geom_label_repel(aes(label = Label)) +
  theme_classic() +
  ylab('pvalue (-log10)') +
  ggtitle('Heat DEG & GSEA') +
  scale_color_manual(values = c('darkred',
                                'darkblue',
                                'grey'))
ggsave('figures/DEG/Heat GO gsea volcanoplot.png')








# Statistical analysis

## 1. Introduction

During this section, we will cover some common statistical analyses that are conducted in eDNA metabarcoding research. The statistical tests employed can vary widely and will mainly depend on the experimental design and research questions. Due to the breadth of statistical tests and the specificity based on the experimental design, we will be providing more of an introduction to the statistical analysis whereby we will focus on some key concepts, rather than providing a full breakdown of the statistical analysis of metabarcoding data.

Before we start with the statistical analysis, it is essential to keep an overview of the experimental design of our data. For the larger study, we want to examine how fish community composition changes across a latitudinal gradient in Japan, with particular interest in assessing evidence for tropicalisation, i.e., the poleward spread of tropical species due to ocean warming. For this tutorial dataset, we have water samples collected at various sites across Amami island, covering 2 distinct regions (north and south). The north region consists of 6 sampling locations, including Akakina, Asan, Imaizaki, Kurasaki, Sakibaru, and Sani. The south region consists of 5 sampling locations, including Ashiken, Edateku, Heda, Taen, and Yadon. Within each sampling location, 3 water samples were collected.

All of this means that we have 3 eDNA samples per location and 11 locations = 33 samples.

For the statistical comparison, we could investigate differences between all locations with samples as replicates. We could also generate a semi-quantitative dataset through an incidence-based approach by combining water samples within a location and compare the north region vs. the south region.

## 2. Summary stats

### 2.1 Read and ASV numbers

One of the first things we will do is to generate some numbers summarising the data. This information is usually provided in the first couple of paragraphs of the results section in a metabarcoding manuscript. Stats that are frequently included are total read count prior and post quality filtering, the average read number per sample and per group, as well as the minimum and maximum values, the number of ASVs, the most abundant and most frequently detected ASVs, read count and ASV number in negative controls, sample drop out, etc.

Most of these numbers have been reported in the previous sections, so I suggest going back to the respective section to look them up if needed. There are, however, a few additional numbers we can report on now that we have the finalised tables.

First, let's prepare the R environment and read in the necessary data. For this, we will open a new R script and save it as **statistical_analysis.R** in the `bootcamp` directory.

```{code-block} R
## set up the working environment
setwd(file.path(Sys.getenv("HOME"), "Desktop", "bootcamp"))
list.files()

## load necessary R packages
for (pkg in c("tidyverse", "ggplot2", "car", "viridis", "UpSetR")) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
  cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
}

## read data into R
metadata <- read.table(file.path("4.results", "metadata_filtered.txt"), header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
count_table <- read.table(file.path("4.results", "count_table_filtered.txt"), header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
taxonomy <- read.table(file.path("4.results", "taxonomy_filtered.txt"), header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
```

One of the first things we can report on is the total number of reads and number of detected ASVs. We can easily do this by summing the values in the count table, as well as printing the number of rows in the count table, respectively.

```{code-block} R
## total number of reads and ASVs
sum(count_table)
nrow(count_table)
```

````{admonition} Output
:class: note, dropdown
```
> sum(count_table)
[1] 2993896
> nrow(count_table)
[1] 1695
```
````

Here, we can see that the final number of reads included in our count table is 2,993,896 reads and the final number of ASVs is equal to 1,695.

Next, we can determine the average number of reads per sample, as well as the standard error. For interest, we can look into the minimum and maximum values as well.

```{code-block} R
## average number of reads per sample
metadata$final_read_count <- colSums(count_table)[rownames(metadata)]
mean(metadata$final_read_count, na.rm = TRUE)
sd(metadata$final_read_count, na.rm = TRUE) / sqrt(sum(!is.na(metadata$final_read_count)))
min(metadata$final_read_count)
rownames(metadata)[metadata$final_read_count == min(metadata$final_read_count)]
max(metadata$final_read_count)
rownames(metadata)[metadata$final_read_count == max(metadata$final_read_count)]
```

````{admonition} Output
:class: note, dropdown
```
> mean(metadata$final_read_count, na.rm = TRUE)
[1] 90724.12
> sd(metadata$final_read_count, na.rm = TRUE) / sqrt(sum(!is.na(metadata$final_read_count)))
[1] 5027.811
> min(metadata$final_read_count)
[1] 52096
> rownames(metadata)[metadata$final_read_count == min(metadata$final_read_count)]
[1] "Sani_3"
> max(metadata$final_read_count)
[1] 161751
> rownames(metadata)[metadata$final_read_count == max(metadata$final_read_count)]
[1] "Asan_3"
```
````

From the output, we can see that the average number of reads across all samples is 90,724 ± 5,028. The minimum number of reads (52,096) was observed for Sani_3, while the maximum number of reads (161,751) was observed for Asan_3.

Next, we can investigate the same numbers for ASVs as well.

```{code-block} R
## average number of ASVs per sample
metadata$final_asv_count <- colSums(count_table > 0)[rownames(metadata)]
mean(metadata$final_asv_count, na.rm = TRUE)
sd(metadata$final_asv_count, na.rm = TRUE) / sqrt(sum(!is.na(metadata$final_asv_count)))
min(metadata$final_asv_count)
rownames(metadata)[metadata$final_asv_count == min(metadata$final_asv_count)]
max(metadata$final_asv_count)
rownames(metadata)[metadata$final_asv_count == max(metadata$final_asv_count)]
```

````{admonition} Output
:class: note, dropdown
```
> mean(metadata$final_asv_count, na.rm = TRUE)
[1] 359.0303
> sd(metadata$final_asv_count, na.rm = TRUE) / sqrt(sum(!is.na(metadata$final_asv_count)))
[1] 17.35287
> min(metadata$final_asv_count)
[1] 188
> rownames(metadata)[metadata$final_asv_count == min(metadata$final_asv_count)]
[1] "Akakina_3"
> max(metadata$final_asv_count)
[1] 521
> rownames(metadata)[metadata$final_asv_count == max(metadata$final_asv_count)]
[1] "Heda_3"
```
````

Here, we can see that the average number of ASVs per sample is 359 ± 17. The lowest number of ASVs (188) was observed in Akakina_3, while the maximum number of ASVs (521) was observed for Heda_3.

Besides read and ASV count per sample, it might also be of interest to group the data per sampling location.

```{code-block} R
## average number of reads and ASVs per location
location_info <- metadata %>%
  group_by(sampleLocation) %>%
  summarise(mean_reads = mean(final_read_count, na.rm = TRUE), se_reads = sd(final_read_count, na.rm = TRUE) / sqrt(n()),
            mean_asvs = mean(final_asv_count, na.rm = TRUE), se_asvs = sd(final_asv_count, na.rm = TRUE) / sqrt(n()))

# reads info
min(location_info$mean_reads)
location_info$sampleLocation[location_info$mean_reads == min(location_info$mean_reads)]
max(location_info$mean_reads)
location_info$sampleLocation[location_info$mean_reads == max(location_info$mean_reads)]

# ASV info
min(location_info$mean_asvs)
location_info$sampleLocation[location_info$mean_asvs == min(location_info$mean_asvs)]
max(location_info$mean_asvs)
location_info$sampleLocation[location_info$mean_asvs == max(location_info$mean_asvs)]
```

````{admonition} Output
:class: note, dropdown
```
> min(location_info$mean_reads)
[1] 63703
> location_info$sampleLocation[location_info$mean_reads == min(location_info$mean_reads)]
[1] "Sani"
> max(location_info$mean_reads)
[1] 139100.3
> location_info$sampleLocation[location_info$mean_reads == max(location_info$mean_reads)]
[1] "Asan"
> # ASV info
> min(location_info$mean_asvs)
[1] 216
> location_info$sampleLocation[location_info$mean_asvs == min(location_info$mean_asvs)]
[1] "Akakina"
> max(location_info$mean_asvs)
[1] 477.3333
> location_info$sampleLocation[location_info$mean_asvs == max(location_info$mean_asvs)]
[1] "Sani"
```
````

The output here shows us that the lowest average number of reads (63,703) and ASVs (216) were observed in Sani and Akakina, respectively, while the highest average number of reads (139,100) and ASVs (477) were observed in Asan and Sani, respectively.

### 2.2 ASV exploration

Besides read count and ASV count per sample, we can also investigate the most abundant (with regards to read count) and most frequently detected (occurrence count) ASVs. First, let's add those data (total read count and occurrence count) into new columns in the taxonomy dataframe.

```{code-block} R
## ASV exploration
# add total read count and occurrence count to taxonomy
taxonomy$final_read_count <- rowSums(count_table)[rownames(taxonomy)]
taxonomy$pct_reads <- round(taxonomy$final_read_count / sum(count_table) * 100, 2)
taxonomy$occurrence_count <- rowSums(count_table > 0)[rownames(taxonomy)]
taxonomy$pct_occurrence <- round(taxonomy$occurrence_count / ncol(count_table) * 100, 2)
```

Now, let's investigate the 10 most abundant and frequently detected ASVs by subsetting the taxonomy dataframe

```{code-block} R
# top 10 abundant ASVs
abundant_asvs <- taxonomy %>%
  slice_max(final_read_count, n = 10, with_ties = TRUE)

# top 10 frequent ASVs
frequent_asvs <- taxonomy %>%
  slice_max(occurrence_count, n = 10, with_ties = TRUE)
```

When investigating these dataframes, we can see that our data is dominated by ASV_1 (974,368 reads; 32.55%) and ASV_2 (842,055 reads; 28.13%) with regards to read count, taking up a total of 1,816,423 (60.68%) reads. The remaining 8 most abundant ASVs (ASV_3, ASV_8, ASV_5, ASV_7, ASV_14, ASV_9, ASV_15, and ASV_41) only take up a combined 368,573 reads or 12.31% within our data. When looking at the taxonomy assignment for ASV_1 and ASV_2, we can see that they are assigned to two different species within the same genus, i.e., *Micromonas commoda* and *Micromonas pusilla*. Micromonas is a widespread, very small-sized, prasinophyta alga that makes up a significant amount of picoplanktonic biomass and productivity in both oceanic and coastal regions. Hence, the abundance of these signals is not unsurprising, given live Micromonas organisms were likely captured on the filter, extracted, amplified using a universal COI primer set, and sequenced. The high abundance of picoplankton in COI metabarcoding data from water samples has been observed frequently and is one of the drawbacks of using a universal primer set that co-amplifies these organisms.

For the most frequently detected ASVs, we can see from the dataframe that we have 26 (1.53%) ASVs that were detected in all samples, including the two most abundant ASVs (ASV_1 and ASV_2). When looking at the taxonomic classifications of these 26 ASVs, we can see that, besides the two Micromonas species, we also detected an ascidian (*Ascidia ahodori*), and two other phytoplankton species (*Phaeocystis globosa* and *Pelagomonas calceolata*). The remaining ASVs only have poorly supported taxonomic classifications (i.e., BLAST percent identity and coverage values below 90) and, hence, are not reliable.

### 2.3 Taxa exploration

When using universal primer sets for metabarcoding research, it is a common occurrence that some percentage of ASVs will be unclassified or receive a low-quality BLAST-hit. To investigate the quality of the BLAST-hits we received for our ASVs in this tutorial data, let's plot the distribution of the percent identity values (note: BLAST is a local alignment, for which we set a threshold of query coverage to 50%, it might be needed to use global alignment values for a more accurate representation of this plot).

```{code-block} R
## Taxa exploration
# generate density plot for BLAST-hits
ggplot(taxonomy, aes(x = pident)) + 
  geom_density(fill = "#377EB8", alpha = 0.3, color = "#377EB8", linewidth = 0.5, adjust = 1, na.rm = TRUE) +
  geom_vline(aes(xintercept = 100), color = "gray40", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = 95), color = "#806491", linetype = "dashed", linewidth = 0.5) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        plot.title = element_text(hjust = 0.5, face = "bold"), axis.line = element_line(color = "black")) + 
  labs(title = "Density Distribution of BLAST-hits", x = "Percent Identity", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  annotate("text", x = 100, y = Inf, vjust = 1, hjust = 1.05,
           label = paste0("Species cutoff (100%):\n", sum(taxonomy$pident >= 100, na.rm = TRUE), " (", 
                          round(sum(taxonomy$pident >= 100, na.rm = TRUE) / nrow(taxonomy) * 100, 2), "%)"), color = "gray40") +
  annotate("text", x = 95, y = Inf, vjust = 2, hjust = 1.05,
           label = paste0("BLAST cutoff (95%):\n", sum(taxonomy$pident >= 95, na.rm = TRUE), " (",
                          round(sum(taxonomy$pident >= 95, na.rm = TRUE) / nrow(taxonomy) * 100, 2), "%)"), color = "#806491")
ggsave(filename = file.path("4.results", "blast_hit_distribution.png"), plot = last_plot(), bg = "white", dpi = 300)
```

```{figure} blast_hit_distribution.png
:name: blast hit distribution

: The distribution of the percentage identity values of the BLAST hits.
```

The figure shows us that we have 123 (7.26%) ASVs with a 100% percent identity BLAST hit and, therefore, could be potentially assigned to species level (though query coverage will need to be taken into account for this as well). Furthermore, when we set a cutoff threshold of 95% for percent identity of BLAST hits, we see that we can assign a taxonomic ID to 301 or 17.76% of ASVs.

Let's explore these high quality BLAST-hits a bit further. For this tutorial, we will subset the data for ASVs reaching a percent identity BLAST-hit of 100% and query coverage of at least 98%.

```{code-block} R
# subset high quality BLAST-hits
taxonomy_species <- taxonomy %>%
  filter(pident >= 100, qcov >= 98)
```

When subsetting the data based on percent identity and query coverage, we see that only 117 ASVs are included. We can now subset the count table as well, and reformat it, so that we can plot a simple stacked bar plot containing the read abundance of each species. Since we have over 100 taxa, it will be better to colour the bars based on a higher taxonomic level (we'll use phylum in our case), and reorder or sort the bars based on the total number of reads to that higher taxonomic level.

```{code-block} R
count_species <- count_table %>%
  filter(rownames(.) %in% rownames(taxonomy_species)) %>%
  mutate(species = taxonomy_species$`matching species IDs`, phylum = taxonomy_species$phylum) %>%
  filter(!is.na(phylum)) %>%
  pivot_longer(cols = -c(species, phylum), names_to = "sample", values_to = "count") %>%
  group_by(phylum) %>%
  mutate(total_reads = sum(count)) %>%
  ungroup() %>%
  mutate(phylum = reorder(phylum, total_reads))

# plot raw bar graph
ggplot(count_species, aes(x = sample, y = count, fill = phylum)) +
  geom_col(position = "stack") +
  labs(x = "Sample", y = "Count", fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = file.path("4.results", "species_raw_bar.png"), plot = last_plot(), bg = "white", dpi = 300)  
```

```{figure} species_raw_bar.png
:name: species raw bar

: An unformatted stacked bar chart showing read count abundance of phyla.
```

Since the number of sequences differ between samples, it might sometimes be easier to use relative abundances to visualise the data.

```{code-block} R
# plot rel abund bar graph
ggplot(count_species, aes(x = sample, y = count, fill = phylum)) +
  geom_col(position = "fill") +
  labs(x = "Sample", y = "Count", fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = file.path("4.results", "species_rel_abund_bar.png"), plot = last_plot(), bg = "white", dpi = 300) 
```

```{figure} species_rel_abund_bar.png
:name: species rel abund bar

: An unformatted relative abundant stacked bar chart showing read count abundance of phyla.
```

While we can start to see that most samples show a similar pattern, except for Heda_1, sometimes it is easier to group the samples to reduce the number of bars in the plot.

```{code-block} R
# combine samples per location
metadata <- metadata %>%
  mutate(sample = rownames(.))
count_location <- count_species %>%
  left_join(metadata %>% select(sample, sampleLocation), by = "sample") %>%
  group_by(sampleLocation, species, phylum) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  group_by(phylum) %>%
  mutate(total_reads = sum(total_count)) %>%
  ungroup() %>%
  mutate(phylum = reorder(phylum, total_reads))

# plot rel abund bar location
ggplot(count_location, aes(x = sampleLocation, y = total_count, fill = phylum)) +
  geom_col(position = "fill") +
  labs(x = "Location", y = "Relative Abundance", fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = file.path("4.results", "location_rel_abund_bar.png"), plot = last_plot(), bg = "white", dpi = 300) 
```

```{figure} location_rel_abund_bar.png
:name: location rel abund bar

: An unformatted relative abundant stacked bar chart showing read count abundance of phyla per location.
```

As you can see, the graphs are easier to interpret when grouped by location, due to the limited number of stacked bars. The high abundance of Chlorophyta, however, is potentially still masking any differences between samples. Furthermore, it might also sstill be difficult to identify which phyla or species are occurring in only one sampling location, due to the number of species detected and the number of colours used. One "trick" we can employ, is to plot the data per phylum or species for all locations. For the tutorial, we will generate the plot per phylum for easy viewing (note: we use the `scales = "free_y"` option to eliminate the issue of the highly abundant Chlorophyta phylum). To generate the species plot, simply change `facet_wrap(~phylum, scales = "free_y")` to `facet_wrap(~species, scales = "free_y")`.

```{code-block} R
# plot per phylum for easy visualisation
ggplot(count_location, aes(fill = sampleLocation, y = total_count, x = sampleLocation)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  ggtitle("Graphing data per phylum") +
  facet_wrap(~phylum, scales = "free_y") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("")
ggsave(filename = file.path("4.results", "per_phylum_bar.png"), plot = last_plot(), bg = "white", dpi = 300) 
```

```{figure} per_phylum_bar.png
:name: per phylum bar

: A formatted bar chart showing phyla abundance per location.
```

Interestingly, when plotting the relative abundance, we will indirectly transform the data to presence-absence, which will aid significantly in distinguishing presence or absence for phyla between sampling locations.

```{code-block} R
# plot rel abund per phylum for easy visualisation
ggplot(count_location, aes(fill = sampleLocation, y = total_count, x = sampleLocation)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  ggtitle("Graphing data per phylum") +
  facet_wrap(~phylum, scales = "free_y") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("")
ggsave(filename = file.path("4.results", "per_phylum_rel_abund_bar.png"), plot = last_plot(), bg = "white", dpi = 300) 
```

```{figure} per_phylum_rel_abund_bar.png
:name: per phylum rel abund bar

: A formatted bar chart showing phyla relative abundance per location.
```

From this plot, we can clearly se that the Chlorophyta, Haptophyta, Chordata, and Porifera phyla were detected at all locations.

````{admonition} Further exploration
:class: tip
What we've shown thus far is, of course, only a very brief overview of some data exploration steps. We suggest you explore this data further at a more interesting taxonomic level compared to phylum, focus on species of interest, or to incorporate a larger proportion of the data by relaxing the BLAST-hit filtering thresholds we've used in this example.
````

## 3. Alpha diversity

In short, alpha diversity is the characterisation of diversity within a site. Researchers use this characterisation to quantify how diverse different systems are. For metabarcoding data, which does not hold true abundance information as the corrleation between eDNA signal strenght and species abundance is weak at best, scientists mainly focus on the simplest approach, i.e., **species richness**. Species richness only considers the presence or absence of a species and disregards differentially weighing abundant and rare species.

Within alpha diversity, metabarcoding studies primarily incorporate **T**axonomic **D**iversity (TD). However, by generating phylogenetic trees, as we've done in the last section, **P**hylogenetic **D**iversity (PD) can also be investigated. Recently, there has been a push to incorporate PD in metabarcoding studies, as intuitively, it makes much more sense that a community comprised of 20 birds within the same family is less diverse than 20 birds from different families. For more information on PD, we recommend these two publications from [Kling et al., 2018](https://royalsocietypublishing.org/doi/epdf/10.1098/rstb.2017.0397) and [Alberdi & Gilbert, 2019](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014).

````{admonition} Hill numbers
:class: tip
We also want to point you to the Hill number framework. Alpha diversity investigations of metabarcoding data can be significantly expanded using this single statistical framework that includes **species richness**, **Shannon index**, and **Simpson index**. This approach, however, requires metabarcoding data to be transformed using a so-called **incidence-based** approach. This approach works by transforming the data set to presence-absence and combining multiple samples from the same site/treatment to obtain the relative abundance information. For more information, we can have a look at the figure below. I also strongly suggest reading the excellent paper by [Alberdi & Gilbert, 2019](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014), which covers all of this in much more detail and in context of a metabarcoding approach. They also provide all the essential citations of the original statistical framework development.

```{figure} hill_numbers_intro.png
:name: hill numbers intro

: Introduction to the Hill number statistical framework. Copyright: [Alberdi & Gilbert, 2019](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014).
```

It should be noted, however, that this data transformation requires an extra level of replication build into the experimental design. For our tutorial data, for example, we would need to transform a single sample to presence absence (0 for absence and 1 for presence). We, then, combine the three samples within a location to create the relative abundance measurement. For example, when an ASV is detected in 3 samples it will receive a value of 3, while a single detection of an ASV across all three samples will receive a value of 1. This approach, however, eliminates the current replication at the location level (we went from 3 replicates to 1 sample per location). If you'd be interested in comparing all locations, it would be statistically impossible due to the lack of replication, i.e., you would be comparing a single sample at each location with each other. A potential solution to this problem for our tutorial data since the locations are grouped into 2 regions, is to statistically compare the regions with each other, as the north region has 6 samples (locations) and south region has 5 samples (locations). If you would like to use this approach while still being able to compare all locations, it would require the collection of multiple grouped samples within a location. Hence, the extra level of replication required for this data transformation.

```{figure} incidence_based_transformation.png
:name: incidence based transformation

: Introduction to incidence-based data transformation. Copyright: [Alberdi & Gilbert, 2019](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014).
```
````

## 4. Beta diversity

Besides the number of species within a site (alpha diversity), species composition or beta diversity is another important aspect of analysis within metabarcoding research. To analyse species composition between sites, we can use ordination and cluster visualisations to depict sample/site/treatment similarity or dissimilarity. PERMANOVA and PERMDISP, on the other hand, are the statistical tests to investigate significant differences between treatments based on species composition between treatments. Lastly, Indicator Species Analysis (ISA) allows us to identify which signals are driving the difference between groups.

Similar to alpha diversity, beta diversity can incorporate **T**axonomic **D**iversity, as well as **P**hylogenetic **D**iversity, something we'll explore during this tutorial.

### 4.1 Visualise ASVs across locations

First, however, it is useful to visualise the overlap in ASV detection between the different locations. Traditionally, Venn diagrams have been used for this purpose. However, drawing Venn diagrams with 11 sections is impractical. Two alternative approaches to Venn diagrams when dealing with a larger number of sets (n > 3) are **UpSet diagrams** and **network graphs**. **UpSet diagrams** are best for showing intersection sizes without drawing overlapping circles and works by using a matrix format with bars to represent intersection frequencies. To draw UpSet diagrams, we can make use of the [UpSetR R package](https://github.com/hms-dbmi/UpSetR). **Network graphs**, on the other hand, are best for highlighting shared elements between groups, whereby nodes are the groups (locations) and edges are the shared elements (ASVs) with thickness showing the overlap size. To draw network graphs, we can make use of the [igraph R package](https://r.igraph.org/) and the [ggraph R package](https://ggraph.data-imaginist.com/).

So, let's generate these plots now to look at the shared and unique ASVs across locations.

Before generating the plots, we need to reformat the data. In this case, we need to combine the metadata file with the count table, as we need to sum the samples by location (which will be our groups). Then, we will transform the data to presence-absence, as these graphs do not include information on ASV abundance or occurrence.

```{code-block} R
## Visualise shared/unique ASVs
# combine metadata and count_table to group per location
location_table <- t(count_table) %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  left_join(metadata %>% select(sample, sampleLocation), by = "sample") %>%
  group_by(sampleLocation) %>%
  summarise(across(-sample, ~ as.integer(sum(.x) > 0))) %>%
  column_to_rownames("sampleLocation") %>%
  t() %>% as.data.frame()
```

Next, we can call the `upset()` function in the [UpSetR R package](https://github.com/hms-dbmi/UpSetR). To keep things tidy, we will only show the 30 most abundant intersections (more on this after we have generated the graph) using the `nintersects = 30` parameter. Please note that we are only showing the default plot and many customisations are possible. For more information, I suggest you read the help documentation (`?upset()`) of the package.

```{code-block} R
# visualise through UpSet
upset(location_table, nsets = 11, nintersects = 30, order.by = "freq", text.scale = 1.2,
      mainbar.y.label = "Number of ASVs", sets.x.label = "Total number of ASVs")
png(filename = file.path("4.results", "upset.png"), width = 10, height = 6, units = "in", res = 300)
upset(location_table, nsets = 11, nintersects = 30, order.by = "freq", text.scale = 1.2,
      mainbar.y.label = "Number of ASVs", sets.x.label = "Total number of ASVs")
dev.off()
```

```{figure} upset.png
:name: UpSet

: UpSet graph visualising shared/unique ASVs across locations.
```

This graph contains several subplots, so let's go over them in detail to help interpret UpSet plots. On the bottom left-hand side, we have the total number of ASVs per location represented as a bar graph. On the top right-hand side, we have the number of ASVs that belong to an intersection represented as a bar graph (only the 30 most abundant intersections are displayed). In the bottom right-hand side, we have the intersections, i.e., to which locations this belongs. For example, in our data set, the most abundant intersection consists of 220 ASVs only belonging to Sani (only coloured dot) and are not found in any other location. 100 ASVs, on the other hand, which is the second most abundant intersection, are found in all 11 locations (all 11 points are coloured).

From this graph, we can clearly see that, while there is a large number of ASVs that are shared between locations, there is also a large number of ASVs that are only present in a single location.

To generate a **network graph**, we have to reformat the data once again. First, we need to calculate the number of shared ASVs between locations and place these values in a matrix.

```{code-block} R
# visualise through network graph
# compute shared ASVs between locations
shared_asvs <- t(location_table) %*% as.matrix(location_table)
```

Next, we need to generate a network graph object from this matrix.

```{code-block} R
# create network graph
net <- graph_from_adjacency_matrix(shared_asvs, mode = "undirected", weighted = TRUE, diag = FALSE)
```

Then, we can generate the plot.

```{code-block} R
# circular layout with optimized spacing
ggraph(net, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(width = weight, alpha = weight), color = "steelblue", strength = 0.2, show.legend = FALSE) +
  geom_node_point(aes(size = 10), fill = "#E1F5FE", shape = 21, color = "white", stroke = 0.5, show.legend = FALSE) +
  geom_node_text(aes(label = name, angle = -((-node_angle(x, y) + 90) %% 180) + 90), size = 3, hjust = "outward", vjust = "outward", repel = FALSE) +
  scale_edge_width(range = c(0.5, 3), name = "Shared ASVs") +
  scale_size_continuous(range = c(4, 12), name = "Degree") +
  theme_void() +
  theme(legend.position = "bottom", plot.margin = unit(c(1,1,1,1), "cm"), panel.background = element_rect(fill = "white", color = NA)) +
  coord_fixed(clip = "off")
ggsave(filename = file.path("4.results", "network_graph.png"), plot = last_plot(), bg = "white", dpi = 300) 
```

```{figure} network_graph.png
:name: network graph

: Network graph visualising shared/unique ASVs across locations.
```

In this circular network graph, the dots represent the 11 sampling locations and the lines represent the connectedness between locations, with width and color linked to the number of ASVs shared.

From this graph, we can see the differing overlap in number of ASVs between locations. However, this plot is also impacted by the total number of ASVs within a location (look at how Akakina, the location with fewest total ASVs, also looks like the least shared ASVs). Finally, this plot also loses the information on unique ASVs, as only shared ASVs are taken into account in this plot.

To include unique ASVs, we can generate a so-called **bipartite network graph**, where we have two types of nodes, one for location and one for ASVs. We can, then, plot the lines to which ASV is linked to which location. For this graph, we need to slightly reformat the data again.

```{code-block} R
# reformat data for Bipartite Network Graph
edge_list <- location_table %>%
  rownames_to_column("ASV") %>%
  pivot_longer(cols = -ASV, names_to = "Location", values_to = "Present") %>%
  filter(Present == 1) %>%
  select(-Present)
```

Next, we can create the network graph object again and plot the graph.

```{code-block} R
# create network
net <- graph_from_data_frame(edge_list, directed = FALSE)

# assign node types and scale
V(net)$type <- ifelse(V(net)$name %in% colnames(location_table), "Location", "ASV")
V(net)$degree <- degree(net)

# visualise using bipartite network graph
ggraph(net, layout = "fr") +
  geom_edge_link(alpha = 0.1, color = "gray80") +
  geom_node_point(aes(size = degree, color = type), show.legend = TRUE) +
  geom_node_text(aes(label = name, color = type), data = . %>% filter(type == "Location"), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Location" = "steelblue", "ASV" = "grey40")) +
  scale_size(range = c(1, 8)) +
  theme_void()
ggsave(filename = file.path("4.results", "bipartite_network_graph.png"), plot = last_plot(), bg = "white", dpi = 300) 
```

```{figure} bipartite_network_graph.png
:name: bipartite network graph

: Bipartite network graph visualising shared/unique ASVs across locations.
```

This plot is a bit difficult to visualise due to the large number of ASVs, but we can clearly see that the Sano location separates in the bipartite network graph, due to the large number of unique ASVs.

````{admonition} Optimise bipartite network graph
:class: tip
To further optimise the bipartite network graph, we could combine all ASVs that share the same location. For example, rather than 220 ASVs unique to Sano, you could have a single node, with a line represented by a weight of 220.
````

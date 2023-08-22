# This file iterates over a master table of fusion transcripts and searches 
# defined regions in a corresponding WGS bam file for gene-spanning discordant read pairs.
# One read should align to the 5' fusion partner and the mate should align to the 3' partner. 

# Input: input/fusion_search_table.txt (a tab delim txt with information on where to search for a particular fusion)

# Output: output/output_spanning_reads.txt w/ standard SAM fields, plus fusion_ID
# query, i.e. whether the read in question came from querying the 5' or 3' fusion partner.
# fiveprime_strand, threeprime strand; in order to define a region to query for softclipped reads



# Load packages and input data

library(tidyverse)
library(optparse)

input_options = list(
  make_option(c("-i", "--input_file"), default = "input/fusion_search_table.txt", help = "Location of input file containing genomic regions"),
  make_option(c("-r", "--ref_genome"), default = "input/ref_genome/hg19.fa", help = "Location of reference genome fasta file"), 
  make_option(c("-d", "--distance_filter"), default = 4000, help = "Minium fragment size of discorand read pairs (default = 4000)")
)

input_file <- parse_args(OptionParser(option_list = input_options))$input_file
ref_genome <- parse_args(OptionParser(option_list = input_options))$ref_genome


sam_fields <- c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE","SEQ","QUAL","OPT")
input <- read_tsv(input_file)

# WGS query and initial filtering
print("Searching WGS data for discordant read pairs...\n")
disc_reads <- input %>% 
  mutate(awk_chr = ifelse(fiveprime_chr == threeprime_chr, "=", threeprime_chr)) %>% 
  mutate(awk_chr2 = ifelse(fiveprime_chr == threeprime_chr, "=", fiveprime_chr)) %>% 
  mutate(cmd = paste0("view -F 512 -F 1024 -F 2048 -q 37 ", path, " -e 'pnext>=", threeprime_search_start, " && pnext<=", threeprime_search_end,"' ", 
    fiveprime_chr, ":",fiveprime_search_start, "-",fiveprime_search_end, " |  awk '($7 == \"", awk_chr, "\")' | awk '!($7 == \"=\" && $9 <= 1000 && $9 >= -1000)'")) %>% 
  mutate(cmd2 = paste0("view -F 512 -F 1024 -F 2048 -q 37 ", path, " -e 'pnext>=", fiveprime_search_start, " && pnext<=", fiveprime_search_end,"' ", 
    threeprime_chr, ":",threeprime_search_start, "-",threeprime_search_end, " |  awk '($7 == \"", awk_chr2, "\")' | awk '!($7 == \"=\" && $9 <= 1000 && $9 >= -1000)'")) %>% 
  mutate(fiveprime_reads = map(cmd, ~system2("samtools", args = .x, stdout = T))) %>% 
  mutate(threeprime_reads = map(cmd2, ~system2("samtools", args = .x, stdout = T)))

# Additional filtering, create output table
print("Applying discorandt read filtering steps...\n")
disc_reads <- disc_reads %>% 
  select(fusion_id, fiveprime_strand, threeprime_strand, fiveprime_reads, threeprime_reads) %>% 
  pivot_longer(cols = c(fiveprime_reads, threeprime_reads), names_to = "query", values_to = "results") %>% 
  unnest(results) %>% 
  separate(results, into = sam_fields, sep = "\t", extra = "merge") %>% 
  distinct() %>%
  # mutate(FLAG = as.numeric(FLAG)) %>% 
  # mutate(to_discard = map_lgl(FLAG, ~any(bitwAnd(.x, 1024), bitwAnd(.x, 512)))) %>% 
  # filter(!as.numeric(MAPQ) < 37) %>% 
  # filter(!to_discard) %>% 
  filter(n() >= 2, .by = c(fusion_id, QNAME))

disc_reads %>% 
  write_tsv("output/discordant_read_pairs.txt")


# Define regions to query soft clipped reads ------------------------------

# We want to find soft-clipped reads close to each discordant read, to see if the soft-clip represents the genomic breakpoint that gives rise to the fusion.
# Discordant reads often form clusters, so instead of extracting soft-clipped reads around every discordant read, we extract them around every *cluster* of discordant reads.
# This minimizes query time and prevents us from fetching the same soft-clipped read more than once. 
print("Searching for breakpoint-supporting reads...\n")
disc_read_clusters <- disc_reads %>% 
  mutate(width = nchar(SEQ), 
    start = as.numeric(POS)) %>% 
  mutate(QNAME_ID = row_number()) %>% 
  select(fusion_id, QNAME_ID, query, width, start, fiveprime_strand, threeprime_strand) %>% 
  mutate(end = start + width) %>% 
  mutate(strand = ifelse(query == "fiveprime_reads", fiveprime_strand, threeprime_strand)) %>% 
  select(-c(fiveprime_strand, threeprime_strand)) %>% 
  mutate(start = ifelse( (query == "fiveprime_reads" & strand == "+") | (query == "threeprime_reads" & strand == "-") , start-50, start-500)) %>% 
  mutate(end = ifelse( (query == "fiveprime_reads" & strand == "+") | (query == "threeprime_reads" & strand == "-") , end+500, end+50)) %>% 
  group_by(fusion_id, query) %>% 
  nest() %>% 
  mutate(iranges = map(data, function(df){
    IRanges::IRanges(df$start, end = df$end) %>% IRanges::reduce() %>% IRanges::as.data.frame() %>% as_tibble()
  })) %>% 
  unnest(iranges) %>% 
  rename(istart = start, iend = end, iwidth = width) %>% 
  unnest(data) %>% 
  ungroup()  %>% 
  mutate(filt = between(start, istart, iend)) %>% 
  filter(filt) %>% 
  select(-c(start, end, width, filt)) %>% 
  rename(start = istart, end = iend, width = iwidth) %>% 
  group_by(fusion_id, query, strand, start, end, width) %>% 
  summarise(n_disc_reads_in_cluster = n(), disc_reads_in_cluster = paste0(QNAME_ID, collapse = ";")) %>% 
  ungroup()

# Assign a cluster ID and associate to linked clusters ----

# Each cluster is given an ID (fusion_id + number). Discordant read pairs link the clusters and tell us where we should align the soft-clipped sequences from each discordant read cluster. 

# Example: A fusion forms four distinct clusters, two on the 5' partner and two on the 3' partner. The 5' clusters are given the IDs 1 and 2, and the 3' given the IDs 3 and 4 respectively. The mates of the reads in cluster 1 all align to cluster 3, and the mates of cluster 2 all align to cluster 4. Soft-clipped reads extracted around cluster 1 are therefore only aligned to the region we extract soft-clipped reads from cluster 3 and vice-versa. Similarly, soft-clips from cluster 2 are aligned to the cluster 4 region and vice-versa. 

# Below code creates a table denoting which clusters are linked together, and what reads beleng to each cluster. 

disc_read_clusters <- disc_read_clusters %>% 
  group_by(fusion_id) %>% 
  nest() %>% 
  mutate(data = map(data, function(df){
    df %>% 
      mutate(cluster_id = 1:nrow(.))
  })) %>% 
  unnest(data) %>% 
  mutate(disc_read_cluster_id = paste0(fusion_id, "_", cluster_id)) %>% 
  select(-cluster_id) %>% 
  ungroup()
  
disc_read_cluster_key <- disc_read_clusters %>% 
  select(fusion_id, query,disc_reads_in_cluster, disc_read_cluster_id) %>% 
  mutate(disc_reads_in_cluster = str_split(disc_reads_in_cluster, ";")) %>% 
  unnest(disc_reads_in_cluster) %>% 
  left_join(
    disc_read_clusters %>% 
      select(fusion_id, query,disc_reads_in_cluster, disc_read_cluster_id) %>% 
      mutate(disc_reads_in_cluster = str_split(disc_reads_in_cluster, ";")) %>% 
      unnest(disc_reads_in_cluster) %>% 
      rename(associated_disc_read_cluster = disc_read_cluster_id,
        associated_query = query), 
    join_by(fusion_id, disc_reads_in_cluster)
    ) %>% 
  filter(query != associated_query) %>% 
  group_by(fusion_id, query, disc_read_cluster_id) %>% 
  summarise(
    disc_reads_in_cluster = paste0(disc_reads_in_cluster, collapse = ";"),
    associated_disc_read_clusters = paste0(unique(associated_disc_read_cluster), collapse = ";")) %>% 
  ungroup()


# Get chromosome information 
disc_read_clusters <- disc_read_clusters %>% 
  left_join(input %>% 
      select(fusion_id, path, fiveprime_chr, threeprime_chr) %>% 
      distinct()) %>% 
  mutate(chr = ifelse(query == "fiveprime_reads", fiveprime_chr, threeprime_chr)) %>% 
  select(-c(fiveprime_chr, threeprime_chr))

# Extract all reads with a high-quality soft clipped end ----
# A softclip supporting a genomic breakpoint can only be on one side of the read, which depends on which fusion partner is being queried and that strand that the gene is on. 

# Start by extracting all reads in the correct location that have a soft clip on the correct side
softclips <- disc_read_clusters %>% 
  mutate(expected_softclip_side = 
      case_when(
        (query == "fiveprime_reads" & strand == "+") | (query == "threeprime_reads" & strand == "-")  ~ "right",
        (query == "fiveprime_reads" & strand == "-") | (query == "threeprime_reads" & strand == "+") ~ "left"
      )) %>% 
  mutate(awk_filter = ifelse(expected_softclip_side == "right", "'$6 ~ /[0-9]+S$/'", "'$6 ~ /^[0-9]+S/'")) %>% 
  mutate(cmd = paste0("view -F 512 -F 1024 -F 2048 -q 1 ", path, " ", chr, ":", start, "-", end, " | awk ", awk_filter)) %>% 
  mutate(softclips = map(cmd, ~system2("samtools", args = .x, stdout = T)))

# Function to calculate the mean sequencing quality of a given string. Converts characters to their ascii values an subtracts 33.
mean_sequencing_quality <- function(string){
  qual = string %>% 
    utf8ToInt() %>% 
    magrittr::subtract(33) %>% 
    mean()
  
  return(qual)
}

# Filter the reads to keep only high quality soft clipped reads (mean seq quality > 15) and a minimum length of 6 nt (the shortest sequence length Novoalign accepts)
softclips <- softclips %>% 
  select(-c(awk_filter, cmd)) %>% 
  unnest(softclips) %>% 
  separate(softclips, into = sam_fields, sep = "\t", extra = "merge") %>% 
  select(-OPT) %>% 
  mutate(softclip_len = str_extract(CIGAR, ifelse(expected_softclip_side == "right", "\\d+S$", "^\\d+S")) %>% str_extract("\\d+") %>% as.numeric()) %>% 
  filter(softclip_len > 5) %>% 
  mutate(softclip_qual_sequence = 
      ifelse(expected_softclip_side == "right", 
        str_sub(QUAL, nchar(QUAL)-softclip_len+1, nchar(QUAL)), 
        str_sub(QUAL, 1, softclip_len))) %>% 
  mutate(softclip_sequencing_quality = map_dbl(softclip_qual_sequence, ~mean_sequencing_quality(.x))) %>% 
  filter(softclip_sequencing_quality > 15) %>% 
  mutate(softclip_sequence = 
      ifelse(expected_softclip_side == "right", 
        str_sub(SEQ, nchar(SEQ)-softclip_len+1, nchar(SEQ)), 
        str_sub(SEQ, 1, softclip_len)))


# Create an BED file containing the coordinates of all the disc read clusters

disc_read_cluster_bed <- disc_read_clusters %>% 
  select(chr, start, end, disc_read_cluster_id) %>% 
  mutate(chr = ifelse(str_detect(chr, "chr"), chr, paste0("chr", chr)))

# Create a fasta file for each disc read cluster
# This fasta file contains the sequences that are to be aligned to this cluster 
# Example: 
# Discordant read cluster 2_1 contains the softclipped sequences GATCGGAAGA and GATCGGAAGAG and 
# is linked to cluster 2_7. This is the only cluster that is linked to cluster 2_7.
# The fasta file for 2_7 would therefore look like this:

# > 2_1_1
# GATCGGAAGA
# > 2_1_2
# GATCGGAAGAG

# The header for each sequence says where the softclip came from. Each softclip that passes filtering is given an ID:

softclips <- softclips %>% 
  mutate(softclip_ID = row_number(), .by = disc_read_cluster_id) %>% 
  mutate(softclip_ID = paste0(disc_read_cluster_id, "_", softclip_ID)) 

dir.create("output/softclip_fasta_files")

softclips %>% 
  select(softclip_ID, softclip_sequence, disc_read_cluster_id) %>% 
  left_join(
    disc_read_cluster_key %>% 
      select(disc_read_cluster_id, associated_disc_read_clusters) %>% 
      mutate(associated_disc_read_clusters = str_split(associated_disc_read_clusters, ";")) %>% 
      unnest(associated_disc_read_clusters),
    multiple = "all"
    ) %>% 
  select(-disc_read_cluster_id) %>% 
  mutate(softclip_ID = paste0(">", softclip_ID)) %>% 
  pivot_longer(cols = c(softclip_ID, softclip_sequence)) %>% 
  select(-name) %>% 
  nest(.by=associated_disc_read_clusters) %>% 
  pwalk(~write_tsv(..2, file = paste0("output/softclip_fasta_files/disc_read_cluster_", ..1, ".fa"), col_names = F))
  
# Write the files to output
disc_read_cluster_key %>% 
  write_tsv("output/discordant_read_cluster_linkage_key.txt")

softclips %>% 
  write_tsv("output/softclipped_reads.txt")

disc_read_clusters %>% 
  write_tsv("output/discordant_read_clusters.txt")

disc_read_cluster_bed %>% 
  write_tsv("output/discordant_read_clusters.bed", col_names = F)


# Use bedtools to create fasta sequences ----

system2(command = "bedtools", args = paste0("getfasta -fi ",ref_genome, " -bed output/discordant_read_clusters.bed -nameOnly > output/discordant_read_cluster_sequences.fa"))

# Split discordant_read_cluster_sequences.fa into multiple smaller fasta files

disc_read_sequences <-  read_lines("output/discordant_read_cluster_sequences.fa") %>% 
  as_tibble()

dir.create("output/discordant_read_cluster_sequences", )

disc_read_sequences %>% 
  mutate(grp = str_extract(value, "(?<=\\>).+")) %>% 
  fill(grp) %>%
  nest(.by = grp) %>% 
  pwalk(~write_tsv(..2, file = paste0("output/discordant_read_cluster_sequences/", ..1, ".fa"), col_names = F ))

# Index these multiple smaller fasta files with novoindex

list.files("output/discordant_read_cluster_sequences/") %>% 
  as_tibble() %>% 
  filter(!str_detect(value, ".ndx$")) %>% 
  mutate(
    cmd = "novoindex", 
    args = paste0("output/discordant_read_cluster_sequences/", value, ".ndx output/discordant_read_cluster_sequences/", value)) %>% 
  pwalk(~system2(command = ..2, args = ..3))
    

# Align soft-clipped sequences ------------------------------------------------



# We now have two files for each discordant read cluster: 
# one in output/discordant_read_cluster_sequences and one in softclip_fasta_files.
# The files in discordant_read_cluster_sequences contain the sequence of the discordant read cluster
# The files in softclip_fasta_files contain the sequences that are to be aligned to that cluster.

# Example: We align the reads in output/softclip_fasta_files/disc_read_cluster_2_1.fa locally to output/discordant_read_cluster_sequences/2_1.fa 
print("Performing alignments...\n")

file.create("output/softclip_alignment.txt")
file.create("output/softclip_alignment_log.txt")


list.files("output/discordant_read_cluster_sequences/") %>% 
  as_tibble() %>% 
  filter(!str_detect(value, ".ndx$")) %>% 
  mutate(cmd = "novoalign",
    args = paste0("-f output/softclip_fasta_files/disc_read_cluster_", value, " -d output/discordant_read_cluster_sequences/", value, ".ndx -r ALL -o SAM -h-1-1 -l6 >> output/softclip_alignment.txt 2>> output/softclip_alignment_log.txt" )) %>% 
  pwalk(~system2(command = ..2, args = ..3))

alignment <- read_tsv("output/softclip_alignment.txt", comment = "@", col_names = sam_fields) %>% 
  filter(MAPQ != 0) %>% 
  select(QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQ)


# Processed aligned soft-clipped reads ------------------------------------

# We now have a list of softclipped reads that successfully aligned in their expected regions. 
# Now we convert these alignment coordinates into full genomic coordinates and identify clusters of softclips.
# Genomic coordinate is calculated as follows:
#  (start coordinate of the search region) + (alignment position of the softclipped read) + (length of aligned sequence, IF the query is 5'/+ strand or 3'/- strand) - (any softclips in the aligned sequence)

alignment <- alignment %>% 
  left_join(
    disc_read_clusters %>% 
      select(disc_read_cluster_id,query, strand, chr, start, end), 
    by = c("RNAME"="disc_read_cluster_id")) %>% 
  mutate(should_shift = (query == "fiveprime_reads" & strand == "+") | (query == "threeprime_reads" & strand == "-")) %>% 
  mutate(
    left_softclip = str_extract(CIGAR, "^\\d+(?=S)") %>% as.numeric() %>% replace_na(0),
    right_softclip = str_extract(CIGAR, "\\d+(?=S)$") %>% as.numeric() %>% replace_na(0),
    insert_sum = str_extract_all(CIGAR, "\\d+(?=I)", simplify = T) %>% as.numeric() %>% map_dbl(~sum(.x) %>% replace_na(0))
    ) %>% 
  mutate(breakpoint_coordinate = start + as.numeric(POS) + nchar(SEQ)*should_shift -  left_softclip*should_shift - right_softclip*should_shift - insert_sum*should_shift) %>% 
  filter(MAPQ >= 20) 


# Use the info of which softclipped sequences aligned to determine the breakpoint in the other fusion partner:

softclips <- softclips %>% 
  mutate(sequence_aligns_to_other_fusion_partner = softclip_ID %in% alignment$QNAME) %>% 
  mutate(breakpoint_coordinate = ifelse(expected_softclip_side == "left", as.numeric(POS) - 1, str_extract(CIGAR, "\\d+(?=M)") %>% as.numeric() %>% magrittr::add(as.numeric(POS))))


# Generate summary tables --------------------------------------------------
# This table includes all fusion IDs used in the input, if they have fusion support, if a breakpoint is identified, number of discordant read pairs supporting the fusion,
# and the coordinates of the most supported breakpoints for the 5' and 3' fusion partners. 

summary_table <- input %>% 
  select(fusion_id) %>% 
  left_join(
    disc_reads %>% 
      group_by(fusion_id) %>% 
      count(fusion_id) %>% 
      mutate(n = n/2) %>% 
      rename(n_discordant_read_pairs = n),
    by = "fusion_id"
  ) %>% 
  mutate(n_discordant_read_pairs = replace_na(n_discordant_read_pairs, 0 )) %>% 
  mutate(discordant_read_support = n_discordant_read_pairs > 0) %>% 
  left_join(
    alignment %>% 
      mutate(fusion_id = str_extract(RNAME, "^\\d+") %>% as.numeric()) %>% 
      distinct(fusion_id) %>% 
      mutate(breakpoint_support = T),
    by = "fusion_id"
  ) %>% 
  mutate(breakpoint_support = replace_na(breakpoint_support, F)) %>% 
  select(fusion_id, discordant_read_support, breakpoint_support)
  
summary_table %>% 
  write_tsv("output/fusion_support_summary_table.txt")

# Create a breakpoint summary table with info from both the softclip file and alignment file

softclip_summary_table <-  alignment %>% 
  mutate(breakpoint_coordinate = paste0(chr, ":", breakpoint_coordinate)) %>% 
  select(softclip_ID = QNAME, query_alignment = query, breakpoint_coordinate_alignment = breakpoint_coordinate) %>% 
  full_join(
    softclips %>% 
      filter(sequence_aligns_to_other_fusion_partner) %>% 
      mutate(breakpoint_coordinate = paste0(chr, ":", breakpoint_coordinate)) %>% 
      select(softclip_ID, query_softclip = query, breakpoint_coordinate_softclip = breakpoint_coordinate)
  ) %>% 
  mutate(
    col1 = paste0(query_alignment, ";", breakpoint_coordinate_alignment),
    col2 = paste0(query_softclip, ";", breakpoint_coordinate_softclip)
    ) %>% 
  select(softclip_ID, col1, col2) %>% 
  pivot_longer(cols = -softclip_ID) %>% 
  select(-name) %>% 
  separate(value, into = c("fusion_partner", "breakpoint_coordinate" ), sep = ";") %>% 
  mutate(fusion_partner = str_replace(fusion_partner, "_reads", "_breakpoint_coordinate")) %>% 
  pivot_wider(names_from = fusion_partner, values_from = breakpoint_coordinate)

softclip_summary_table <- softclip_summary_table %>% 
  left_join(
    softclips %>% 
      select(fusion_id, softclip_ID, softclip_sequence )
    ) %>% 
  select(fusion_id, everything())

softclip_summary_table %>% 
  write_tsv("output/fusion_summary_table_breakpoint_supporting_reads.txt")

print("Done!")
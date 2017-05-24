#!/usr/bin/env python2.7

# M.D. Warman
# May 12, 2017

# This script finds stop codon positions in ordered gff3 files and associates 
# one stop codon position with one gene. It then makes a window around the stop
# codon and outputs one window for each stop codon in a tab-delimited BED6 file:
#
# <1>   Chrom
# <2>   Start 
# <3>   End
# <4>   Name (Gene name)
# <5>   Score (1 is used as a placeholder)
# <6>   Strand
#
# Usage:
# stop_codon_finder.py [genome_annotations.gff3]

import sys
import io
import re

# Checks to see if any files were supplied
if len(sys.argv) == 1:
    print("\nNo files supplied.\n\nUsage: stop_codon_finder.py [FILE]\n")
    quit()

# Sets the window size. Here upstream referes to towards to coding region of
# the gene, while downstream refers to the 3' UTR region:
window_upstream_size = 200
window_downstream_size = 200

# Creates some lists to store important info. The script processes each line
# individually, but doesn't know which CDS will be the end CDS until it
# processes all the CDSs for that gene, so it has to store the information for
# all the CDSs in each gene until it reads through that entire gene, then it
# can pull out the gene name from the list and the stop site, and delete the
# contents of the lists.
current_gene_name_list = list()
current_CDS_start_list = list()
current_CDS_end_list = list()
current_strand_list = list()
current_chr_list = list()
list_pos_counter = 0

# Opens the gff3 file and first separates into a list based on \t
fandle = io.open(sys.argv[1], "rU")
out_fandle = io.open("stop_window_output.bed", "wb")

for line in fandle:
    linestripped = line.strip()
    line_list = linestripped.split("\t")
    if len(line_list) >= 8:                 # bypasses comment lines, etc
        if line_list[2] == "CDS":
            sub_list = line_list[8].split(";") # splits the last field
            current_gene_name_list.append(sub_list[0])
            current_CDS_start_list.append(line_list[3])
            current_CDS_end_list.append(line_list[4])
            current_strand_list.append(line_list[6])
            current_chr_list.append(line_list[0])
            if len(current_gene_name_list) != 1:
                # This if statement checks to see if the current gene has
                # changed since the last time. Once the gene changes, it's time
                # to look at the previous lines and the strand and figure out
                # which CDS contains the stop codon.
                if current_gene_name_list[list_pos_counter] != \
                   current_gene_name_list[list_pos_counter - 1]:
                    # Checks to see if it's the first isoform. For this script,
                    # second isoforms and up will be ignored.
                    if re.search(r"_P01$", \
                     current_gene_name_list[list_pos_counter - 1]):
                        if current_strand_list[list_pos_counter - 1] == "+":
                            current_window_start = (int(current_CDS_end_list[ \
                             list_pos_counter - 1]) - \
                             int(window_upstream_size)) 
                            current_window_end = (int(current_CDS_end_list[ \
                             list_pos_counter - 1]) + \
                             int(window_downstream_size)) 
                            out_fandle.write(current_chr_list[ \
                             list_pos_counter - 1] + "\t")
                            out_fandle.write(str(current_window_start) + "\t")
                            out_fandle.write(str(current_window_end) + "\t")
                            out_fandle.write(current_gene_name_list[ \
                             list_pos_counter - 1] + "\t")
                            out_fandle.write("1" + "\t")
                            out_fandle.write("+" + "\n")
                        if current_strand_list[list_pos_counter - 1] == "-":
                            current_window_start = (int(current_CDS_start_list[ \
                             0]) - int(window_downstream_size)) 
                            current_window_end = (int(current_CDS_start_list[ \
                             0]) + int(window_upstream_size)) 
                            out_fandle.write(current_chr_list[ \
                             list_pos_counter - 1] + "\t")
                            out_fandle.write(str(current_window_start) + "\t")
                            out_fandle.write(str(current_window_end) + "\t")
                            out_fandle.write(current_gene_name_list[ \
                             list_pos_counter - 1] + "\t")
                            out_fandle.write("1" + "\t")
                            out_fandle.write("-" + "\n")
                    # This stores the "current" gene, since the list is looking
                    # back one gene. If it's not stored then added to the blank
                    # list, then it will delete everything and skip genes that
                    # have only 1 CDS. Does the same thing for CDS start and
                    # end and strand.
                    gene_name_to_add = \
                     current_gene_name_list[len(current_gene_name_list) - 1]
                    CDS_start_to_add = \
                     current_CDS_start_list[len(current_CDS_start_list) - 1]
                    CDS_end_to_add = \
                     current_CDS_end_list[len(current_CDS_end_list) - 1]
                    strand_to_add = \
                     current_strand_list[len(current_strand_list) - 1]
                    chr_to_add = \
                     current_chr_list[len(current_chr_list) - 1]
                    current_gene_name_list = list()
                    current_gene_name_list.append(gene_name_to_add)
                    current_CDS_start_list = list()
                    current_CDS_start_list.append(CDS_start_to_add)
                    current_CDS_end_list = list()
                    current_CDS_end_list.append(CDS_end_to_add)
                    current_strand_list = list()
                    current_strand_list.append(strand_to_add)
                    current_chr_list = list()
                    current_chr_list.append(chr_to_add)
                    list_pos_counter = 0
            list_pos_counter = list_pos_counter + 1

fandle.close()
out_fandle.close()

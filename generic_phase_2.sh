#' ---
#' title: "Generic - phase 2"
#' author: "Sami Kilpinen"
#' date: "April 21st, 2020"
#' ---


# In Bash in Linux
cd splitted
rm -R count_reads_summit_peaks_output/
rm -R macs2_results/
rm barcodes*
rm clade*.bam
tail -n+2 clades.csv | grep -Po '\d+' | sort | uniq | while read i;
    do sed -En 's/^([ATCG]+).*'"$i"'$/\1/p' clades.csv | awk '{print "TAG_CB_"$1"-1.bam.filt.sorted.bam"}' > barcodes$i;
    done;
    
find ./ -ipath "./barcodes*" | parallel -P 8 --progress --eta "bamtools merge -list {} -out 'clade_id_'{/.}'.bam'"

# Run MACS2
find ./ -ipath "./clade_id_*.bam" | parallel -P 8 --progress --eta 'macs2 callpeak --nomodel --format BAMPE --name E14_vM --gsize 1.87e9 --keep-dup all --extsize 200 --shift -100 -t {} --outdir ./macs2_results/{}'

# Append
find ./macs2_results/ -iname '*narrowPeak' -exec cat {} >> combined.summit.beds \;
# Sort
bedtools sort -i combined.summit.beds > combined.summit.sorted.bed
# Merge
bedtools merge -i combined.summit.sorted.bed > merged.summits.bed

# PAck for further analysis
#tar -zcvf peaks.for.clades.tar.gz macs2_results/

#' Remove blacklisted regions from peak detection
bedtools subtract -a merged.summits.bed -b ../mm10-blacklist.v2.bed > merged.summits.filt.bed

# Parallel run of bedtools to check coverage of reads of each cell against the summit bins
mkdir count_reads_summit_peaks_output
find ./ -iname '*.sorted.bam' | parallel -P 8 --progress --eta --plus 'bedtools coverage -a ./merged.summits.filt.bed -b {} > ./count_reads_summit_peaks_output/{/.}.summit.peaks.txt'

# Check coverage of each unique cell in regions defined by the bed file in marker genes bed
#find ./ -iname '*.sorted.bam' | parallel -P 10 --progress --eta --plus 'bedtools coverage -a ~/data_scr/mm10/marker.coords.bed -b {} > ./count_reads_marker_genes_output/{/.}.peaks.txt'

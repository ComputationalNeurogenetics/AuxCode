# Preparation for scATAC-seq

# Adaptation of Cusanovich et al. method - the TF-IDF SVD double layer peak calling for scATACseq - to use Snap-data objects


## Snap atac tools part


```bash
# Copy bam file in CSC to small-cluster with NVME drive (and 10 CPUs with 256GB memory)
cp ~/data_scr/E14_vMB_1_ATAC_v1-2/orig/outs/e14vmb1.possorted_bam.bam ./
cp ~/data_scr/E14_vMB_1_ATAC_v1-2/orig/outs/singlecell.csv ./e14vmb1.singlecell.csv

cp ~/data_scr/E14_vR1_2_ATAC_v1-2/orig/outs/possorted_bam.bam ./e14vr1.possorted_bam.bam &
cp ~/data_scr/E14_vR1_2_ATAC_v1-2/orig/outs/singlecell.csv ./e14vr1.singlecell.csv
cp ~/data_scr/E14_vR1_2_ATAC_v1-2/orig/outs/filtered_peak_bc_matrix.h5 ./e14vr1.filtered_peak_bc_matrix.h5
cp ~/data_scr/E14_vR1_2_ATAC_v1-2/orig/outs/fragments.tsv.gz ./e14vr1.fragments.tsv.gz
cp ~/data_scr/E14_vR1_2_ATAC_v1-2/orig/outs/fragments.tsv.gz.tbi ./e14vr1.fragments.tsv.gz.tbi

cp /scratch/project_2001539/kaia/E14_DI_1_ATAC_v1-2_new/outs/possorted_bam.bam ./e14di.possorted_bam.bam
cp /scratch/project_2001539/kaia/E14_DI_1_ATAC_v1-2_new/outs/singlecell.csv ./e14di.singlecell.csv

cp ~/data_scr/E15_vMB_vR1_2_ATAC_v1-2/orig/outs/possorted_bam.bam e15vmbvr1.possorted_bam.bam &
cp ~/data_scr/E15_vMB_vR1_2_ATAC_v1-2/orig/outs/singlecell.csv ./e15vmbvr1.singlecell.csv
cp ~/data_scr/E15_vMB_vR1_2_ATAC_v1-2/orig/outs/filtered_peak_bc_matrix.h5 ./e15vmbvr1.filtered_peak_bc_matrix.h5
cp ~/data_scr/E15_vMB_vR1_2_ATAC_v1-2/orig/outs/fragments.tsv.gz ./e15vmbvr1.fragments.tsv.gz
cp ~/data_scr/E15_vMB_vR1_2_ATAC_v1-2/orig/outs/fragments.tsv.gz.tbi ./e15vmbvr1.fragments.tsv.gz.tbi


# Copy mm10 genome file
cp ~/data_scr/mm10/mm10.genome.filt ./
```


```bash
# Activate conda environment with needed tools
conda activate /projappl/project_2001539/bioconda3_env/sami_env
```


```bash
# extract the header file
samtools view ./e14vmb1.possorted_bam.bam -H > ./e14vmb1.possorted_bam.bam.header.sam
samtools view ./e14vr1.possorted_bam.bam -H > ./e14vr1.possorted_bam.bam.header.sam
samtools view ./e14di.possorted_bam.bam -H > ./e14di.possorted_bam.bam.header.sam
samtools view ./e15vmbvr1.possorted_bam.bam -H > ./e15vmbvr1.possorted_bam.bam.header.sam

# Modify BAM so that CB tags are embedded in the read name
cat <( cat e14vmb1.possorted_bam.bam.header.sam) <( samtools view ./e14vmb1.possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) | samtools view -bS - > ./e14vmb1.possorted_bam.snap.bam &
cat <( cat e14vr1.possorted_bam.bam.header.sam) <( samtools view ./e14vr1.possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) | samtools view -bS - > ./e14vr1.possorted_bam.snap.bam &
cat <( cat e14di.possorted_bam.bam.header.sam) <( samtools view ./e14di.possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) | samtools view -bS - > ./e14di.possorted_bam.snap.bam &
cat <( cat e15vmbvr1.possorted_bam.bam.header.sam) <( samtools view ./e15vmbvr1.possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) | samtools view -bS - > ./e15vmbvr1.possorted_bam.snap.bam &
```


```bash
# Check
samtools view e14vmb1.possorted_bam.snap.bam | cut -f 1 | head
```


```bash
# Sort the bam file
samtools sort -n -@ 10 -m 12G e14vmb1.possorted_bam.snap.bam > e14vmb1.possorted_bam.snap.sorted.bam &
samtools sort -n -@ 4 -m 8G e14vr1.possorted_bam.snap.bam > e14vr1.possorted_bam.snap.sorted.bam &
samtools sort -n -@ 2 -m 8G e14di.possorted_bam.snap.bam > e14di.possorted_bam.snap.sorted.bam &
samtools sort -n -@ 4 -m 8G e15vmbvr1.possorted_bam.snap.bam > e15vmbvr1.possorted_bam.snap.sorted.bam &
```


```bash
# Generate snap file, with filtering TODO: (write open)
snaptools snap-pre  \
--input-file=e14vmb1.possorted_bam.snap.sorted.bam  \
--output-snap=e14vmb1.possorted_bam.snap \
--genome-name=mm10  \
--genome-size=mm10.genome.filt  \
--min-mapq=30  \
--min-flen=0  \
--max-flen=1000  \
--keep-chrm=TRUE  \
--keep-single=TRUE  \
--keep-secondary=False  \
--overwrite=True  \
--min-cov=100  \
--verbose=True &

snaptools snap-pre  \
--input-file=e14vr1.possorted_bam.snap.sorted.bam  \
--output-snap=e14vr1.possorted_bam.snap \
--genome-name=mm10  \
--genome-size=mm10.genome.filt  \
--min-mapq=30  \
--min-flen=0  \
--max-flen=1000  \
--keep-chrm=TRUE  \
--keep-single=TRUE  \
--keep-secondary=False  \
--overwrite=True  \
--min-cov=100  \
--verbose=True &

snaptools snap-pre  \
--input-file=e14di.possorted_bam.snap.sorted.bam  \
--output-snap=e14di.possorted_bam.snap \
--genome-name=mm10  \
--genome-size=mm10.genome.filt  \
--min-mapq=30  \
--min-flen=0  \
--max-flen=1000  \
--keep-chrm=TRUE  \
--keep-single=TRUE  \
--keep-secondary=False  \
--overwrite=True  \
--min-cov=100  \
--verbose=True &

snaptools snap-pre  \
--input-file=e15vmbvr1.possorted_bam.snap.sorted.bam  \
--output-snap=e15vmbvr1.possorted_bam.snap \
--genome-name=mm10  \
--genome-size=mm10.genome.filt  \
--min-mapq=30  \
--min-flen=0  \
--max-flen=1000  \
--keep-chrm=TRUE  \
--keep-single=TRUE  \
--keep-secondary=False  \
--overwrite=True  \
--min-cov=100  \
--verbose=True &

```


```bash
# Binary mat to snap object, with window sizes of 1,5 and 10kbp
snaptools snap-add-bmat \
    --snap-file=e14vmb1.possorted_bam.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True &

snaptools snap-add-bmat \
    --snap-file=e14vr1.possorted_bam.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True &

snaptools snap-add-bmat \
    --snap-file=e14di.possorted_bam.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True &

snaptools snap-add-bmat \
    --snap-file=e15vmbvr1.possorted_bam.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True &

```

## Custom tools part


```bash
# In CSC
cat ./orig/outs/singlecell.csv | sed -rn 's/^([ACGT]+-1).*(cell_[0m-9]+).*/\1,\2/p' >> uniq_cell_ids.csv
```

Filter BAM file with FilterBAMByCB.py (see http://localhost:8888/notebooks/Research/Jupyter%20Notebooks/FilterBAMByCB.ipynb) to reduce the number of resulting BAM files.

Split singular BAM files to a BAM file for each unique barcode that has been filtered from 10x Genomics output file singlecell.csv. Note, relies on 10x Genomics to identify which barcodes have enough reads to be categorized as cell. BAM files are pre-filtered to remove reads that do not map to these cell_ids as otherwise following bamtools operation would generate tens-hundreds of thousands extra files.

Please note that due to the ulimit -n settings and restrictions following operation is very hard or impossible to run on like Mac laptop or similar capacity machine. Hard Mac Os ulimit -n is quite close to number of resulting files.



```bash
bamtools split -in ./possorted_bam_CB_filt.bam -tag CB
```


```bash
# Filter reads per each BAM so that we accept only mapped, paired and no duplicates, with mapping quality of at least 30 sort as well
find ./*-1.bam | parallel -P 6 --progress --eta --plus 'samtools view -F 1540 -f 3 -q 30 -u {} | samtools sort > {/^.*filt\./}.filt.sorted.bam'
```

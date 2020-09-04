import pysam
import csv
import os
import sys

class Printer():
    """Print things to stdout on one line dynamically"""
    def __init__(self,data):
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()

barcode_dict = {}
with open('./uniq_cell_ids.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        barcode_dict[row[0]] = row[1]

barcodes = set(x for x in barcode_dict.keys())

fin = pysam.AlignmentFile("./orig/outs/possorted_bam.bam", "rb",threads=2)
out = pysam.AlignmentFile("./possorted_bam_CB_filt.bam", "wb", template = fin)


n = 0
for read in fin:
    n=n+1
    Printer(n)
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    else: 
        continue
    cell_id = barcode_dict.get(cell_barcode)
    if cell_id:
        out.write(read)


fin.close()
out.close()
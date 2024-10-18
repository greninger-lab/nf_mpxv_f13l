#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import math
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_ref_record(ref_path):
    ref_records = list(SeqIO.parse(open(ref_path, "r"), "fasta"))
    
    if len(ref_records) == 0:
        logger.error(f"The given input file {args.ref} contains no valid fasta records!")
        sys.exit(2)
    elif len(ref_records) > 1:
        logger.error(f"The given input file {args.ref} contains more than 1 valid fasta record!")
        sys.exit(2)
    else:
        return ref_records[0]
        

def get_aa_seq(start,end,orientation,record):
    if orientation == '-':
        return record.seq[start:end].reverse_complement().translate(to_stop=True)
    else:
        return record.seq[start:end].translate(to_stop=True)
        
        
def parse_gff(gff_path, ref_path):
    cds_dict = {}
    
    ref_record = get_ref_record(ref_path)

    with open(gff_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            orientation = fields[6]

            attributes = fields[8].split(';')
            attr_dict = {}
            for attr in attributes:
                attr = attr.strip()
                if attr:
                    key, value = attr.split('=')
                    attr_dict[key] = value

            if feature_type == 'CDS':
                product = attr_dict.get('product', '')
                cds_dict[(start, end)] = {  'name': product,
                                            'orientation': orientation,
                                            'nt': str(ref_record.seq[start-1:end]),
                                            'aa': str(get_aa_seq(start-1,end,orientation,ref_record))
                                         }

    return cds_dict

def find_matching_cds_range(locus, cds_dict):
    matching_cds_ranges = []
    for key in cds_dict.keys():
        if isinstance(key, tuple) and len(key) == 2:
            if locus >= key[0] and locus <= key[1]:
                matching_cds_ranges.append(key)

    if len(matching_cds_ranges) > 0:
        #return the first matching cds range, there should only be one except for some regions in M2
        return matching_cds_ranges[0]
    else:
        return np.nan

def find_gene_product(locus, cds_dict):
    cds_range=find_matching_cds_range(locus, cds_dict)
    if isinstance(cds_range, tuple):
        return cds_dict[cds_range]['name']
    else:
        if math.isnan(cds_range):
            return None
    
    
def get_alt_aa(alt, alt_aa, gene_product):
    if alt[0] == '-' and gene_product:
        return '-'
    elif alt[0] == '+' and gene_product:
        return '+'
    else:
        return alt_aa
        

def get_snpid(alt, position, ref):
    # SNPID special format for indels e.g. TAT39410del - TAT is deleted at position 39410
    if alt[0] == "-":
        return alt[1:] + position + "del"
    elif alt[0] == "+":
        return alt[1:] + position + "ins"
    # else SNPID reference nt position query nt e.g. T39411C - T at position 39411 is a C in query nt
    else:
        return ref + position + alt
    

def get_pos_ref_aa(position, alt_aa, pos_aa, ref_aa, alt, cds_dict):
    if alt_aa == '-' or alt_aa == '+':
        cds_range = find_matching_cds_range(position, cds_dict)
        if isinstance(cds_range, tuple):
            orientation = cds_dict[cds_range]['orientation']
            aa_seq = cds_dict[cds_range]['aa'] 
            # Check orientation and add deletion length to nt coordinate if CDS is in reverese orientation.
            # This should provide first affected aa.
            if orientation == '-' and alt_aa == '-':
                nt_length = cds_range[1] - (position + len(alt) - 1)
            # If orientation is reverse and there is an insertion, position + 1 (other side of insert)
            elif orientation == '-' and alt_aa == '+':
                nt_length = cds_range[1] - position          
            else:
                nt_length = position - cds_range[0] + 1
            
            # Get the first nt of the indel by increasing sequence length by 1
            nt_length = nt_length + 1
            
            if nt_length > 0:
                # Get length in aa space
                remainder = nt_length % 3
                if remainder > 0:
                    nt_length = nt_length - remainder
                    aa_pos = int((nt_length/3) + 1)
                else:
                    aa_pos = int(nt_length/3)
                
                if aa_pos == len(aa_seq) + 1:
                    aa = '*'
                elif aa_pos > len(aa_seq) + 1:
                    aa = '*<'
                else:
                    aa = aa_seq[aa_pos-1]
            else:
                aa_pos = -1
                aa = '?'
                
            #print(f"{position} {aa} {nt_length} {len(aa_seq)} {len(alt)} {remainder} {aa_pos} {orientation} {alt_aa} {cds_range}")
            return (aa_pos, aa)
    else:
        return (pos_aa, ref_aa)    


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Edit iVar variants tsv output to specific format",
        epilog="Example: python edit_ivar_variants.py ivar_variants.tsv ref.gff ref.fasta ivar_variants.edited"
    )
    parser.add_argument(
        "ivar_variants",
        metavar="IVAR_VARIANTS",
        help="iVar variants tsv",
        type=Path
    )
    parser.add_argument(
        "gff",
        metavar="GFF",
        help="Reference gff file",
        type=Path
    )
    parser.add_argument(
        "ref",
        metavar="REF",
        help="Reference fasta file",
        type=Path
    )    
    parser.add_argument(
        "file_out_name",
        metavar="FILE_OUT_NAME",
        help="Output file name",
        type=str
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )   

    return parser.parse_args(argv)

def main(argv=None):
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.ivar_variants.is_file():
        logger.error(f"The given input file {args.ivar_variants} was not found!")
        sys.exit(2)
    if not args.gff.is_file():
        logger.error(f"The given input file {args.gff} was not found!")
        sys.exit(2)
    if not args.ref.is_file():
        logger.error(f"The given input file {args.ref} was not found!")
        sys.exit(2)

    sample_name=args.ivar_variants.name.split('.')[0]

    cds_dict=parse_gff(args.gff, args.ref)
    df_variants=pd.read_csv(args.ivar_variants,sep='\t')
    if not df_variants.empty:
        df_variants['GENE_PRODUCT']=df_variants['POS'].apply(lambda x: find_gene_product(x,cds_dict))
        df_variants['POS_AA']=df_variants['POS_AA'].astype('Int64')
        df_variants['ALT_AA']=df_variants.apply(lambda x: get_alt_aa(x.ALT,x.ALT_AA,x.GENE_PRODUCT), axis=1)
        df_variants['SNPID']=df_variants.apply(lambda x: get_snpid(x.ALT,str(x.POS),x.REF), axis=1)
        
        df_variants['POS_AA']=df_variants.apply(lambda x: get_pos_ref_aa(x.POS,x.ALT_AA,x.POS_AA,x.REF_AA,x.ALT,cds_dict)[0], axis=1)
        df_variants['REF_AA']=df_variants.apply(lambda x: get_pos_ref_aa(x.POS,x.ALT_AA,x.POS_AA,x.REF_AA,x.ALT,cds_dict)[1], axis=1)

        df_variants_new=df_variants[['GENE_PRODUCT','POS_AA','REF_AA','ALT_AA','TOTAL_DP', 'ALT_DP','ALT_FREQ','POS','SNPID']]
        df_variants_new.insert(0,'SAMPLE',sample_name)
        df_variants_new.to_csv(f"{args.file_out_name}.tsv",sep='\t',index=False)
    else:
        with open(f"{args.file_out_name}.tsv", "w") as no_data_out:
            no_data_out.write('SAMPLE\tGENE_PRODUCT\tPOS_AA\tREF_AA\tALT_AA\tTOTAL_DP\tALT_DP\tALT_FREQ\tPOS\tSNPID\n')
            
    #df_variants_nonsynonamous= \
    #    df_variants_new[(df_variants_new['REF_AA'] != df_variants_new['ALT_AA']) \
    #    & df_variants_new['REF_AA'].notnull() \
    #    & df_variants_new['ALT_AA'].notnull()]
    #df_variants_nonsynonamous.to_csv(f"{sample_name}.variants.nonsynonamous.tsv",sep='\t',index=False)

if __name__ == "__main__":
    sys.exit(main())


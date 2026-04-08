import gzip
import pandas as pd
from collections import defaultdict
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import warnings
from Bio import BiopythonWarning
import os

warnings.simplefilter('ignore', BiopythonWarning)

# --- CONFIGURATION ---
INPUT_VCF = "input_variants.vcf"
REF_GTF = "Homo_sapiens.GRCh38.105.gtf.gz"
REF_FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
# ---------------------

def get_net_charge(seq):
    if not seq or len(seq) == 0: return 0.0
    clean_seq = seq.replace('X','A').replace('B','D').replace('Z','E').replace('*','')
    if len(clean_seq) == 0: return 0.0
    try: return ProteinAnalysis(clean_seq).charge_at_pH(7.4)
    except: return 0.0

def build_tx_map(gtf_path):
    tx_exons, tx_cds, tx_strand, tx_gene = defaultdict(list), defaultdict(list), {}, {}
    open_func = gzip.open if gtf_path.endswith('.gz') else open
    
    with open_func(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            feat = parts[2]
            if feat not in ['exon', 'CDS']: continue
            
            chrom, strand = parts[0].replace('chr', ''), parts[6]
            tx_id, gene_name, biotype = None, None, None
            for attr in parts[8].split(';'):
                if 'transcript_id' in attr: tx_id = attr.split('"')[1].split('.')[0]
                elif 'gene_name' in attr: gene_name = attr.split('"')[1]
                elif 'gene_biotype' in attr: biotype = attr.split('"')[1]
                    
            if biotype != 'protein_coding' or not tx_id or not gene_name: continue
            tx_strand[tx_id], tx_gene[tx_id] = (chrom, strand), gene_name
            start, end = int(parts[3]), int(parts[4])
            if feat == 'exon': tx_exons[tx_id].append((start, end))
            elif feat == 'CDS': tx_cds[tx_id].append((start, end))

    canonical_txs, tx_final_junction_pos = set(), {}
    for tx, exons in tx_exons.items():
        canonical_txs.add(tx)
        strand = tx_strand[tx][1]
        if len(exons) == 1: tx_final_junction_pos[tx] = 0 
        else:
            sorted_exons = sorted(exons, key=lambda x: x[0], reverse=(strand=='-'))
            tx_final_junction_pos[tx] = sum((abs(e[1] - e[0]) + 1) for e in sorted_exons[:-1])
            
    return tx_exons, tx_cds, tx_strand, canonical_txs, tx_gene

def extract_genomic_seq(fasta, chrom, intervals, strand):
    chrom_query = chrom if chrom in fasta else f"chr{chrom}"
    if chrom_query not in fasta: return None
    seq = "".join(str(fasta[chrom_query][start-1:end]) for start, end in sorted(intervals, key=lambda x: x[0]))
    return str(Seq(seq).reverse_complement()) if strand == '-' else seq

def parse_variants():
    if not os.path.exists(INPUT_VCF): 
        print(f"Error: {INPUT_VCF} not found.")
        return
        
    tx_exons, tx_cds, tx_strand, canonical_txs, tx_gene = build_tx_map(REF_GTF)
    genome = Fasta(REF_FASTA)
    results = []
    
    with open(INPUT_VCF, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 8: continue
            
            chrom, pos, ref, alt, info = parts[0].replace('chr', ''), int(parts[1]), parts[3], parts[4], parts[7]
            if ',' in alt: continue 
            
            ann_blocks = [x[4:] for x in info.split(';') if x.startswith('ANN=')]
            if not ann_blocks: continue
            
            target_tx = None
            for block in ann_blocks[0].split(','):
                b_parts = block.split('|')
                if len(b_parts) > 6 and 'frameshift_variant' in b_parts[1]:
                    target_tx = b_parts[6].split('.')[0]
                    break
            
            if not target_tx or target_tx not in tx_exons or target_tx not in tx_cds: continue
            gene, strand = tx_gene[target_tx], tx_strand[target_tx][1]
            
            wt_pos_list = []
            if strand == '+':
                for start, end in sorted(tx_exons[target_tx], key=lambda x: x[0]): wt_pos_list.extend(range(start, end + 1))
                start_genomic = min([x[0] for x in tx_cds[target_tx]])
            else:
                for start, end in sorted(tx_exons[target_tx], key=lambda x: x[0], reverse=True): wt_pos_list.extend(range(end, start - 1, -1))
                start_genomic = max([x[1] for x in tx_cds[target_tx]])
                
            wt_seq = extract_genomic_seq(genome, chrom, tx_exons[target_tx], strand)
            if not wt_seq or len(wt_seq) != len(wt_pos_list): continue
            
            genomic_ref_positions = list(range(pos, pos + len(ref)))
            indices = [wt_pos_list.index(p) for p in genomic_ref_positions if p in wt_pos_list]
            if not indices: continue
            
            min_idx, max_idx = min(indices), max(indices)
            alt_seq = str(Seq(alt).reverse_complement()) if strand == '-' else alt
            if (len(alt_seq) - (max_idx - min_idx + 1)) % 3 == 0: continue
                
            mut_seq = wt_seq[:min_idx] + alt_seq + wt_seq[max_idx+1:]
            
            try:
                wt_start_idx = wt_pos_list.index(start_genomic)
                mut_start_idx = wt_pos_list[:min_idx].index(start_genomic) if start_genomic in wt_pos_list[:min_idx] else min_idx
            except ValueError: continue
            
            wt_protein = str(Seq(wt_seq[wt_start_idx:wt_start_idx + len(wt_seq[wt_start_idx:]) - (len(wt_seq[wt_start_idx:]) % 3)]).translate(to_stop=True))
            
            mut_orf = mut_seq[mut_start_idx:]
            mut_protein_full = ""
            for i in range(0, len(mut_orf)-2, 3):
                aa = str(Seq(mut_orf[i:i+3]).translate())
                if aa == '*': break
                mut_protein_full += aa
                
            diverge_idx = 0
            for i in range(min(len(wt_protein), len(mut_protein_full))):
                if wt_protein[i] != mut_protein_full[i]:
                    diverge_idx = i; break
            else: diverge_idx = min(len(wt_protein), len(mut_protein_full))
                
            mut_tail = mut_protein_full[diverge_idx:]
            wt_tail = wt_protein[diverge_idx:diverge_idx+len(mut_tail)] if len(wt_protein) > diverge_idx else ""
            
            local_start = mut_start_idx + (diverge_idx * 3)
            local_end = mut_start_idx + (len(mut_protein_full) * 3)
            local_window = mut_seq[local_start:local_end]
            
            if len(local_window) > 0 and len(mut_tail) > 0:
                local_gc = ((local_window.count('G') + local_window.count('C') + local_window.count('g') + local_window.count('c')) / len(local_window)) * 100
                charge_diff = get_net_charge(mut_tail) - get_net_charge(wt_tail)
                
                results.append({
                    'Variant_ID': f"{gene}_{pos}",
                    'Gene': gene,
                    'Transcript_ID': target_tx,
                    'WT_Length_AA': len(wt_protein),
                    'De_Novo_Length_AA': len(mut_tail),
                    'Tail_GC_Percent': round(local_gc, 2),
                    'Charge_Difference': round(charge_diff, 2)
                })

    df = pd.DataFrame(results)
    df.to_csv("variant_seq_metrics.csv", index=False)
    print(f"Sequence parsing complete. Outputs saved to variant_seq_metrics.csv")

if __name__ == "__main__":
    parse_variants()
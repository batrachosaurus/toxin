import os
import sys
import shutil
import json
import zlib
import time
import traceback
import urllib.request

import numpy as np
import pandas as pd
from Bio import Entrez

from eaglib._utils.logging import eagle_logger
from eaglib.alignment import SeqsProfileInfo, SeqsProfile
from eaglib.seqs import SeqsDict, load_fasta_to_dict

np.random.seed(42)

ass_sum_path = "../input/bacteria_assembly_summary.txt"

def get_taxonomy(tax_id):
    tax_keys = ["superkingdom", "phylum", "clade", "class", "order", "family", "genus", "species"]
    tax_dict = {tax_key: None for tax_key in tax_keys}
    
    for _ in range(5):
        try:
            record = Entrez.efetch(db="taxonomy", id=tax_id, retmode='xml')
            tax_info = Entrez.read(record)[0]
        except:
            time.sleep(5)
            continue
        else:
            tax_dict["species"] = tax_info['ScientificName']
            for lin_tax in tax_info['LineageEx']:
                if lin_tax['Rank'] in tax_dict:
                    tax_dict[lin_tax['Rank']] = lin_tax['ScientificName']
                    
            return [tax_dict[tax_key] for tax_key in tax_keys]
    
    eagle_logger.info("\n\nfailed to retrieve taxid: " + str(tax_id) + "\n" + str(rownum) + " " + name)
    eagle_logger.info(traceback.format_exc())
    sys.exit(1)

toxn_profile = SeqsProfile(SeqsProfileInfo.load_from_dict({
    'name': 'ToxN',
    'path': '../input/PF13958.hmm',
    'type': 'protein',
    'weight': 1.0,
    'method': 'hmmer'
}))


Entrez.email = "selifonov2002@gmail.com"
db_dir = "bacteria"
tcds_path = os.path.join(db_dir, "translated_cds.faa")


processed_ac = list()
bact_df = pd.read_csv(ass_sum_path, sep="\t", low_memory=False).query("assembly_level=='Complete Genome'")
#bact_df_raref = bact_df.sample(frac=1).groupby("organism_name").head(3).reset_index(drop=True)
eagle_logger.info(" %s genomes to prepare" % len(bact_df))
for rownum, row in bact_df.iterrows():
    ac = row['assembly_accession']  # id field in genomes_table
    asm = row['asm_name']
    name = row['organism_name'] + ("" if pd.isna(row['infraspecific_name']) else " " + row['infraspecific_name'])
    taxonomy = get_taxonomy(row['species_taxid'])
    ftp_prefix = (row['ftp_path'] + "/" + ac + "_" + asm).replace(" ", "_")
    fna_seq = [ftp_prefix+"_genomic.fna.gz"]
    btc_seqs = [None]
    btc_seqs_dict = dict()
    target_name, target_len, ali_from, ali_to, seq_score = None, -1, -1, -1, -1
    status = 0  # 0 for no hits; 1 for 1/more hits; 2 for error
    tback = None
    
    try:
        with open(tcds_path, 'wb') as tcds_f:
            tcds_f.write(zlib.decompress(urllib.request.urlopen(ftp_prefix+"_translated_cds.faa.gz").read(), 15+32))
        psr_toxn_df = toxn_profile.search(seqdb=tcds_path, threads=3)
        if not psr_toxn_df.empty:
            status = 1
            max_score_toxn = psr_toxn_df.loc[psr_toxn_df["domain_score"].idxmax()]
            target_name, target_len = max_score_toxn["target name"], max_score_toxn["tlen"]
            ali_from, ali_to = max_score_toxn["ali_from"], max_score_toxn["ali_to"]
            seq_score = max_score_toxn["seq_score"]

            btc_seqs_dict[toxn_profile.name] = load_fasta_to_dict(fasta_path=tcds_path)[target_name][ali_from-1: ali_to]
            btc_seqs = [os.path.join(db_dir, ac+"_btc.fasta")]
            SeqsDict.load_from_dict(btc_seqs_dict).dump(btc_seqs[0])
    except:
        status = 2
        tback = traceback.format_exc()
        print(tback)

    hit_data = {
        "row": int(rownum), "id": ac, "name": name, "taxonomy": taxonomy,
        "btc_seqs": btc_seqs,
        "target_name": target_name, "target_len": int(target_len), "ali_coords": [int(ali_from), int(ali_to)],
        "seq_score": float(seq_score),
        "fna_seq": fna_seq,
        "status": status, "traceback": tback
    }
    processed_ac.append(hit_data)
    eagle_logger.info(" " + ac + "\t" + str(status) + "\t" + name)
    # if _ >= 10: break ###

with open("processed_bacteria.json", "w") as processed_ac_f:
    json.dump(processed_ac, processed_ac_f, indent=2)

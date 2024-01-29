import urllib.request

ass_sum_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
ass_sum_path = "bacteria_assembly_summary.txt"
with open(ass_sum_path, 'w') as ass_sum_f:
    for i, line in enumerate(urllib.request.urlopen(ass_sum_url)):
        if i == 0 and line.startswith(b'#') and b'assembly_accession' not in line:
            continue
        elif line.startswith(b'#'):
            ass_sum_f.write(line.strip(b'#').decode('utf-8'))
        else:
            ass_sum_f.write(line.decode('utf-8'))
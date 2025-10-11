from Bio import SeqIO
import os 


for file in os.listdir():
    if ("uniprot" in file) and ("viruses" not in file):
        with SeqIO.parse(file, "swiss") as org:
            with open("NonViral.fasta", "a") as nonviral:
                for record in org:
                    nonviral.write(f">Nonviral_{record.id}\n")
                    nonviral.write(f"{record.seq}\n")
    if ("uniprot" in file) and ("viruses" in file):
        with SeqIO.parse(file, "swiss") as org:
            with open("Viral.fasta", "a") as viral:
                for record in org:
                    viral.write(f">Viral_{record.id}\n")
                    viral.write(f"{record.seq}\n")




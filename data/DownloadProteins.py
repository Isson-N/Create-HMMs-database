from Bio import SeqIO
import os 
import pandas as pd

# Как оказалось, в uniprot семейства вирусов располагаются по разному (в разделе taxonomy нет единого индекса), поэтому придётся за этим следить...
df = pd.read_excel("ICTV_Master_Species_List_2024_MSL40.v2.xlsx", sheet_name="MSL")
families = df["Family"].unique()

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
                    family_name = ""
                    for i in record.annotations.get("taxonomy"):
                        if i in families:
                            family_name = i
                            break
                    if family_name == "":
                        family_name = "Unrecognized"                 
                    viral.write(f">{family_name}_{record.id}\n")
                    viral.write(f"{record.seq}\n")




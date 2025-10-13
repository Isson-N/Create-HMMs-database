params.script1 = "data/DownloadProteins.py"
params.script2 = "data/Stat.py"
params.model = ""


process ExtractProteins {
    input:
      path script
    
    output:
      path "*.fasta"

    script:
    """
    # Archaea
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.dat.gz
    gzip -d uniprot_sprot_archaea.dat.gz
    # Bacteria
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
    gzip -d uniprot_sprot_bacteria.dat.gz
    # Fungi
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_fungi.dat.gz
    gzip -d uniprot_sprot_fungi.dat.gz
    # Human
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz
    gzip -d uniprot_sprot_human.dat.gz
    # Invertebrates
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_invertebrates.dat.gz
    gzip -d uniprot_sprot_invertebrates.dat.gz
    # Mammals
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_mammals.dat.gz
    gzip -d uniprot_sprot_mammals.dat.gz
    # Plants
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz
    gzip -d uniprot_sprot_plants.dat.gz
    # Rodents
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_rodents.dat.gz
    gzip -d uniprot_sprot_rodents.dat.gz
    # Vertebrates
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_vertebrates.dat.gz
    gzip -d uniprot_sprot_vertebrates.dat.gz

    #Viruses
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz
    gzip -d uniprot_sprot_viruses.dat.gz

    python $script
    """
}


process DownloadDatabase {
    output:
      path "VOGs.hmm*"

    script:
    """
    wget https://fileshare.csb.univie.ac.at/vog/vog232/vog.hmm.tar.gz
    tar -xzf vog.hmm.tar.gz
    cat hmm/*.hmm > VOGs.hmm
    rm -r hmm
    hmmpress VOGs.hmm
    """
}

process HmmSearch {
    input:
      path hmm_database
      path proteins
    
    output:
      path "*.txt"

    script:
    """
    hmmsearch --tblout Viral_results.txt VOGs.hmm Viral.fasta
    hmmsearch --tblout NonViral_results.txt VOGs.hmm NonViral.fasta
    """
}

process Stat {
    input:
      path res
      path script

    output:
      path "*"

    script:
    """
    python $script
    """
}


workflow {
log.info """
===============================================================================
This pipeline use HMMs from VOGDB and test it on viral and non viral proteins
===============================================================================
"""
script1 = file(params.script1)
script2 = file(params.script2)


extractProteins_ch = ExtractProteins(script1)

if (!params.model) {
    downloadDatabase_ch = DownloadDatabase()
}
else {
    downloadDatabase_ch = Channel.fromPath("${params.model}/*")
}

resultsTxt_ch = HmmSearch(downloadDatabase_ch.collect(), extractProteins_ch.collect())
stat_ch = Stat(resultsTxt_ch, script2)
}

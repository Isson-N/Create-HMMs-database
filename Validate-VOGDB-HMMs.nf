params.script1 = "data/DownloadProteins.py"
params.script2 = "data/ViralNonviralStat.py"
params.script3 = "data/FamiliesStat.py"
params.msl = "data/ICTV_Master_Species_List_2024_MSL40.v2.xlsx"
params.model = ""
params.output = "/home/bessonov_id/work/Create-HMMs-database"
params.testFamilies = false

process ExtractProteins {
    input:
      path script
      path msl
    
    output:
      path "Viral.fasta", emit: viral
      path "NonViral.fasta", emit: nonviral

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
      path protein
    
    output:
      path "${protein.baseName}_results.txt"

    script:
    def baseName = protein.baseName
    """
    hmmsearch --tblout "${baseName}_results.txt" *.hmm $protein
    """
}

process ViralNonviralStat {
publishDir "${params.output}/ViralNonviralStat", mode: 'copy'
    input:
      path script
      path results

    output:
      path "*"

    script:
    """
    python $script
    """
}

process FamiliesStat {
publishDir "${params.output}/FamiliesStat", mode: 'copy'
    input:
      path script
      path viral_res

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

If default mode (only output was specified)
1. ExtractProteins - makes viral and nonviral fasta files from uniprot ftp
2. DownloadDatabase - download hmm models for each VOG from VOGDB
3. HmmSearch - hmmsearch models agains fasta files
4. ViralNonviralStat - produce distributions and ROC curves, how model divide viral and noviral sequences

If model was specified
1. ExtractProteins - makes viral and nonviral fasta files from uniprot ftp
2. DownloadDatabase - use specified model
3. HmmSearch - hmmsearch models agains fasta files
4. ViralNonviralStat - produce distributions and ROC curves, how models divide viral and noviral sequences

If model was specified and testFamilies = true
1. ExtractProteins - makes viral and nonviral fasta files from uniprot ftp. Next processes
use only viral fasta files. It is important, that family name specified in each record id
2. DownloadDatabase - use specified model
3. HmmSearch - hmmsearch models agains fasta files
4. FamiliesStat - produce distributions and ROC curves, how models divide each families
===============================================================================
"""
script1 = file(params.script1)
script2 = file(params.script2)
script3 = file(params.script3)
msl1 = file(params.msl)



if (!params.model && params.testFamilies == false)  {
    extractProteins_ch = ExtractProteins(script1, msl1)
    downloadDatabase_ch = DownloadDatabase()
    multi_ch = extractProteins_ch.viral.mix(extractProteins_ch.nonviral)
    resultsTxt_ch = HmmSearch(downloadDatabase_ch.collect(), multi_ch)
    stat_ch = ViralNonviralStat(script2, resultsTxt_ch.collect())
}
else if (params.model && params.testFamilies == false) {
    extractProteins_ch = ExtractProteins(script1, msl1)
    downloadDatabase_ch = Channel.fromPath("${params.model}/*")
    multi_ch = extractProteins_ch.viral.mix(extractProteins_ch.nonviral)
    resultsTxt_ch = HmmSearch(downloadDatabase_ch.collect(), multi_ch)
    stat_ch = ViralNonviralStat(script2, resultsTxt_ch.collect())
}
else if (params.model && params.testFamilies == true) {
    extractProteins_ch = ExtractProteins(script1, msl1)
    downloadDatabase_ch = Channel.fromPath("${params.model}/*")
    resultsTxt_ch = HmmSearch(downloadDatabase_ch.collect(), extractProteins_ch.viral)
    stat_ch = FamiliesStat(script3, resultsTxt_ch)

}

}

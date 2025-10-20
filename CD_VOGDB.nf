params.splitByFamilies = "./data/UsingVOGDB.py"
params.api_key = "b187bc11e44bbe5f49b434dee3b63ddeaf09"
params.configuration = "./data/configuration.cnf"
params.output = ""



// This process download protein sequences of groups
process DownloadVOGs {
    output:
      path "faa/*.faa"
      
    script:
    """
    wget https://fileshare.csb.univie.ac.at/vog/vog232/vog.faa.tar.gz
    tar -xzf vog.faa.tar.gz
    rm vog.faa.tar.gz
    """
}


// Thie process use python script to split each VOG to different families
process SplitByFamilies {
errorStrategy 'retry'
maxRetries 4
maxForks 5
    input:
      path VOG
      path script
      val api_key
      
    output:
      path "*.fasta", optional: true
      
    script:
    """
    python $script -f $VOG -a $api_key
    """
}


// This process align each protein fasta file 
process Alignment {
errorStrategy 'ignore'
    input:
      path seq
      
    output:
      path "${seq.baseName}_aligned.fasta"

    script:
    def baseName = seq.baseName
    """
    muscle -super5 $seq -output ${baseName}_aligned.fasta
    """
}


// This process run tabajara and create HMMs
process RunTabajara {
    input:
      path aligned
      path conf
      
    output:
      path "${aligned.baseName}.hmm", optional: true

    script:
    def baseName = aligned.baseName
    """
    cp $conf new_conf
    echo "input_file=$aligned" >> new_conf
    echo "output=cons" >> new_conf
    tabajara.pl -conf new_conf
    if [ -z "\$(ls -A cons/hmms/valid_HMMs)" ]; then
        echo "Папка пуста"
    else
        cat cons/hmms/valid_HMMs/* > ${baseName}.hmm
    fi
    """
}


// This process press all hmms to database
process PressHMMs {
publishDir "${params.output}/VOGDB_own_hmm", mode: 'copy'
    input:
      path hmms
      
    output:
      path "*"
      
    script:
    """
    cat $hmms > hmm_database.hmm
    hmmpress hmm_database.hmm
    """
}


workflow {
log.info """
=========================================
DESCRIPTION: Create HMMs database using proteis sequences of different viral families using VOGDB
=========================================
"""



script1 = file(params.splitByFamilies )
config1 = file(params.configuration)


downloadVogs_ch = DownloadVOGs()
splitByFamilies_ch = SplitByFamilies(downloadVogs_ch.flatten(), script1, params.api_key)
alignment_ch = Alignment(splitByFamilies_ch.flatten())
runTabajara_ch = RunTabajara(alignment_ch, config1)
pressHMMs_ch = PressHMMs(runTabajara_ch.collect())
}




params.splitByFamilies = "./data/UsingVOGDB.py"
params.api_key = "b187bc11e44bbe5f49b434dee3b63ddeaf09"
params.configuration = "./data/configuration.cnf"
params.output = ""

log.info """
=========================================
DESCRIPTION: Create HMMs database using proteis sequences of different viral families using VOGDB
=========================================
"""



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
maxForks 5
errorStrategy 'retry'
maxRetries 4
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
      path "*.hmm", optional: true

    script:
    """
    cp $conf new_conf
    echo "input_file=$aligned" >> new_conf
    echo "output=cons" >> new_conf
    tabajara.pl -conf new_conf
    if [ -d "cons/hmms/valid_HMMs" ]; then
        mv cons/hmms/valid_HMMs/ ./
    fi
    """
}


// This process press all hmms to database
process PressHMMs {
publishDir "${params.output}", mode: 'copy'
    input:
      path shmm
      
    output:
      "*"
      
    script:
    """
    cat $hmms > hmm_database.hmm
    hmmpress hmm_database.hmm
    rm hmm_database.hmm
    """
}


workflow {
script1 = file(params.splitByFamilies )
config1 = file(params.configuration)


downloadVogs_ch = DownloadVOGs()
splitByFamilies_ch = SplitByFamilies(downloadVogs_ch.flatten(), script1, params.api_key)
alignment_ch = Alignment(splitByFamilies_ch.flatten())
runTabajara_ch = RunTabajara(alignment_ch, config1)
pressHMMs_ch = PressHMMs(runTabajara_ch.collect())


}




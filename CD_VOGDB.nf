params.splitByFamilies = "/data/UsingVOGDB.py"
params.email = "ilabessonov9@gmail.com"
params.configuration = "data/configuration.cnf"
params.output = ""

log.info """
=========================================
DESCRIPTION: Create HMMs database using proteis sequences of different viral families using VOGDB
=========================================
"""



// This process download protein sequences of groups
process DownloadVOGs {
    output:
      path "*.faa"
      
    script:
    """
    wget https://fileshare.csb.univie.ac.at/vog/vog220/vog.faa.tar.gz
    tar -xzf vog.faa.tar.gz
    rm vog.faa.tar.gz
    """
}


// Thie process use python script to split each VOG to different families
process SplitByFamilies {
    input:
      path VOG
      path script
      val email
      
    output:
      path ".fasta", optional: true
      
    sceipt:
    """
    python $script -f $VOG -e $email
    """
}


// This process align each protein fasta file 
process Alignment {
    input:
      path seq
      
    output:
      path "{basename}_aligned.fasta"

    script:
    def baseName = seq.baseName
    """
    muscle -align $seq -output ${baseName}_aligned.fasta
    """
}


// This process run tabajara and create HMMs
process RunTabajara {
    input:
      path aligned
      path conf
      
    output:
      path "*hmm"

    script:
    """
    echo "input_file=$aligned" >> $conf
    echo "output=cons" >> $conf
    tabajara -conf $conf
    cd cons/hmms/valid_HMMs/
    """
}


// This process press all hmms to database
process PressHMMs {
publishDir "${params.output}", mode: 'copy'
    input:
      path hmm
      
    output:
      "*"
      
    script:
    """
    cat hmm > hmm_database.hmm
    hmmpress hmm_database.hmm
    rm hmm_database.hmm
    """
}


workflow {
downloadVogs_ch = DownloadVOGs()
splitByFamilies_ch = SplitByFamilies(downloadVogs_ch.flatten(), params.splitByFamilies, params.email)
alignment_ch = Alignment(splitByFamilies_ch)
runTabajara_ch = RunTabajara(alignment_ch, params.configuration)
pressHMMs_ch = PressHMMs(runTabajara_ch.collect())


}





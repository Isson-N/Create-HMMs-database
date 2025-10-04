params.parse_msl = '/home/kali/Desktop/DownloadsWork/data/ParseMSL.py'
params.msl = '/home/kali/Desktop/DownloadsWork/data/ICTV_Master_Species_List_2024_MSL40.v2.xlsx'


log.info """
=========================================
DESCRIPTION: Create HMMs database using proteis sequences of different viral families using GenBank and ICTV
=========================================
"""


// This process use python script, for extracting all familes from MSL file from ICTV
process ParseMSL {
    input:
      path script
      path msl
    
    output:
      path "Families.csv"
      
    script:
    """
    python $script
    """
}


// This process use Family name to make file with string "{Family_name}; {List_with_proteins}"
process GetProteinsNames {
    input:
      val family
      
    output:
      path "*"
      
    script:
    """
    ### Сделай 1> и 2> А потом python скрипт. По умолчанию обычный поток идёт на аутпут, а поток ошибок в папку publishdir. Надо ещё обрабатывать пустые файлы с потоком ошибок
    ### Скачивай сразу через download. Если не скачивает и это видит питон скрипт, то ничего не делает. Если же не пустой, то тогда парсит.
    """
}











workflow {



families_ch = ParseMSL(params.parse_msl, params.msl).splitCsv().map{ item -> item[0] }
strings_ch = GetProteinsNames(families_ch)

}




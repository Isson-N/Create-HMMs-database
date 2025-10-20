import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, RocCurveDisplay
import matplotlib.pyplot as plt
import numpy as np


# Функция, чтобы парсить таблицы, полученные из hmmsearch
def parser(file_path):
    col_names = ['target_name', 'target_accession', 'query_name', 'query_accession', 
                 'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
                 'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description']
    df = pd.read_csv(file_path, comment='#', sep=r'\s+', names=col_names, header=None)
    return df

def add_family(df):
    real_families = []
    detected_families = []
    new_df = df
    for name in df["target_name"].to_list():
        real_families.append(name.split("_")[0].replace(">", ""))
    for name in df["query_name"].to_list():
        detected_families.append(name.split("_")[0])
    new_df.insert(len(new_df.columns), 'real_family', real_families)
    new_df.insert(len(new_df.columns), 'detected_family', detected_families)
    return new_df 
        

# Применяем, а ещё добавляем инфу о том, вирусный или нет
viral_df = parser("Viral_results.txt")
viral_df = add_family(viral_df)
print(viral_df["detected_family"].unique())



def parse_model(args):
	model, viral_df = args
	print(model, end="\r")
	dfm = combined_df[combined_df["query_name"] == model]
 
	scores = dfm['full_score']
	labels = dfm['is_viral']
	fpr, tpr, thresholds = roc_curve(labels, scores)
	roc_auc = auc(fpr, tpr)
	if len(dfm[dfm["is_viral"] == 0]) == 0:
		roc_auc = 1
	
	line = []
	line.append(model)
	line.append(roc_auc)

	try:
		viral_max_full_score = np.max(dfm[dfm["is_viral"] == 1]["full_score"])
	except:
		viral_max_full_score = ""
	line.append(viral_max_full_score)
    
	try:
		nonviral_max_full_score = np.max(dfm[dfm["is_viral"] == 0]["full_score"])
	except:
		nonviral_max_full_score = ""
	line.append(nonviral_max_full_score)
    
	try:
		relation = viral_max_full_score/nonviral_max_full_score
	except:
		relation = ""
	line.append(relation)
    
	try:
		mean_only_viral = np.mean(dfm[(dfm["is_viral"] == 1) & (dfm["full_score"] >= nonviral_max_full_score)]["full_score"])
	except:
		mean_only_viral = ""
	line.append(mean_only_viral) 

	try:
		viral_count = len(dfm[dfm["is_viral"] == 1])
	except:
		viral_count = ""
	line.append(viral_count)

	try:
		nonviral_count = len(dfm[dfm["is_viral"] == 0])
	except:
		nonviral_count = ""
	line.append(nonviral_count)

	try:
		only_viral = len(dfm[(dfm["is_viral"] == 1) & (dfm["full_score"] >= nonviral_max_full_score)])
	except:
		only_viral = ""
	line.append(only_viral)
	return line


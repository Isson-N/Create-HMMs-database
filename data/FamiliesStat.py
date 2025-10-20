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
    new_df['real_family'] = real_families
    new_df['detected_family'] = detected_families
    return new_df 
        
def parse_model(args):
	family, viral_df = args
	dfm = viral_df[viral_df["real_family"] == family].copy()

	dfm["positive"] = (dfm["real_family"] == dfm["detected_family"]).astype("int")

	scores = dfm['full_score']
	labels = dfm['positive']
	fpr, tpr, thresholds = roc_curve(labels, scores)
	roc_auc = auc(fpr, tpr)
	if len(dfm[dfm["positive"] == 0]) == 0:
		roc_auc = 1

	plt.plot(fpr, tpr, color='blue', linewidth=0.5)

	line = []
	line.append(family)
	line.append(roc_auc)

	models_count = len(dfm)
	line.append(models_count)

	true_detect_count = len(dfm[dfm["positive"] == 1])
	line.append(true_detect_count)

	false_detect_count = len(dfm[dfm["positive"] == 0])
	line.append(false_detect_count)
 
	max_correct_full_score = np.max(dfm[dfm["positive"] == 1]["full_score"])
	line.append(max_correct_full_score)

	max_uncorrect_full_score = np.max(dfm[dfm["positive"] == 0]["full_score"])
	line.append(max_uncorrect_full_score)

	return line
 
	
	
	

# На этот раз обойдёмся без multiprocessing...
# Применяем, а ещё добавляем инфу о том, вирусный или нет
viral_df = parser("Viral_results.txt")
viral_df = add_family(viral_df)

results_df = pd.DataFrame(columns = ["Family", "AUC", "Models_count", "True_detect_count", "False_detect_count", "Max_correct_full_score", "Max_uncorrect_full_score"])

fam_cnt = len(viral_df["real_family"].unique())
k = 0
for family in viral_df["real_family"].unique().tolist():
	ln = parse_model((family, viral_df))
	results_df.loc[len(results_df)] = ln
	k+=1
	print(f"\rВыполнил для семьи: {family}. Это {k} семья из {fam_cnt}")
 
 
results_df.to_csv("Families_models.csv")
plt.savefig("Families_ROC_curves.png")
plt.show()

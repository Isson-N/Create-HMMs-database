import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, RocCurveDisplay
import numpy as np


# Функция, чтобы парсить таблицы, полученные из hmmsearch
def parser(file_path):
    col_names = ['target_name', 'target_accession', 'query_name', 'query_accession', 
                 'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
                 'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description']
    df = pd.read_csv(file_path, comment='#', sep=r'\s+', names=col_names, header=None)
    return df

# Применяем, а ещё добавляем инфу о том, вирусный или нет
viral_df = parser("Viral_results.txt")
nonviral_df = parser("NonViral_results.txt")
viral_df["is_viral"] = 1
nonviral_df["is_viral"] = 0


# Строим гистограмму для full_score
figure, axis = plt.subplots(1, 2, figsize=(15, 6))
axis[0].hist([viral_df["full_score"], nonviral_df["full_score"]],
        bins=100,
        alpha=0.5,
        label=['Viral', 'NonViral'],
        color=['#50e991', '#EE204D'])
axis[0].set_yscale('log')
axis[0].set_title('Distribution of full_score in Viral sequences and NonViral sequences')
axis[0].legend(['Viral', 'NonViral'])
axis[0].set_ylabel('Number of sequences')
axis[0].set_xlabel('Full score')
axis[0].grid(True, alpha=0.3)

# Строим ROC-кривую для full_score
combined_df = pd.concat([viral_df, nonviral_df], ignore_index=True)
scores = combined_df['full_score']
labels = combined_df['is_viral']
fpr, tpr, thresholds = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

axis[1].plot(fpr, tpr, color='#0bb4ff', linewidth=2, 
                label=f'full_score (AUC = {roc_auc:.3f})')
axis[1].plot([0, 1], [0, 1], 'k--', alpha=0.5)

optimal_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[optimal_idx]
optimal_fpr = fpr[optimal_idx]
optimal_tpr = tpr[optimal_idx]
axis[1].scatter(optimal_fpr, optimal_tpr, color='red', s=100, zorder=5,
           label=f'Optimal threshold = {optimal_threshold:.3f}\n(TPR = {optimal_tpr:.3f}, FPR = {optimal_fpr:.3f})')

axis[1].set_xlabel('False Positive Rate')
axis[1].set_ylabel('True Positive Rate')
axis[1].set_title('ROC Curve for Full Score')
axis[1].legend()
axis[1].grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("full_score_distribution_and_ROC.png")
plt.show()
print("Сохранил график full_score_distribution_and_ROC")


# Построим ROC кривые для каждой модели
models = pd.DataFrame(columns = ['Model_name', "AUC", "Viral_max_full-score", "NonViral_max_full-score",  "Relation_Viral-NonViral", "Mean_only_viral_full-score" ])
k = 0
ln = len(combined_df["query_name"].unique())
for model in combined_df["query_name"].unique():
	dfm = combined_df[combined_df["query_name"] == model]
	scores = dfm['full_score']
	labels = dfm['is_viral']
	fpr, tpr, thresholds = roc_curve(labels, scores)
	roc_auc = auc(fpr, tpr)
 
	line = []
	line.append(model)
	line.append(roc_auc)
	
	try:
		viral_max_full_score = np.max(dfm[dfm["is_viral"] == 1]["full_score"])
	except:
		viral_max_full_score = 0
	line.append(viral_max_full_score)
    
	try:
		nonviral_max_full_score = np.max(dfm[dfm["is_viral"] == 0]["full_score"])
	except:
		nonviral_max_full_score = 0
	line.append(nonviral_max_full_score)
    
	try:
		relation = viral_max_full_score/nonviral_max_full_score
	except:
		relation = 0
	line.append(relation)
    
	try:
		mean_only_viral = np.mean(dfm[(dfm["is_viral"] == 1) and (dfm["full_score"] >= nonviral_max_full-score)]["full_score"])
	except:
		mean_only_viral = 0
	line.append(mean_only_viral) 
	
	models.loc[len(models)] = line
  	
	plt.plot(fpr, tpr, color='#0bb4ff', linewidth=1, 
                label=f'full_score (AUC = {roc_auc:.3f})')
	
	k+=1
	print(f"\rОбработал модель: {model}. Это {k} модель из {ln}", end="")

models.to_csv("models.csv", index=False)
filt_models = models.dropna(subset=['AUC'])
filt_models = models[models["AUC" >= 0.75]]
filt_models.to_csv("good_models.csv", index=False)

plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for Full Score')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("ROC_curves_for_models.png")
plt.show()
    
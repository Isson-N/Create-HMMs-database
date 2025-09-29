import pandas as pd

msl = pd.read_excel('ICTV_Master_Species_List_2024_MSL40.v2.xlsx', sheet_name='MSL')
families = pd.DataFrame(msl['Family'].dropna().unique())

families.to_csv('Families.csv', index=False, header=False)

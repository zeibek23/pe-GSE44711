import pandas as pd

expression = pd.read_csv("results/01_expression.csv", index_col=0)
diagnosis = pd.read_csv("results/01_diagnosis.csv", index_col=0).squeeze()

print(expression.shape)
print(expression.head(5))
print(type(diagnosis))
print(len(diagnosis))
print(all(diagnosis.index == expression.columns))


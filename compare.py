from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv('ref.dat',delim_whitespace=True)
df2 = pd.read_csv('comparar.dat')
df1 = pd.read_csv('curva.dat')
df.columns=['ang','pot']
df2.columns=['ang','pot']
df1.columns=['ang','pot']
df['dataset']='REF'
df2['dataset']='TEST'
df1['dataset']='CURVE'
df3=pd.concat([df,df2])
df3=pd.concat([df3,df1])
print(df2)
sns.lineplot(data=df3,x='ang',y='pot',hue='dataset')
plt.show()

#! /usr/bin/env python3
import pandas as pd

table=[]
with open("./d16-200.table-cs.hDijet.txt","r") as fi:
    for i in fi.readlines():
        #row=[j.strip() for j in i.split()]
        row= i.split()
        table.append(row)
df=pd.DataFrame(table[1:len(table)], columns=table[0])
#df=pd.read_csv("../d16-200.table-cs.hDijet.txt",sep='\t')

minimum=df[['q2min','q2max','Pt_min', 'Pt_max','Sigma', 'tot+(%)','tot-(%)']]

#print(minimum['q2min'])
#print(minimum['q2min'].iloc[0])
for i,q2min in enumerate(minimum['q2min'].unique()):
    print(i,"  ", q2min)
    #sub=minimum.query('q2min=={}'.format(q2min))
    sub=minimum[minimum['q2min']==q2min]
    print(sub.head())
    
    sub['Sigma/dq2dpt']=sub.apply(lambda row: float(row['Sigma'])/((float(row['q2max'])-float(row['q2min']))*(float(row['Pt_max'])-float(row['Pt_min']))),axis=1)
    print(sub.head())
    #sub[['Pt_min', 'Pt_max','Sigma', 'tot+(%)','tot-(%)','Sigma/dq2dpt']]
    with open('subtable-{}.txt'.format(i+1),"w") as fi:
        fi.write(sub[['Pt_min', 'Pt_max','Sigma', 'tot+(%)','tot-(%)','Sigma/dq2dpt']].to_string(header=False,index=False))


#print(df.columns)
#print(df.query('float(q2min)==5.5').head())
#for i in pd.unique(df['q2min']):
#    print(df.query('q2min=={}'.format(i)).head())

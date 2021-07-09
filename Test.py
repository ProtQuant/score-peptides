import heapq
import os
import pandas as pd

dirname = '.\\1_findOutput'
# for maindir, subdir, file_name_list in os.walk(dirname):
#     print(maindir)
#     print(subdir)
#     print(file_name_list)
#     print("******************")

# filelist=os.listdir(dirname)

# print(filelist)

# list = [1,2,3]
# print([list])
from joblib.numpy_pickle_utils import xrange

# df = pd.DataFrame({"pepIndex":[[1,2], [2, 3, 4]],
#                     'count':[0,0]})
# print(df)
# df['count'] = df.apply(lambda x: len(x['pepIndex']), axis=1)
# print(df)
# df.loc[0,'pepIndex'].append(8)
# print(df)
# # df['new'] = [[3,5],[4]]
# new = [[] for _ in xrange(len(df))]
# print(new)
# df['new'] = new
# df.loc[0,'new'].append(9)
# print(df)

# a = [3,4,5,5]
#
# print(sum(heapq.nlargest(3,a)))

# df1 = pd.DataFrame(columns=['pep' ,'in'])
# print(df1)
# # df2 = pd.DataFrame({'pep': ['aa'], 'in': [44]})
# # print(df2)
# # df3 = df2.append(df1, ignore_index=True)
# # print(df3)
# df4 = df1.append([{'pep': "ss", 'in': 6}], ignore_index=True)
# print(df4)
# df5 = df4.append([{'pep': "ss", 'in': 6}], ignore_index=True)
# print(df5)

# d1 = {"a":[1], "b":[7]}
# d2 = {"a":[2]}
# d1.update(d2)
# print(d1)

# pep = ['a', 'c', 'b', 'c', 'c', 'c', 'b']
# ind = [50, 50, 50, 100, 100, 100, 100]
# df = pd.DataFrame({'pep': pep, 'ind': ind})
# print(df.value_counts())
# print(df.value_counts().index[0])
# print(df.groupby(['pep']))
# df1 = df.groupby(['pep'])['ind'].agg(lambda x: x.value_counts().index[0])  #sort=False
# print(df1)

# df_a = pd.DataFrame({'pep': ['aa','ss','bb'], 'ind': [0, 0, 0]})
# df_b = pd.DataFrame({'pep': ['ss', 'cc'], 'ind': [7, 8]})
# print(df_a)
# print(df_b)
# print(df_a['pep'].isin(df_b['pep']))
# df_a_filter = df_a[~ df_a['pep'].isin(df_b['pep'])]
# print(df_a_filter)
# print(df_b)
# df_c = pd.DataFrame(df_a, columns=['pep', 'ind'])
# print(df_c)

proteinInfoFolder = ".\\t1\\t2"
if not os.path.exists(proteinInfoFolder):
    os.makedirs(proteinInfoFolder)

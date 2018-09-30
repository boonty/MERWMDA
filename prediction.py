# encoding: UTF-8
import xlrd
import cmath
import math
from math import e
import numpy as np
import numpy.linalg as LA

N = 495
M = 383
'''

def Transition_Matrix(Relation_Matrix):
    row = np.shape(Relation_Matrix)[0]
    col = np.shape(Relation_Matrix)[1]
    Transition_P = np.zeros((row,col))
    value, vector = LA.eig(Relation_Matrix)
    value = list(value)
    value = [item.real for item in value]
    max_lambda = max(value)
    max_index = value.index(max(value))
    max_vector = vector[:, max_index]
    max_vector = [u.real for u in list(max_vector)]
    temp = sum(max_vector)
    norm_max_vector = [x ** 0.5 / temp ** 0.5 for x in max_vector]
    for i in range(row):
        for j in range(col):
            if norm_max_vector[i] == 0.:
                Transition_P[i][j] = 0
            else:
                Transition_P[i][j] = (Relation_Matrix[i][j] * norm_max_vector[j]) / (max_lambda * norm_max_vector[i])
        Transition_P[i] = Transition_P[i]/(np.sum(Transition_P[i]))
    return Transition_P
def Prediction_MiRNA(Interaction,Transtion):
    row = np.shape(Interaction)[0]
    col = np.shape(Interaction)[1]
    Predictionmatrix=np.zeros((row,col))
    Final_Prediction = np.zeros((N,M))
    for k in range(row):
        error = 1000000
        der = sum(Interaction[k])
        if der == 0:
            der += 1
        else:
            der = der
        initial = Interaction[:,k] / der
        count=0
        while error > 1e-6:
            temp_initial=np.dot(np.transpose(Transtion),initial)
            error = LA.norm(temp_initial-initial,1)
            initial = temp_initial
            count+=1
        Predictionmatrix[k] = initial
        P_1 = Predictionmatrix[0:383,383:col]
        P_2 = Predictionmatrix[383:row,0:383]
        Final_Prediction = (np.transpose(P_1) + P_2) * 0.5
    return Final_Prediction

Adj_Matrix = np.zeros((N,M))     # 邻接矩阵
with open('knowndiseasemirnainteraction.txt') as File1:
    for line in File1:
        arr = map(int, line.split())
        arr = list(arr)
        if arr[1] != 0:
            Adj_Matrix[arr[0] - 1][arr[1] - 1] = 1
Mirna_Matrix = []            # RNA功能联系矩阵
with open('weight_mirna.txt') as File2:
    lines = File2.readlines()
    for j in range(N):
        line = [int(x) for x in lines[j].split()]
        Mirna_Matrix.append(line)
    Mirna_Matrix = np.array(Mirna_Matrix)
    for k in range(N):
        Mirna_Matrix[k][k] = 0.
Disease_Matrix = []           # 疾病语义联系矩阵
with open('weight_disease.txt') as File3:
    lines = File3.readlines()
    for j in range(M):
        line = [int(x) for x in lines[j].split()]
        Disease_Matrix.append(line)
    Disease_Matrix = np.array(Disease_Matrix)
    for k in range(M):
        Disease_Matrix[k][k] = 0.
temp1 = np.hstack((Disease_Matrix, np.transpose(Adj_Matrix)))
temp2 = np.hstack((Adj_Matrix, Mirna_Matrix))
Hyper_Matrix = np.vstack((temp1,temp2))    # 异质网络邻接矩阵


Transition_P_matrix = Transition_Matrix(Hyper_Matrix)
print(Transition_P_matrix)
print(np.sum(Transition_P_matrix,axis =1))
Prediction_as_MiRNA = Prediction_MiRNA(Hyper_Matrix,Transition_P_matrix)
print (Prediction_as_MiRNA)
print (np.shape(Prediction_as_MiRNA))

data =xlrd.open_workbook(r'diseasenumbers.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
disease_data=[]
for row_index in range(table.nrows):
    disease_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])
print (disease_data)

data =xlrd.open_workbook(r'miRNAnumbers.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)

miRNA_data=[]
for row_index in range(table.nrows):
    miRNA_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])
print (miRNA_data)
f7w=open(r'results_for_all.txt','w')
str3=''
for i in range(N):
    for j in range(M):
        if Adj_Matrix[i][j]==0:
            str3+=str(disease_data[j][1])+'\t'+str(miRNA_data[i][1])+'\t'+str(Prediction_as_MiRNA[i][j].real)+'\n'
f7w.write(str3)
f7w.close()

print ('ok')
'''
data =xlrd.open_workbook(r'results_for_all.xlsx')
table=data.sheets()[0]
x=table.nrows
for row_index in range(table.nrows):
    f=open(r'%s.txt'%table.row(row_index)[0].value.encode('utf-8'),'a')
    f.write(str(table.row(row_index)[0].value.encode('utf-8'))+'\t'+str(table.row(row_index)[1].value.encode('utf-8'))+'\t'+str(table.row(row_index)[2].value)+'\n')
print ('okky')
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 19:08:33 2018

@author: Konstantin

"""
from sympy import *
import re
from sympy.parsing.sympy_parser import parse_expr

# f - main function 
# fSign - sign of main function 
# limsList - list of limiting inequalities
# x - tuple of symbolic variables 
# diffs - dict of f's diffs {var : f.diff(var)}       
# diffs2 - 2-nd diff (see 'diff')
# statPoints - stationary points of f function
# matrixList - a list of matrix that actually determine Guesse matrix 
# detsGes - a list of Guesse matrix determinants
# lagFun - Lagrange function
# w - tuple of additional nonnegative variables of lagVun
# lagVars - tuple of all variables of lagFun
# lagDiffs - dict of lagFun's diffs
# kunSys - a list of Kuhn–Tucker conditions
# ansList - a list of points that solve knuSys
# correctAnsList - a list of ponts from ansList that fulfil list of limiting inequalities


inputStrFun = input('Введите целевую функцию: f=')
print('Для выбора максимизации/минимизации\nвведите max или min соответственно: \n')
while True:                                 # ввод  максимизации/минимизации функции
    fSign = input()
    if fSign == 'min' or fSign == 'max':
        break
    print('Input Error')

limsList = []
print('Введите ограничеия:\n(ограничение) >= 0.\nПустая строка - окончание ввода.')
while True:
    curLim = input()                    #ввод ограничевающих неравенств
    if curLim == '' or curLim == ' ':
        break
    limsList.append(parse_expr(curLim))
    
    
varSet = set(re.findall('x\d+', inputStrFun))    # извлекаем список переменных с помощью регулярных выражений
x = symbols(tuple(varSet))                       # создаем картеж переменных символьного (алгебраического) типа 

if fSign == 'max':
    f = parse_expr(inputStrFun)                      # парсим функцию (приводим у нужному типу и упрощаем)
else:
    f = -1*parse_expr(inputStrFun)

try:
    diffs = {v: f.diff(v) for v in x}                # дифференциируем функцию по всем переменным и кладем полученные занчения в словарь
except:
    print('Функция не дифференцируема.')
    raise SystemExit(1)

print('\nЦелевая функция: ',f,'\nПервые производные: ')     # вывод производных и целевой функции
for key in diffs.keys():
    print('df/d', key, ' = ', diffs[key], sep = '')

print('\nВторые производные: ')
diffs2 = {}
for v1 in diffs.keys():
    for v2 in diffs.keys():
        diffs2[S(str(v1) + str(v2))] = f.diff(v1,v2)        # создаем словарь вторых производных 
        print('d2f/'+'d'+str(v1)+'d'+str(v2)+' = ',f.diff(v1,v2) ,sep='')   # выводим его

try:
    statPoints = solve(diffs.values(), x)       # приравнимаем производные к 0 -> находим стацонарные точки
except:
    print('Нет стационарных точек.')
    raise SystemExit(1)

try:
    statPoints = [tuple(statPoints.values())]   # преобразование ответов к нужной стректуре
except:
    pass

print('\nСтационарные точки с координатами в формате т.',tuple(diffs.keys()),': ',sep='')  # вывод стационарных точек
for i in statPoints:
    print('т.',i,sep='')

varList = list(varSet)
varList.sort()

matrixList = [[diffs2[S(varList[0] + varList[0])]]]   
for i in range(2, len(varList)+1):                  # Формирование списка матриц
    ges = []
    for j in range(1, i+1):
        row = []
        for k in range(1, i+1):
            row.append(diffs2[S(varList[j-1] + varList[k-1])])
        ges.append(row)
    matrixList.append(ges)
            
matrixList = [Matrix(m) for m in matrixList]

detsGes = [x.det() for x in matrixList]     # список определителей матрицы Гессе
print('\n',matrixList[-1], sep = '')
print('Список определителей: ', detsGes)

detsGesInt = []
try:
    detsGesInt = [int(x) for x in detsGes]          # вычисление знака матрицы
    isPos = all(val > 0 for val in detsGesInt)
    isNeg = all(val > 0 for val in detsGesInt[1::2]) and all(val < 0 for val in detsGesInt[::2]) 
    
    if isPos == True:
        print('Матрица положительно определённая.\nФункция выпуклая.')
    elif isNeg == True:
        print('Матрица отрицательно определённая.\nФункция вогнутая.')
    else:
        print('Матрица не знакоопределённая')
except:
    print('Не могу определить знак.')

lagIneqs = []
w = symbols('w:'+ str(len(limsList)), nonnegative = True)      # список доп. переменных при ограничениях для ф-ии Лагранжа
lagVars = tuple(list(x) + list(w))        # полный список переменных для ф-ии Лагранжа
lagFun = f

for i in range(len(w)):
    curIneq = w[i]*limsList[i]
    lagIneqs.append(curIneq)
    lagFun += curIneq           # составление ф-ии Лагранжа
print('\nФункия Лагранжа: ', lagFun)

lagDiffs = {v: lagFun.diff(v) for v in x }  # нахождение произвдных функции Лагранжа 
kunSys =  list(lagDiffs.values()) + lagIneqs
ansList = solve(kunSys,lagVars)

ineqs = [S(str(v) + ' >= 0') for v in limsList]
correctAnsList = []

for ans in ansList :        # Проверка условия: содержатся ли точки в заданной области ограничения
    correctAns = True
    for ineq in ineqs:
        correctAns = correctAns and ineq.subs([(x[i],ans[i]) for i in range(len(x))])
    if correctAns == True:
        correctAnsList.append(ans)
    
print('Ответ',lagVars,': ',correctAnsList, sep='')









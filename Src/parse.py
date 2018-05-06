import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr,spearmanr
import numpy as np
import math


def calculateContact():
    bin_size = 500000  # Megabyte
    data = pd.read_table('data/ch7.txt', sep=' ', header=None)
    # print(data.head())

    connections = data[[2, 4]]
    connections.columns = [1, 2]

    sorted_data = connections.sort_values(by=[1])
    # bin the data
    binned_data = sorted_data // bin_size
    # start at bin 1
    binned_data = binned_data + 1
    # remove connection to self
    binned_data = binned_data[binned_data[1] != binned_data[2]]
    # remove connects that are already found
    binned_data = binned_data[binned_data[1] < binned_data[2]]
    # print(binned_data.head())
    out_data = binned_data.groupby([1, 2]).size()
    out_data.to_csv('chr_7_500kb.txt', sep='\t')
    # print(sorted_data.head())
    exit()


def Evaluation(wishMatrix, recontructedDistance):

    small_number = 0.1
    dim = len(wishMatrix)
    wishDistance = [dim][dim]
    correlation = []
    while small_number < 5:

        for i in range(dim):
            for j in range(dim):
                wishDistance[i][j] = 1/(wishMatrix[i][j] ** small_number)

        correlation.append(spearmanr(wishDistance, recontructedDistance, axis=None))
        small_number += .2
    return correlation


def normalization(matric):
    dim = len(matric)
    total = sum([sum(i) for i in matric])
    rowSum = [sum(matric[i]) for i in range(dim)]
    result = [dim][dim]
    for i in range(dim):
        for j in range(dim):
            result[i][j] = result[j][i] = matric[i][j]/(rowSum[i] * rowSum[j])*total

    return result


def CaDistanceMatrix(fileName):
    A = []
    for line in open(fileName):
        coordinateList = line.split()
        id = coordinateList[0]
        if id == 'ATOM':
            type = coordinateList[2]
            if type == 'CA':
                A.append([float(coordinateList[5]), float(coordinateList[6]), float(coordinateList[7])])

    sizePoints = len(A)
    distMatrix = np.zeros((sizePoints, sizePoints))
    for i, eachFir in enumerate(A):
        for j, eachSec in enumerate(A):
            print(PointDistance(eachFir, eachSec))
            distMatrix[i][j] = PointDistance(eachFir, eachSec)
            distMatrix[i][j] = distMatrix[j][i]

    return distMatrix


def PointDistance(pointFir, pointSec):
    new_array = [(pointFir[0]-pointSec[0])**2,(pointFir[1]-pointSec[1])**2, (pointFir[2]-pointSec[2])**2]
    dist = math.sqrt(new_array[0]+new_array[1]+new_array[2])
    return dist


reconstructedDis = CaDistanceMatrix("../LorDG/output/chromo_7_1MB_1525635742467.pdb")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == '__main__':
    binned_data = pd.read_csv('../Results/chromo_7_1MB.txt', sep='\t', header=None)
    x = binned_data.iloc[:, 0].copy().values
    y = binned_data.iloc[:, 1].copy().values
    matrix_data = binned_data.iloc[:, 2].copy().values

    size = max(max(x),max(y))
    print(size+1)

    list_matrix = [[0] * (size+1) for i in range(size+1)]


    plt.imshow(list_matrix, cmap='hot', interpolation='nearest')
    plt.show()
    i = 0

    while i < len(x):
        string = str(x[i]) + '\t' +str(y[i]) + '\t' + str(matrix_data[i])
        print(string)
        list_matrix[x[i]][y[i]] = matrix_data[i]

        i += 1

    plt.imshow(list_matrix, cmap='hot', interpolation='nearest')
    plt.show()


    #initial_matrix = binned_data.values

    #print(initial_matrix)


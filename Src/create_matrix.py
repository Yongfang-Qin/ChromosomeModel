import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from functools import reduce

import seaborn as sb

def normalization(matric):
    dim = len(matric)
    rowSum = [reduce(lambda x, y: x + y, item) for item in matric]
    #rowSum = sum(sum(matric[i]) for i in range(dim))
    total = sum([sum(i) for i in matric])
    print(total)
    result = [[None] * dim for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            #print(i, '\t', j, '\t', dim)

            if rowSum[i] != 0 and rowSum[j] != 0:
                result[i][j] = (matric[i][j]/(rowSum[i] * rowSum[j]))
            else:
                result[i][j] = 0

    return result

if __name__ == '__main__':
    binned_data = pd.read_csv('../Results/chromo_7_1MB.txt', sep='\t', header=None)
    x = binned_data.iloc[:, 0].copy().values
    y = binned_data.iloc[:, 1].copy().values
    matrix_data = binned_data.iloc[:, 2].copy().values

    size = max(max(x),max(y))
    print(size+1)

    list_matrix = [[0] * (size+1) for i in range(size+1)]

    i = 0

    while i < len(x):
        string = str(x[i]) + '\t' +str(y[i]) + '\t' + str(matrix_data[i])
        print(string)
        list_matrix[x[i]][y[i]] = matrix_data[i]

        i += 1





    for i in range(0, (int(len(list_matrix)))):
        for j in range(i + 1, len(list_matrix)):
            list_matrix[j][i]= list_matrix[i][j]

    out_matrix = np.asanyarray(list_matrix)
    # Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, icefire, icefire_r, inferno, inferno_r, jet, jet_r, magma, magma_r, mako, mako_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, rocket, rocket_r, seismic, seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r, vlag, vlag_r, winter, winter_r

    plt.imshow(out_matrix, cmap='BrBG')
    plt.savefig('../Results/contact_matrix_before_norm.png')
    #plt.show()
    plt.close()
    print(out_matrix)


    norm_matrix = normalization(out_matrix)

    plt.imshow(norm_matrix, cmap='BrBG')
    plt.savefig('../Results/contact_matrix_after_norm.png')
    #plt.show()
    plt.close()
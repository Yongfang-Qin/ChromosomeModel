import pandas as pd
from tqdm import tqdm

bin_size = 500000  # Megabyte
data = pd.read_table('data/ch7.txt', sep=' ', header=None)
#print(data.head())

connections = data[[2, 4]]
connections.columns = [1, 2]


sorted_data = connections.sort_values(by=[1])
# bin the data
binned_data = sorted_data//bin_size
# start at bin 1
binned_data = binned_data + 1
# remove connection to self
binned_data = binned_data[binned_data[1]!=binned_data[2]]
# remove connects that are already found
binned_data = binned_data[binned_data[1] < binned_data[2]]
#print(binned_data.head())
out_data = binned_data.groupby([1,2]).size()
out_data.to_csv('chr_7_500kb.txt', sep='\t')
#print(sorted_data.head())
exit()

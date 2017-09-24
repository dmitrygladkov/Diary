import numpy as np
import urllib.request

url = "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"

# Setup your proxy if any
#proxy_support = urllib.request.ProxyHandler({"https":"..."})
#opener = urllib.request.build_opener(proxy_support)
#urllib.request.install_opener(opener)

raw_data = urllib.request.urlopen(url)
X = np.loadtxt(raw_data, delimiter = ",", usecols = range(4)) # Download the first 4 columns
raw_data = urllib.request.urlopen(url)
y = np.loadtxt(raw_data, delimiter = ",", dtype = str, usecols = [4]) # Download last column
N, d = X.shape

print(N, d) # Number of objects and number of features
print(X)
print(y)
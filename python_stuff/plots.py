import matplotlib.pyplot as plt 
import sys


if len(sys.argv) > 1:
    file_name = sys.argv[1]
else:
    file_name = "performance.txt"
    
    
X = list()
Y = list()
with open(file_name, 'r') as f:
    for line in f:
        x, y = line.split()
        X.append(int(x))
        Y.append(float(y))
plt.plot(X, Y)
plt.xlabel('Число процессов, шт')
plt.ylabel('Время выполнения, c')
plt.savefig("performance.png")
plt.show()
    
    
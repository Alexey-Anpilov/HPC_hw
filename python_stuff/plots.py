import matplotlib.pyplot as plt 

file_name = '../Gram-Schmidt/performance.txt'
X = list()
Y = list()

if __name__ == "__main__":
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
    
    
import sys
import matplotlib.pyplot as plt

def parse(filename):
    N = []
    TIME = []
    MEM = []
    K = []
    for line in open(filename):
        tokens = line.split(',')
        TIME.append(float(tokens[1]))
        MEM.append(float(tokens[2]) / 1000) # GB
        K.append(int(tokens[3]))
        N.append(int(tokens[4]))
    return N, TIME, MEM, K

def do_scatter_plot(X1, Y1, X2, Y2, label1, label2, xlabel, ylabel, title, filename):
    plt.scatter(X1, Y1, marker='x', label = label1)
    plt.scatter(X2, Y2, marker='x', label = label2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Start axes from origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    # Annotete points
    for i in range(len(linear_n)):
        plt.annotate(linear_k[i], (linear_n[i] - 2e7, 0), color="red")
    plt.annotate("k =", (1e7, 0), color="red")

    plt.legend()
    plt.title(title)

    plt.savefig(filename)
    plt.show()

linear_n, linear_time, linear_mem, linear_k = parse("linear_out.txt")
basic_n, basic_time, basic_mem, basic_k = parse("basic_out.txt")

do_scatter_plot(linear_n, linear_time, basic_n, basic_time, "Linear", "Basic", "SBWT columns", "Time (s)", "Coli3682", "time.png")
do_scatter_plot(linear_n, linear_mem, basic_n, basic_mem, "Linear", "Basic", "SBWT columns", "Memory (GB)", "Coli3682", "mem.png")
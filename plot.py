import sys
import matplotlib.pyplot as plt
import os

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
    plt.figure()
    plt.scatter(X1, Y1, marker='x', label = label1)
    plt.scatter(X2, Y2, marker='x', label = label2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Start axes from origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    # Annotate points
    #for i in range(len(linear_n)):
    #    plt.annotate(linear_k[i], (linear_n[i] - 2e7, 0), color="red")
    #plt.annotate("k =", (1e7, 0), color="red")

    plt.legend()
    plt.title(title)

    print("Saving to", filename)
    plt.savefig(filename)

def plot_human(linear_infile, basic_infile, dataset_name, outfile_prefix):
    linear_n, linear_time, linear_mem, linear_k = parse(linear_infile)
    basic_n, basic_time, basic_mem, basic_k = parse(basic_infile)

    do_scatter_plot(linear_n, linear_time, basic_n, basic_time, 
                    "Linear", "Basic", "SBWT columns", "Time (s)", 
                    "Time ({})".format(dataset_name), 
                    "{}_time_by_n.pdf".format(outfile_prefix))
    
    do_scatter_plot(linear_n, linear_mem, basic_n, basic_mem, 
                    "Linear", "Basic", "SBWT columns", "Memory (GB)", 
                    "Memory ({})".format(dataset_name), 
                    "{}_mem_by_n.pdf".format(outfile_prefix))

    do_scatter_plot(linear_k, linear_time, basic_k, basic_time, 
                    "Linear", "Basic", "k", "Time (s)", 
                    "Time ({})".format(dataset_name), 
                    "{}_time_by_k.pdf".format(outfile_prefix))
    
    do_scatter_plot(linear_k, linear_mem, basic_k, basic_mem, 
                    "Linear", "Basic", "k", "Memory (GB)", 
                    "Memory ({})".format(dataset_name), 
                    "{}_mem_by_k.pdf".format(outfile_prefix))
    
if not os.path.exists("plots"):
   os.makedirs("plots")

plot_human("data_for_plots/linear_human.txt", "data_for_plots/basic_human.txt", "Human genome", "plots/human")
plot_human("data_for_plots/linear_coli.txt", "data_for_plots/basic_coli.txt", "E. coli genomes", "plots/coli")
plot_human("data_for_plots/linear_metagenome.txt", "data_for_plots/basic_metagenome.txt", "Metagenome reads", "plots/metagenome")
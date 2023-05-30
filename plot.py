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
    
    do_scatter_plot(basic_time, basic_mem, linear_time, linear_mem, 
                    "Linear", "Basic", "Memory (GB)", "Time (s)", 
                    "Space and time ({})".format(dataset_name), 
                    "{}_mem_time_tradeoff.pdf".format(outfile_prefix))

def read_kmer_curve(filename):
    K, SIZE = [], [] 
    for line in open(filename):
        tokens = line.split(",")
        k, size = int(tokens[-2]), int(tokens[-1])
        K.append(k)
        SIZE.append(size)
    return K, SIZE

def make_kmer_plot(human_file, coli_file, metagenome_file, outfile):

    human_K, human_SIZE = read_kmer_curve(human_file)
    coli_K, coli_SIZE = read_kmer_curve(coli_file)
    metagenome_K, metagenome_SIZE = read_kmer_curve(metagenome_file)

    plt.figure()

    plt.plot(human_K, human_SIZE, label="Human genome", marker="x")
    plt.plot(coli_K, coli_SIZE, label="E. coli genomes", marker="x")
    plt.plot(metagenome_K, metagenome_SIZE, label="Metagenome reads", marker="x")

    plt.xlabel("k")
    plt.ylabel("n")

    # Start axes from origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.legend()
    plt.title("Number of sets in SBWT")

    print("Saving to", outfile)
    plt.savefig(outfile)


if not os.path.exists("plots"):
   os.makedirs("plots")

make_kmer_plot("data_for_plots/linear_human.csv", "data_for_plots/linear_coli.csv", "data_for_plots/linear_metagenome.csv", "plots/kmers.pdf")
plot_human("data_for_plots/linear_human.csv", "data_for_plots/basic_human.csv", "Human genome", "plots/human")
plot_human("data_for_plots/linear_coli.csv", "data_for_plots/basic_coli.csv", "E. coli genomes", "plots/coli")
plot_human("data_for_plots/linear_metagenome.csv", "data_for_plots/basic_metagenome.csv", "Metagenome reads", "plots/metagenome")
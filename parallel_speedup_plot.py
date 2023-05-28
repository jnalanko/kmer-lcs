import sys
import matplotlib.pyplot as plt

def parse(filename):
    times = []
    for line in open(filename):
        tokens = line.split(',')
        times.append(float(tokens[1]))
    return times

times = parse(sys.argv[1])
threads = list(range(1, len(times) + 1))

speedups = [times[0] / t for t in times]
plt.plot(threads, speedups)

# Plot dashed line at y = x
plt.plot(threads, threads, linestyle='--', color='black')
plt.xlabel("Threads")
plt.ylabel("Speedup")

plt.show()

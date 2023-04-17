import matplotlib.pyplot as plt
import pandas as pd
import sys
import subprocess

path = "./builddir/parallel-sorting"

scaling = [1, 2, 4, 8]
N = [80000000, 40000000, 20000000, 10000000]
# N = [100, 100, 100, 100]
repeat = 3

measurements = []

for benchmark in range(len(scaling)):
    processes = scaling[benchmark]
    items = N[benchmark]

    cmd = "mpirun -n {} {} 0 {} 2 {}".format(processes, path, items, items)
    print(cmd)
    runtime_avg = 0
    for run in range(repeat):
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        runtime = float(result.stdout.split("Elapsed time: ")[1])
        runtime_avg += runtime
    runtime_avg /= repeat
    print("Runtime: {}".format(runtime_avg))
    measurements.append([processes, runtime_avg, 1])

for benchmark in range(1, len(scaling)):
    # Calculate speedup
    measurements[benchmark][2] = (
        measurements[benchmark - 1][1] / measurements[benchmark][1]
    )

speedup_df = pd.DataFrame(measurements, columns=["# of Processes", "Time", "Speedup"])

print(speedup_df)
speedup_df.to_pickle("speedup.pkl")

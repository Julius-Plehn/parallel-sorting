import sys
import random
import sorting
import timeit
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    random.seed(10)
    rs = [5, 7, 9, 11]
    sizes = np.linspace(0, 200, 200, endpoint=False, dtype=int)
    # Large test
    # sizes = np.linspace(0, 1000, 1000, endpoint=False, dtype=int)

    sys.setrecursionlimit(sizes[len(sizes) - 1])

    repeat = 5

    timings_quicksort = []
    timings_sorted = []
    timings_median = {}
    for r in rs:
        timings_median[r] = []

    for size in sizes:
        input = list(range(1, size))
        random.shuffle(input)

        timing_sorted = timeit.timeit(
            stmt="sorting.insertion_sort(input.copy(), 0, len(input))",
            globals=globals(),
            number=repeat,
        )
        timings_sorted.append(timing_sorted / repeat)
        print("Time: {}".format(timing_sorted / repeat))

        # Quicksort using last element as pivot

        timing_quicksort = timeit.timeit(
            stmt="sorting.quicksort(input.copy(), 0, len(input)-1)",
            globals=globals(),
            number=repeat,
        )
        timings_quicksort.append(timing_quicksort / repeat)
        print("Time: {}".format(timing_quicksort / repeat))

        for r in rs:
            print("Size: {}, r={}".format(size, r))
            # sorted_input = sorting.quicksort_median(input, 0, len(input)-1, r=r)
            timing_median = timeit.timeit(
                stmt="sorting.quicksort_median(input.copy(), 0, len(input)-1, r=r, n=50)",
                globals=globals(),
                number=repeat,
            )
            timings_median[r].append(timing_median / repeat)
            print("Time: {}".format(timing_median / repeat))

    fig, ax = plt.subplots()
    plt.xlabel("Size n")
    plt.ylabel("Time in Sec.")
    fig.set_size_inches(14, 5)

    # ax.set_xticks(sizes)
    ax.plot(
        sizes,
        timings_sorted,
        marker="o",
        markersize=1,
        # linewidth=2.0,
        label="Insertion",
    )
    ax.legend()

    ax.plot(
        sizes,
        timings_quicksort,
        marker="o",
        markersize=1,
        # linewidth=2.0,
        label="Quicksort (Random Pivot)",
    )
    ax.legend()

    for r in rs:
        ax.plot(
            sizes,
            timings_median[r],
            marker="o",
            markersize=1,
            # linewidth=2.0,
            label="MoM r: {}".format(r),
        )
        ax.legend()

    plt.savefig("benchmark.pdf", dpi=200)

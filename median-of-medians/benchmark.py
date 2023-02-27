import sys
import random
import sorting
import timeit
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # r = int(sys.argv[1])
    # size = int(sys.argv[2])
    random.seed(10)
    rs = [5, 7, 9, 11]
    sizes = [100, 1000, 10000, 100000]
    # sizes = [100, 1000, 10000]
    sys.setrecursionlimit(sizes[len(sizes) - 1])

    repeat = 1

    timings_quicksort = []
    timings_sorted = []
    timings_median = {}
    for r in rs:
        timings_median[r] = []

    for size in sizes:
        input = list(range(1, size))
        random.shuffle(input)

        timing_sorted = timeit.timeit(
            stmt="sorting.insertion_sort(input.copy())",
            globals=globals(),
            number=repeat,
        )
        timings_sorted.append(timing_sorted / repeat)
        print("Time: {}".format(timing_sorted / repeat))

        # Quicksort using last element as pivot
        """
        timing_quicksort = timeit.timeit(stmt="sorting.quicksort(input, 0, len(input)-1)", globals=globals(), number=repeat)
        timings_quicksort.append(timing_quicksort/repeat)
        print("Time: {}".format(timing_quicksort/repeat))
        """
        for r in rs:
            print("Size: {}, r={}".format(size, r))
            # sorted_input = sorting.quicksort_median(input, 0, len(input)-1, r=r)
            timing_median = timeit.timeit(
                stmt="sorting.quicksort_median(input.copy(), 0, len(input)-1, r=r)",
                globals=globals(),
                number=repeat,
            )
            timings_median[r].append(timing_median / repeat)
            print("Time: {}".format(timing_median / repeat))

    print(timings_quicksort)
    print(timings_sorted)
    print(timings_median)

    fig, ax = plt.subplots()
    plt.xlabel("Size")
    plt.ylabel("Time in Sec.")
    fig.set_size_inches(14, 5)

    ax.set_xticks(sizes)
    ax.plot(
        sizes,
        timings_sorted,
        marker="o",
        markersize=3,
        # linewidth=2.0,
        label="Insertion",
    )
    ax.legend()

    """
    ax.plot(
        sizes,
        timings_quicksort,
        marker="o",
        markersize=3,
        # linewidth=2.0,
        label="Quicksort",
    )
    ax.legend()
    """
    for r in rs:
        ax.plot(
            sizes,
            timings_median[r],
            marker="o",
            markersize=3,
            # linewidth=2.0,
            label="MoM r: {}".format(r),
        )
        ax.legend()

    plt.savefig("benchmark.png", dpi=200)

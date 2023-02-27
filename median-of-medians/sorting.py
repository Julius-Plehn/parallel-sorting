import sys
import math
import random


def compare(A, B):
    if sorted(A) == B:
        print("Correctly sorted")
    else:
        print("Error")


# Earlier experiment, not used
def partition_old(array, low, high, median=None):
    if median is None:
        pivot_element = array[high]
        pivot_index = high
    else:
        pivot_element = median
        pivot_index = array.index(median)

    array[pivot_index] = array[high]
    array[high] = pivot_element
    pivot_index = high

    idx = low
    moved_high = 1
    for _ in array[low:high]:
        if idx == pivot_index:
            idx += moved_high
        element = array[idx]
        if element > pivot_element:
            move = array.pop(array.index(element))
            array.insert(high, move)
            moved_high += 1
            pivot_index -= 1
        else:
            idx += 1

    return (array, pivot_index)


def partition(array, low, high, median=None):
    # Random pivot
    if median is None:
        median_idx = random.randint(low, high)
        median = array[median_idx]

    for i in range(low, high + 1):
        if array[i] == median:
            array[i], array[high] = array[high], array[i]
            break

    i = low - 1
    for j in range(low, high):
        if array[j] <= array[high]:
            i += 1
            array[i], array[j] = array[j], array[i]
    array[i + 1], array[high] = array[high], array[i + 1]
    return array, i + 1


def quicksort(array, low, high):
    if low < high:
        array, partitioning_index = partition(array, low, high)
        quicksort(array, low, partitioning_index - 1)
        quicksort(array, partitioning_index + 1, high)
    return array


def median_of_medians(array, low, high, k, r=5):
    median_idx = r // 2
    medians = []
    block_range = high - low

    if block_range > 0:
        for b in range(low, high + 1, r):
            end_idx = min(b + r - 1, high)
            if end_idx == high:
                median_idx = (end_idx - b) // 2
            medians.append(sorted(array[b : end_idx + 1])[median_idx])
    else:
        medians.append(low)
    MoM = 0
    if len(medians) == 1:
        return medians[0]
    else:
        MoM = median_of_medians(medians, 0, len(medians) - 1, len(medians) // 2, r)

    _, q = partition(array, low, high, MoM)
    i = q - low + 1
    if i == k:
        return array[q]
    elif i > k:
        return median_of_medians(array, low, q - 1, k, r)
    else:
        return median_of_medians(array, q + 1, high, k - i, r)


def quicksort_median(array, low, high, r=5, n=100):
    if low < high:
        if high - low < n:
            return insertion_sort(array, low, high + 1)
        m = median_of_medians(array, low, high, (high - low) // 2, r=r)
        array, p = partition(array, low, high, m)
        quicksort_median(array, low, p - 1, r=r, n=n)
        quicksort_median(array, p + 1, high, r=r, n=n)
    return array


def insertion_sort(array, low, high):
    for key in range(low, high):
        item = array[key]
        cmp = key - 1
        while cmp >= 0 and item < array[cmp]:
            array[cmp + 1] = array[cmp]
            cmp -= 1
        array[cmp + 1] = item
    return array


if __name__ == "__main__":
    # sys.setrecursionlimit(10)
    if len(sys.argv) < 2:
        print("Provide a file as input")
        exit()
    r = 5
    if len(sys.argv) == 3:
        r = int(sys.argv[2])
    f = open(sys.argv[1], "r")
    input = list(map(int, f.read().split()))

    sorted_input = quicksort_median(input, 0, len(input) - 1, r=r)
    print(sorted_input)

    # Various Tests below
    """
    # Check:
    compare(input, sorted_input)

    input = list(range(1, 1000))
    random.shuffle(input)

    sorted_input = quicksort_median(input, 0, len(input) - 1, r=r)
    compare(input, sorted_input)

    input = list(range(1, 100))
    random.shuffle(input)
    # sorted_input = quicksort(input, 0, len(input) - 1)
    # compare(input, sorted_input)

    test = [9, 5, 1, 4, 3]
    sorted_test = insertion_sort(test, 0, len(test))

    compare(test, sorted_test)
    """

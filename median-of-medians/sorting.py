import sys
import math
import random


def compare(A, B):
    if sorted(A) == B:
        print("Correctly sorted")
    else:
        print("Error")


# Version with last element as pivot element
def partition(array, low, high):
    pivot_element = array[high]
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


def quicksort(array, low, high):
    if low < high:
        array, partitioning_index = partition(array, low, high)
        quicksort(array, low, partitioning_index - 1)
        quicksort(array, partitioning_index + 1, high)
    return array


def median_of_medians(array, low, high, r=5):
    median_idx = math.ceil(r / 2) - 1
    medians = []
    block_range = high - low
    blocks = math.ceil(block_range / r)
    for b in range(blocks):
        start_idx = low + b * r
        end_idx = start_idx + r - 1
        if b == blocks - 1:
            end_idx = high
            median_idx = math.ceil((end_idx - start_idx) / 2) - 1
        medians.append(sorted(array[start_idx : end_idx + 1])[median_idx])

    pivot_element = sorted(medians)[math.ceil(len(medians) / 2) - 1]
    pivot_index = low + array[low : high + 1].index(pivot_element)

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


def quicksort_median(array, low, high, r=5):
    if low < high:
        array, partitioning_index = median_of_medians(array, low, high, r=r)
        quicksort_median(array, low, partitioning_index - 1, r=r)
        quicksort_median(array, partitioning_index + 1, high, r=r)
    return array


def insertion_sort(array):
    for key in range(1, len(array)):
        item = array[key]
        cmp = key - 1
        while cmp >= 0 and item < array[cmp]:
            array[cmp + 1] = array[cmp]
            cmp -= 1
        array[cmp + 1] = item
    return array


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Provide a file as input")
        exit()
    r = 5
    if len(sys.argv) == 3:
        r = int(sys.argv[2])
    f = open(sys.argv[1], "r")
    input = list(map(int, f.read().split()))
    # print("Input:")
    print(input)
    sorted_input = quicksort_median(input, 0, len(input) - 1, r=r)

    # Check:
    compare(input, sorted_input)

    input = list(range(1, 1000))
    random.shuffle(input)

    sorted_input = quicksort_median(input, 0, len(input) - 1, r=r)

    compare(input, sorted_input)

    input = list(range(1, 100))
    random.shuffle(input)
    sorted_input = quicksort(input, 0, len(input) - 1)
    compare(input, sorted_input)

    test = [9, 5, 1, 4, 3]
    sorted_test = insertion_sort(test)

    compare(test, sorted_test)

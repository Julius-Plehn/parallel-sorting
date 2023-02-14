import sys

def partition(array, low, high):
    print("Partition [{}][{}]".format(low, high))
    # Version with last element as pivot element
    pivot_element = array[high]
    pivot_index = high
    idx = 0
    for _ in array[:high]:
        #print(idx)
        element = array[idx]
        if element > pivot_element:
            move = array.pop(array.index(element))
            array.insert(high,move)
            pivot_index -= 1
        else:
            idx +=1
        #print(array)
    print(array)
    #exit()
    return (array, pivot_index)

def quicksort(array, low, high, partition_f):
    if low < high:
        array, partitioning_index = partition_f(array, low, high)
        #exit()
        quicksort(array, low, partitioning_index-1, partition_f)
        # maybe start with +1 in the following
        quicksort(array, partitioning_index+1, high, partition_f)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Provide a file as input")
        exit()
    f = open(sys.argv[1], "r")
    input = list(map(int, f.read().split()))
    print("Input:")
    print(input)
    quicksort(input, 0, len(input)-1, partition)


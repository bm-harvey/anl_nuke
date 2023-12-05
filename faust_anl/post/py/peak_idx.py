import random
import matplotlib.pyplot as plt


def main():
    peaks_1 = [random.uniform(20, 150) for _ in range(11)]
    peaks_2 = [random.uniform(20, 150) for _ in range(4)]
    peaks_3 = [random.uniform(20, 150) for _ in range(3)]
    peaks_4 = [random.uniform(20, 150) for _ in range(2)]

    peaks_1.sort()
    peaks_2.sort()
    peaks_3.sort()
    peaks_4.sort()

    indices_1 = [i for i in range(11)]
    indices_2 = [find_index(peaks_1, val) for val in peaks_2]
    indices_3 = [find_index(peaks_1, val) for val in peaks_3]
    indices_4 = [find_index(peaks_1, val) for val in peaks_4]

    plt.scatter(x=indices_1, y=peaks_1)
    plt.scatter(x=indices_2, y=peaks_2)
    plt.scatter(x=indices_3, y=peaks_3)
    plt.scatter(x=indices_4, y=peaks_4)

    plt.savefig("peak_idx.png")


def find_index(ref_data: list[float], value) -> int:
    min_idx = 0
    min_val = abs(ref_data[0] - value)

    for idx in range(1, len(ref_data)):
        val = abs(ref_data[idx] - value)
        if val < min_val:
            min_val = val
            min_idx = idx
    return min_idx


if __name__ == "__main__":
    main()

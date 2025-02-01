import numpy as np

def generate_matrix(nx, ny, total_sum):
    """
    Generate a matrix of the specified size (nx, ny) randomly,
    and adjust the sum of the elements to total_sum.
    """
    matrix = np.random.rand(ny, nx)  # generate random matrix
    matrix *= total_sum / matrix.sum()  # scale the matrix
    return matrix

def main():
    nx, ny = 80, 80
    total_sum = 0.009
    matrix = generate_matrix(nx, ny, total_sum)

    for row in matrix:
        print(" ".join(f"{val:.8f}" for val in row))

if __name__ == "__main__":
    main()

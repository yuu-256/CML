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

    nx, ny = 80, 40
    w_sum = 0.009
    w_matrix = generate_matrix(nx, int(ny*2), w_sum)
    w_l = w_matrix[:, :40]
    w_v = w_matrix[:, 40:]

    E_sum = 4.0*nx*ny
    E_matrix = generate_matrix(nx, ny, E_sum)
    
    u_sum = 0.5*nx*ny
    v_sum = 0.5*nx*ny
    u_matrix = generate_matrix(nx, ny, u_sum)
    v_matrix = generate_matrix(nx, ny, u_sum)
    # マイナスの値を持たせる
    u_matrix = u_matrix - 0.5 
    v_matrix = v_matrix - 0.5

    fnm = f"config/config_{nx}_{ny}_{w_sum}.txt"
    with open(fnm, "w") as f:
        f.write("# Simulation Configuration File\n")
        f.write(f"# Grid Size: {nx}x{ny}\n")
        f.write("# Energy E(x,y)\n")
        f.write("E:\n")
        np.savetxt(f, E_matrix, fmt="%.8f")
        f.write("# velocity u(x,y)\n")
        f.write("u:\n")
        np.savetxt(f, u_matrix, fmt="%.8f")
        f.write("# velocity v(x,y)\n")
        f.write("v:\n")
        np.savetxt(f, v_matrix, fmt="%.8f")
        f.write("# Liquid water content w_l(x, y)\n")
        f.write("w_l:\n")
        np.savetxt(f, w_l, fmt="%.8f")
        f.write("# Ice water content w_v(x, y)\n")
        f.write("w_v:\n")
        np.savetxt(f, w_v, fmt="%.8f")

if __name__ == "__main__":
    main()

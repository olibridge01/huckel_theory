import numpy as np
import scipy
import networkx as nwx
import sys

def h_linear(n: int):
    """
    Returns the Hückel Hamiltonian matrix for a linear polyene of length n.
    
    Input: number of atoms
    Output: np.ndarray Hamiltonian matrix for the molecule
    """

    if n < 1:
        raise ValueError("n must be at least 1!")
    H = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if j == i - 1 or j == i + 1:
                H[i,j] = -1
    return H


def h_cyclic(n: int):
    """
    Returns the Hückel Hamiltonian matrix for a cyclic polyene of length n.
    
    Input: number of atoms
    Output: np.ndarray Hamiltonian matrix for the molecule
    """

    if n < 3:
        raise ValueError("n must be at least 1!")
    H = h_linear(n)
    
    H[n-1, 0] = -1
    H[0, n-1] = -1
    
    return H


def platonic(n):
    """
    Returns the Hückel Hamiltonian matrix for an sp2 hybridised platonic solid.
    Uses Python library NetworkX, a useful library for graph theory and network analysis.
    NetworkX can generate graph objects for the platonic solids, and the Hamiltonian matrix is simply -1 * adjacency matrix for the graph.
    
    Input: n = 4, 6, 8, 12, 20
    Output: ndarray Hamiltonian matrix
    """

    graph = nwx.Graph()
    
    if n == 4:
        graph = nwx.tetrahedral_graph()
    elif n == 6:
        graph = nwx.octahedral_graph()
    elif n == 8: 
        graph = nwx.cubical_graph()
    elif n == 12:
        graph = nwx.icosahedral_graph()
    elif n == 20:
        graph = nwx.dodecahedral_graph()
    else:
        raise ValueError("Error: n must be 4, 6, 8, 12 or 20")
        
    H = -1 * nwx.to_numpy_array(graph)
    
    return H


def get_evals(n):
    """
    Takes a matrix as a np.ndarray and returns a sorted np.ndarray of eigenvalues
    """

    evals = np.round(np.linalg.eigvals(n),3)
    return sorted(evals)


def degeneracies(evals: np.ndarray):
    """
    Take a sorted list of eigenvalues and produce a list of the eigenvalues with their corresponding degeneracies.
    
    Input: Numpy array of eigenvalues
    Output: Paired list of eigenvalues and degeneracies
    """

    counter = 1
    marker = ""
    degeneracies = []
    
    # Iterate through eigenvalues and check for repeated values, producing a paired list of e_val and degeneracy.
    for e_val in evals:
        if e_val == marker:
            counter += 1
        else:
            if marker != "":
                degeneracies.append((marker, counter))
            marker = e_val
            counter = 1
    
    if counter > 0:
        degeneracies.append([marker, counter])
        
    return degeneracies


def print_energies(degen_list):
    """
    Takes a paired list of eigenvalues and degeneracies and prints them in a tabular format.
    """

    # Counter for total number of orbitals
    count = 0
    
    # Printing the table
    print()
    print("{:<10} {:<10}".format(" Energy","Degeneracy"))
    print("---------------------")
    
    # Printing energies/degeneracies from highest to lowest energy
    for i in reversed(degen_list):
        energy, degeneracy = i
        count += degeneracy
        if energy < 0:
            print("-{:<10}    {:<10}".format(f'{abs(energy):.3f}',degeneracy))
        else:
            print(" {:<10}    {:<10}".format(f'{abs(energy):.3f}',degeneracy))
    
    print("---------------------")
    print(f"System has {count} orbitals")
    print()


def user_help():
    """
    Print a help guide for the user to know what to input into the command line.
    """
    print()
    print("User input: python3 huckel.py [linear/cyclic/platonic] [number of atoms]\n"
         "NB: for platonic solids n = 4, 6, 8, 12, 20")
    print()


# Code to call the different functions based on a command line input

platonic_nums = ['4','6','8','12','20']

if len(sys.argv) != 3 or  not sys.argv[2].isdigit():
    user_help()
else:
    if sys.argv[1] == "linear":
        H = h_linear(int(sys.argv[2]))
        print_energies(degeneracies(get_evals(H)))
    
    elif sys.argv[1] == "cyclic":
        H = h_cyclic(int(sys.argv[2]))
        print_energies(degeneracies(get_evals(H)))

    elif sys.argv[1] == "platonic" and sys.argv[2] in platonic_nums:
        H = platonic(int(sys.argv[2]))
        print_energies(degeneracies(get_evals(H)))
    
    else:
        user_help()
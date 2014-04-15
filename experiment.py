from qutip import *
from scipy import *

from multiprocessing import Pool
from random import randrange
import sys

###############################################################################
# Calculate the minimum eigenvalue gap
# The creation of this function is documented through the IPython notebooks
#
# Parameters:
#   h_b - the base Hamiltonian
#   h_p - the problem Hamiltonian
#
def calculate_min_gap(h_b, h_p):
  # Get the number of qubits
  assert len(h_b.dims[0]) == len(h_b.dims[1]) == len(h_p.dims[0]) == len(h_p.dims[1])
  n = len(h_b.dims[0])
  
  # Increase taumax to make the sweep more adiabatic
  taumax = 5.0
  taulist = linspace(0, taumax, 100)
  
  # The time dependent function
  h_t = [[h_b, lambda t, t_max : (t_max-t)/t_max],
          [h_p, lambda t, t_max : t/t_max]]
  psi0 = tensor([basis(2,0) for _ in range(n)])
  
  # Must be an array because callbacks need pass by reference
  min_gap = [+inf]

  # The callback function
  def process_rho(tau, psi):
    H = qobj_list_evaluate(h_t, tau, taumax)
    # eigvals=2 returns just the two smallest eigenvalues
    evals, ekets = H.eigenstates(eigvals=2)
            
    min_gap[0] = min(min_gap[0], evals[1]-evals[0])
  
  mesolve(h_t, psi0, taulist, [], process_rho, taumax)
          
  return min_gap[0]

###############################################################################
# Create a base Hamiltonian whose ground state is known
# This follows the general formula found in Farhi et al.
#
# Parameters:
#   dims - the number of qubits to use
#
def base(dims):
  si = qeye(2)
  sx = sigmax()
  sx_list = []
  
  # Create all the components to make up this basis
  for n in range(dims):
    op_list = []
    for m in range(dims):
      op_list.append(si)

    op_list[n] = sx
    sx_list.append(tensor(op_list))
  
  # Sum them together, after converting from {0,1} to {-1,1}
  h_b = 0
  for n in range(dims):
    h_b += 0.5 * (1 - sx_list[n])
  
  return h_b

###############################################################################
# Convert a QUBO problem to an Ising spin glass Hamiltonian
#
# Parameters:
#   J - 2D Array of connection biases
#   h - 1D Array of individual biases
#
def convert_ising(J, h):
  assert(len(J) == len(h))
  
  n = len(J)
  
  si = qeye(2)
  sz = sigmaz()
  
  h_ising = 0
  
  # Convert the J matrix
  for i in range(0, n-1):
    for j in range(i+1, n):
      op_list = []
      
      for _ in range(n):
        op_list.append(si)
      
      op_list[i] = 0.5 * (1 - sz)
      op_list[j] = 0.5 * (1 - sz)
      
      h_ising += J[i][j] * tensor(op_list)
  
  # Convert the h matrix
  for i in range(0, n):
    op_list = []
        
    for _ in range(n):
      op_list.append(si)
    
    op_list[i] = 0.5 * (1 - sz)
    
    h_ising += h[i] * tensor(op_list)
    
  return h_ising

###############################################################################
# Create a random QUBO for a particular number of qubits
# then calculate its corresponding minimum gap
#
# Parameters:
#   n - the number of qubits to use
#
def random_sim(n):
  J = numpy.zeros(n*n).reshape((n, n))
  h = numpy.zeros(n)
  
  # Connection biaseshave a range of -100 to 100
  for i in range(0, n-1):
    for j in range(i+1, n):
      J[i][j] = randrange(-100,100)
  
  # Individual biases have a range of -100 to 100
  for i in range(0, n):
    h[i] = randrange(-100,100)
    
  return calculate_min_gap(base(n), convert_ising(J, h))

# Some parameters for the experiment
num_trials = 1000
base_dir = "data/correct_experiment/" # with trailing slash

###############################################################################
# Run num_trials experiments for a particular number of qubits
# then write the minimum values to a file in base_dir/num_qubits.out
#
# Parameters:
#   num_qubits - the number of qubits to use
#
def output_result(num_qubits):
  count = 0

  while count < num_trials:
    res = random_sim(num_qubits)
    
    with open(base_dir + str(num_qubits) + ".out", "a") as outfile:
      outfile.write(str(res) + "\n")

    count += 1

    sys.stderr.write(str(num_qubits) + ": Completed = " + str(count) + "\n")

# Use multiple processors if possible
if __name__ == '__main__':
  pool = Pool(processes=8)
  pool.map(output_result, range(1,13))
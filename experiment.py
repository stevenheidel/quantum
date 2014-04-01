from qutip import *
from scipy import *
from random import randrange

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
  
  min_gap = [+inf]
  def process_rho(tau, psi):
    H = qobj_list_evaluate(h_t, tau, taumax)
    evals, ekets = H.eigenstates(eigvals=2)
            
    min_gap[0] = min(min_gap[0], evals[1]-evals[0])
  
  mesolve(h_t, psi0, taulist, [], process_rho, taumax)
          
  return min_gap[0]


def base(dims):
  si = qeye(2)
  sx = sigmax()
  sx_list = []
  
  for n in range(dims):
    op_list = []
    for m in range(dims):
      op_list.append(si)

    op_list[n] = sx
    sx_list.append(tensor(op_list))
  
  h_b = 0
  for n in range(dims):
    h_b += 0.5 * (1 - sx_list[n])
  
  return h_b

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
  
  # Scale so that zero is minimum
  #h_ising = h_ising - h_ising.diag().min()
  
  return h_ising

def random_sim(n):
  J = numpy.zeros(n*n).reshape((n, n))
  h = numpy.zeros(n)
  
  for i in range(0, n-1):
    for j in range(i+1, n):
      J[i][j] = randrange(100)
  
  for i in range(0, n):
    h[i] = randrange(100)
  
  #print J, h
  
  return calculate_min_gap(base(n), convert_ising(J, h))

num_qubits = int(sys.argv[1])
num_trials = int(sys.argv[2])

count = 0

while count < num_trials:
  res = random_sim(num_qubits)
  
  # Get rid of boring results
  if res > 0.0 + 1e-9 and res < 1.0 - 1e-9:
    print res
    count += 1
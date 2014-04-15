from pylab import *

# Some parameters
base_dir = "data/correct_experiment/" # with trailing slash
n = 12 # max qubits to use
one = 1 # min qubits to use
epsilon = 1e-9

# Arrays needed for later
X = []
Y = []
medians = []
means = []

# Return true if number is not within epsilon of 1
def filter_ones(x):
  return x < 1.0 - epsilon

# Return true if number is not within epsilon of 0
def filter_zeros(x):
  return x > 0.0 + epsilon

# Convert to time
def t(list):
  return map(lambda x: x**-2, list)

# Log a list of values
def l(list):
  return map(lambda x: log(x), t(list))

# For each number of quibts, get all the data and process it
for i in range(one,n+1):
  with open(base_dir + str(i) + ".out") as f:
    floats = map(float, f)

    # Get rid of ones and/or zeros
    floats = filter(filter_zeros, floats)
    #floats = filter(filter_ones, floats)

    # Only use a subset of the total samples
    floats = floats[:128]

    print str(len(floats)) + " samples of " + str(i)

    # Calculate the medians
    medians.append(median(floats))
    means.append(mean(floats))

    # Store the scatter
    X.extend(tile(i, len(floats)))
    Y.extend(floats)

# Find a degree 2 poly best fit line
p = poly1d(polyfit(X,Y,2))
xp = linspace(1,n,1000)

############################################################################
# CHART 1: Minimum Eigenvalue Gaps

# Plot the scatter
p1, = plot(X, Y, '.')

# Plot the best fit line
p2, = plot(xp, p(xp), '-')

# Plot the medians and means
p3, = plot(range(one,n+1), medians)
p4, = plot(range(one,n+1), means)

# Now some things to make chart look nice
xlim(0, n+1) # add padding left and right
ylim(0, 1)
xticks(range(1, n+1)) # show all qubit markers
xlabel("# of Qubits")
ylabel("Minimum Eigenvalue Gap\n(proportional to original eigenvalue gap)", multialignment='center')

# Legend
legend([p1,p2,p3,p4], ["Samples", "Best fit", "Medians", "Means"], loc=3)

# Plot
show()

############################################################################
# CHART 2: Time Required according to Adiabatic Theorem

# Plot the best fit line
p1, = plot(xp, p(xp)**-2, '-')

# Plot the medians and means
p2, = plot(range(one, n+1), t(medians))
p3, = plot(range(one, n+1), t(means))

# Now some things to make chart look nice
xlim(0, n+1) # add padding left and right
xticks(range(1, n+1)) # show all qubit markers
xlabel("# of Qubits")
ylabel("Time Required")

# Legend
legend([p1,p2,p3], ["Best fit", "Medians", "Means"], loc=2)

# Plot
show()

############################################################################
# CHART 3: Log Scatter Plot of Time Required

# Plot the scatter
p1, = plot(X, l(Y), '.')

# Plot the medians and means
p2, = plot(range(one, n+1), l(medians))
p3, = plot(range(one, n+1), l(means))

# Now some things to make chart look nice
xlim(0, n+1) # add padding left and right
xticks(range(1, n+1)) # show all qubit markers
xlabel("# of Qubits")
ylabel("Time Required (log)")

# Legend
legend([p1,p2,p3], ["Samples", "Medians", "Means"], loc=2)

# Plot
show()
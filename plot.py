from pylab import *

n = 10
X = []
Y = []

def filter_ones(x):
  return x < 1.0 - 1e-9
def filter_zeros(x):
  return x > 0.0 + 1e-9

for i in range(1,n+1):
  with open("data/" + str(i) + ".out") as f:
    floats = map(float, f)

    # Get rid of ones and zeros
    floats = filter(filter_zeros, floats)
    floats = filter(filter_ones, floats)

    print str(len(floats)) + " samples of " + str(i)

    X.extend(tile(i, len(floats)))

    Y.extend(floats)

p = poly1d(polyfit(X,Y,3))

xp = linspace(1,n,1000)

plot(X,Y,'.', xp,p(xp),'-')
xlim(0,n+1)
xticks(range(1,n+1))
#ylim(0,1)
show()
from pylab import *

base_dir = "data/integers_negpos100/" # with trailing slash
n = 12
one = 1
X = []
Y = []
medians = []
means = []

def filter_ones(x):
  return x < 1.0 - 1e-9
def filter_zeros(x):
  return x > 0.0 + 1e-9

for i in range(one,n+1):
  with open(base_dir + str(i) + ".out") as f:
    floats = map(float, f)

    # Get rid of ones and zeros
    #floats = filter(filter_zeros, floats)
    #floats = filter(filter_ones, floats)

    floats = floats[:150]
    #floats = floats[:1000]

    print str(len(floats)) + " samples of " + str(i)

    medians.append(median(floats))
    means.append(mean(floats))

    X.extend(tile(i, len(floats)))

    Y.extend(floats)

p = poly1d(polyfit(X,Y,2))

xp = linspace(1,n,1000)

plot(X,Y,'.', xp,p(xp),'-')
plot(range(one,n+1),medians)
plot(range(one,n+1),means)
xlim(0,n+1)
xticks(range(1,n+1))

#plot(xp, p(xp)**-2, '-')
show()
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
    floats = filter(filter_zeros, floats)
    #floats = filter(filter_ones, floats)

    floats = floats[:130]
    #floats = floats[:1000]

    print str(len(floats)) + " samples of " + str(i)

    medians.append(median(floats))
    means.append(mean(floats))

    X.extend(tile(i, len(floats)))

    Y.extend(floats)

p = poly1d(polyfit(X,Y,2))

xp = linspace(1,n,1000)

xlim(0,n+1)

ylim(0,1)
plot(X,Y,'.', xp,p(xp),'-')
plot(range(one,n+1),medians)
plot(range(one,n+1),means)
xticks(range(1,n+1))
'''

# Convert to time
def t(list):
  return map(lambda x: x**-2, list)

# Log time list
def l(list):
  return map(lambda x: log(x), t(list))

ylim(0,10)
plot(X, t(Y), '.')
plot(range(one,n+1), t(medians))
plot(range(one,n+1), t(means))
#plot(xp, p(xp)**-2, '-')
'''
show()
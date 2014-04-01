from pylab import *

n = 10
X = []
Y = []

for i in range(1,n+1):
  with open("data/" + str(i) + ".out") as f:
    floats = map(float, f)

    X.extend(tile(i, len(floats)))

    Y.extend(floats)

p = poly1d(polyfit(X,Y,4))

xp = linspace(1,n,1000)

plot(X,Y,'.', xp,p(xp),'-')
xlim(0,n+1)
xticks(range(1,n+1))
#ylim(0,1)
show()
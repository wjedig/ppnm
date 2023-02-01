def f(x):
    return 1/x

s = 0.0
for i in range(1,1000):
    s+=f(i)

print(s)

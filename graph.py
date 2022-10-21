import matplotlib.pyplot as plt

# pip install matplotlib


file1 = open('results.txt', 'r')
Lines = file1.readlines()

x = []
y = []

for line in Lines:
    x1, y1, x2, y2 = line.split(" ")
    x.append(int(x1))
    y.append(int(y1))
    x.append(int(x2))
    y.append(int(y2))



plt.plot(x, y)
  
plt.xlabel('x - axis')

plt.ylabel('y - axis')

plt.title('Polygon')
  
plt.show()
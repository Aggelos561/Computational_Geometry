import matplotlib.pyplot as plt

# pip install matplotlib


file1 = open('results.txt', 'r')
Lines = file1.readlines()

x = []
y = []

for line in Lines:
    x1, y1, x2, y2 = line.split(" ")
    
    x1 = x1.replace("e+06", "")
    x1 = x1.replace(".", "")
    y1 = y1.replace("e+06", "")
    y1 = y1.replace(".", "")
    x2 = x2.replace("e+06", "")
    x2 = x2.replace(".", "")
    y2 = y2.replace("e+06", "")
    y2 = y2.replace(".", "")

    x.append(int(x1))
    y.append(int(y1))
    x.append(int(x2))
    y.append(int(y2))



plt.plot(x, y)
  
plt.xlabel('x - axis')

plt.ylabel('y - axis')

plt.title('Polygon')
  
plt.show()
import matplotlib.pyplot as plt

X = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
Y = [0, 0, 0, 1, 1, 10, 42, 305, 917, 13937]


fig = plt.figure(1)
fig.add_subplot(111)
plt.plot(X, Y, label = "Naive", color = 'g')
plt.xlabel("Длина числа")
plt.ylabel("Время (с)")
plt.tight_layout()

plt.show()
import matplotlib.pyplot as plt

generation1 = []
Fitness1 = []

f = open('SerialData.txt','r')
for row in f:
	row = row.split('    ')
	generation1.append(int(row[0]))
	Fitness1.append(float(row[1]))
generation2 = []
Fitness2 = []


f = open('ParallelData.txt','r')
for row in f:
	row = row.split('    ')
	generation2.append(int(row[0]))
	Fitness2.append(float(row[1]))


plt.plot(generation1,Fitness1 , color = 'g', label = 'Fitness of Serial Genetic Algorithm')
plt.plot(generation2,Fitness2 , color = 'r', label = 'Fitness of Prallel Genetic Algorithm')


plt.xlabel('Generation Number', fontsize = 12)
plt.ylabel('Fitness Value', fontsize = 12)

plt.title('Serial Vs Parallel Genetic Algorithm', fontsize = 20)
plt.legend()
plt.show()

import random
file = open('code_search_data.txt',mode='w')
for i in range(1000000):
	arr = []
	for j in range(3):
		arr.append(random.gauss(1,1))
	arr.append(-1 * random.gauss(100,10))
	arr.append('This is a String\n')
	file.write(' '.join([str(val) for val in arr]))
	


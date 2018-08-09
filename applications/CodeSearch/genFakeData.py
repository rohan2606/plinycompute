# import random
# file = open('code_search_data.txt',mode='w')
# for i in range(1000000):
# 	arr = []
# 	for j in range(3):
# 		arr.append(random.gauss(1,1))
# 	arr.append(-1 * random.gauss(100,10))
# 	arr.append('This is a String\n')
# 	file.write(' '.join([str(val) for val in arr]))
#

import random
import json
#file = open('code_search_data.txt',mode='w')

Programs = {}
bigGuy = []

for i in range(1000000):
        prob = []
        for j in range(3):
                prob.append(random.gauss(1,1))
        probY = -1 * random.gauss(100,10)
        Prog = 'This is a String'
        bigGuy.append({'ProbVec':prob, 'Priors':probY, 'Prog':Prog})

Programs['Programs'] = bigGuy


with open("/home/ubuntu/code_search_data.json",'w') as f:
    json.dump(Programs, f, indent=2)

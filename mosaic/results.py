import re
import matplotlib.pyplot as plt

with open("/home/zeng/Desktop/log_test.log") as f:
    contents = f.read()

regex = re.compile("After this iteration current coverage is: .*")
results = re.findall(regex, contents)

results = [float(re.search('0.[0-9]*', x).group()) for x in results]

plt.title("Coverage v.s. iteration")
plt.plot(results)
plt.savefig("results.jpg")

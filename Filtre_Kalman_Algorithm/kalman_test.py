
# this code is for plotting the signal value with Noise + Signal after
# Using Kalman Filter 

import math
from random import random
import matplotlib.pyplot as plt
import numpy as np


input_file = open("noisy.txt", "r")
result_file = open("noisy_resultat.txt","r")

y1 = []
values_count = 0
for value in input_file:
    #print(value)
    y1.append(float(value.strip()))
    values_count+=1
print("Signal + Noise")

print(y1)

y2 = []
for value in result_file:
  #print(value)
  y2.append(float(value.strip()))

print("Signal after using Kalman Filter")
print(y2)

step = 1/(values_count/5)
x1 = np.arange(0, 5, step).tolist()
print(x1)

# plotting the points
plt.plot(x1, y1,label='signal_avec_noise')
plt.plot(x1, y2,label='filtered')
plt.title("Kalman Filter Example")
plt.ylabel("Temperature Sensor Reading (°C)")
plt.xlabel("Time (s)")
plt.legend()
plt.show()

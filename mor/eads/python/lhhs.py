import random
iterations = 10
segmentSize = 1.0 / iterations
for i in range(iterations):
    segmentMin = i * segmentSize
    point = segmentMin + (random.random() * segmentSize)
    pointValue = (point * (150 - 0.02)) + 0.02
    print " S = ", pointValue

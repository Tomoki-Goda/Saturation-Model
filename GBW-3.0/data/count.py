#!/usr/bin/env python3

import argparse
import sys


#-------------------------------------------------------------------------------
def main():
  data = sys.stdin.readlines()

  counter = 1
  for line in data:
    numbers = line.split()
    numbers = [float(x) for x in numbers]

    # conditions
    q2 = numbers[0]
    x = numbers[1]
    if q2 <= 10 and x <= 0.01:
      print(*numbers)
      counter += 1

  print(counter)

main()

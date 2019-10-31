def power_sum(*numbers):
    sum = 0
    for n in numbers:
        sum = sum + n*n
    return sum


a = power_sum(1, 2, 3)
print(a)

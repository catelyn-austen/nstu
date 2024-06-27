import random

def mod_pow(x, y, m):
    result = 1
    while y > 0:
        if y % 2 == 1:
            result = (result * x) % m
        x = (x * x) % m
        y //= 2
    return result

def diffie_hellman():
    hex_string_p = input("p: ")
    p = int(hex_string_p, 16)
    g = int(input("g: "))
    print("Here is your decimal p: ", p)

    a = random.randint(1, p - 1)
    print("a = ", a)

    A = mod_pow(g, a, p)
    hex_A = format(A, "x")
    hex_A = str(hex_A).upper()

    print("A = ", A)
    print("Here is your hexagon A: ", hex_A)

    B = int(input("B: "))

    S = mod_pow(B, a, p)
    hex_S = format(S, "x")
    hex_S = str(hex_S).upper()
    print("S =", S)
    print("Here is your hexagon S: ", hex_S)

diffie_hellman()

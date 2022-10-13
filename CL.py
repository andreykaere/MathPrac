import time
start_time = time.time()

def is_square(x):
    # ? correct range?
    for i in range(x):
        if (i * i) == x:
            return True
    return False

def is_square_mod(x, a):
    # ? check?
    if a == 0:
        return True
    # ? correct range?
    for i in range(1, abs(a) + 2):
        if ((i * i) - x) % a == 0:
            return True
    return False



def CL(a, b, c):
    return is_square_mod(b*c, a) and is_square_mod(a*c, b) and is_square_mod(-a*b, c)

def find_solution(a, b, c):

    X = int((b * c) ** (1/2)) + 1
    Y = int((a * c) ** (1/2)) + 1
    Z = int((a * b) ** (1/2)) + 1
    # mb not enough precision
    for x0 in range(X):
        for y0 in range(Y):
            # ? mb improve
            for z0 in range(Z):
                if ((a*x0 * x0 + b*y0*y0 - c * z0*z0) == 0) and ((x0 != 0) or (y0 != 0) or (z0 != 0)):
                    return (x0, y0, z0)
    

def main():
    print("input a, b, c")
    a = int(input())
    b = int(input())
    c = int(input())
    start_time = time.time()
    if (not CL(a, b, c)):
        print("No solution")
    
    else:
        print("------ %s seconds ------" % (time.time() - start_time))
        print("There is a solution, wait....................")
        
        print(find_solution(a, b, c))
    
    print("-------- %s seconds -------" % (time.time() - start_time))
    print()
        

while True:
    main()










from sympy import jacobi_symbol
from sympy import isprime
from math import isqrt

def selfridge(n, R = 1):
  '''
  Choose Q and D = 1 - 4 * Q * R such that j = (D|n) = -1.
  Return j = 0 if a D with 1 < GCD(D, n) < n is found.
  '''
  
  j = 1
  Q = 1 if R == 1 else 0
  while j != -1: 
    Q -= 1 if Q <= 0 else 0
    Q *= -1
    D = 1 - 4 * R * Q
    j = jacobi_symbol(D, n)
    if j == 0 and n != abs(D):
      return Q, D, 0
  return Q, D, j


def mod_mat_pow_1st_col(M11, M12, M21, M22, k, n):
  '''
  Given the entries of the matrix M = [[M11, M12], [M21, M22]],
  an exponent k and a modulus n, evaluate M^k mod n and return
  the first column.
  '''
  
  a, b, c, d = M11, M12, M21, M22
  digits = bin(k)[2 : ]
  subscript = 1
  for digit in digits[1 : ]:
    M11, M12, M21, M22 = (pow(M11, 2, n) + M12 * M21) % n, \
                         M12 * (M11 + M22) % n, \
                         M21 * (M11 + M22) % n, \
                         (pow(M22, 2, n) + M12 * M21) % n
    subscript *= 2
    if digit == '1':
      M11, M12, M21, M22 = (a * M11 + b * M21) % n, \
                           (a * M12 + b * M22) % n, \
                           (c * M11 + d * M21) % n, \
                           (c * M12 + d * M22) % n
      subscript += 1
  return M11, M21
  

def pseudoprimes(primes):
  '''
  Given a list of integers (obtained from a primality test),
  check the primality of each one using the isprime() module.
  If the test returns False, the input is added to the list
  of pseudoprimes, which is the output.
  Examples:
    LUCAS TEST
      pseudoprimes(lucas(range(10**5))) = OEIS A217120
    DOUBLE LUCAS TEST
      pseudoprimes(double_lucas(range(10**5))) = OEIS A212423
    GENERALIZED LUCAS TEST
      pseudoprimes(generalized_lucas(range(10**3),R=-1)) = [119, 161, 209, 221, 299, 329, 341, 371, 539, 551, 581, 611, 629, 671, 689, 731, 749, 779, 791, 851, 869, 899, 959, 989]
      pseudoprimes(generalized_lucas(range(10**5),R=1)) = OEIS A212423
      pseudoprimes(generalized_lucas(range(2**38),R=2)) = []
    GENERALIZED PELL TEST
      pseudoprimes(generalized_pell(range(10**5),X=2,Y=1)) = [1853, 5473, 5777, 10877, 23323, 27403, 75077]
      pseudoprimes(generalized_pell(range(2**38),X=3,Y=2) = []
  '''

  pseudos = []
  for n in primes:
    if isprime(n) == False:
      pseudos.append(n)
  return pseudos


def lucas(n, P = 1, Q = ''):
  '''
  *** LUCAS TEST ***
  INPUTS:
        / integer to be tested
    n = - path of a file with a number in each line
        \ range(min, max, step) containing the integers to be tested
        
    P, Q  = non-zero integer parameters for Lucas sequence.
      or
    empty = Selfridge method will be applied.
          
  OUTPUTS:
    True or False if n is a single integer
      or
    list of primes if n is a range of integers or a file path

  Examples:
    lucas(range(3,8,2)) = test 3, 5, 7 with Lucas-Selfridge
    lucas(11,P=3,Q=2) = test 11 with Lucas, for P = 3, Q = 2
  '''

  # n is a single integer
  if isinstance(n, int):                
    
    # exclude 0, 1, 2, 3, even and square numbers
    if n == 0 or n == 1:                  
      return False
    if n == 2 or n == 3:
      return True
    if n % 2 == 0 or isqrt(n)**2 == n:
      return False
    
    # obtain P, Q, D and j = (D|n)
    if isinstance(Q, int):  # P, Q given
      D = P**2 - 4 * Q
      j = jacobi_symbol(D, n)
    else:                   # no parameters in input, use Selfridge method
      Q, D, j = selfridge(n)
      
    # preliminary controls
    if P != 0 and Q != 0 and D == 0:    
      return 'D = 0, change P and/or Q.'
    if D % n == 0 and D != 0:
      return 'D is a multiple of ' + str(n) + ', so it does not pass the test.'
    if j == 0:
      return False
    if Q % n == 0:
      return 'Q is a multiple of ' + str(n) + ', so it does not pass the test.'
    
    # start the test
    U, _ = mod_mat_pow_1st_col(P, -Q, 1, 0, n - j - 1, n)
    if U == 0:
      return True
    return False

  # n is a txt file
  elif isinstance(n, str):                
    
    primes = []
    with open(n, 'r') as file:
      for line in file:
        i = int(line[ : -1])
        res = lucas(i, P = P, Q = Q)
        if res:
          primes.append(i)
    return primes

  # n is a range
  elif isinstance(n, range):                                 
    
    primes = []
    for i in n:
      res = lucas(i, P = P, Q = Q)
      if res:
        primes.append(i)
    return primes


def double_lucas(n, P = 1, Q = ''):
  '''
  *** DOUBLE LUCAS TEST ***
  INPUTS:
        / integer to be tested
    n = - path of a file with a number in each line
        \ range(min, max, step) containing the integers to be tested
        
    P, Q  = non-zero integer parameters for Lucas sequence
      or
    empty = Selfridge method will be applied
          
  OUTPUTS:
    True or False if n is a single integer
      or
    list of primes if n is a range of integers or a file path

  Examples:
    double_lucas(range(3,8,2)) = test 3, 5, 7 with double Lucas-Selfridge
    double_lucas(11,P=1,Q=2) = test 11 with double Lucas, for P = 1, Q = 2
  '''

  # n is a single integer
  if isinstance(n, int):

    # exclude 0, 1, 2, 3, even and square numbers
    if n == 0 or n == 1: 
      return False
    if n == 2 or n == 3:
      return True
    if n % 2 == 0 or isqrt(n)**2 == n:
      return False

    # obtain P, Q, D and j=(D|n)
    if isinstance(Q, int):  # P, Q given
      D = P**2 - 4 * Q
      j = jacobi_symbol(D, n)
    else:                   # no parameters in input, use Selfridge method
      Q, D, j = selfridge(n)

    # preliminary controls
    if P != 0 and Q != 0 and D == 0:    
      return 'D = 0, change P and/or Q.'
    if D % n == 0 and D != 0:
      return 'D is a multiple of ' + str(n) + ', so it does not pass the test.'
    if j == 0:
      return False
    if Q % n == 0:
      return 'Q is a multiple of ' + str(n) + ', so it does not pass the test.'
    
    # start the test
    U1, U0 = mod_mat_pow_1st_col(P, -Q, 1, 0, n - j, n)
    if j == -1 and U0 == 0 and U1 == Q % n:
      return True
    if j == 1 and U0 == 0 and U1 == 1:
      return True
    return False

  # n is a txt file
  elif isinstance(n, str):                
    
    primes = []
    with open(n, 'r') as file:
      for line in file:
        i = int(line[ : -1])
        res = double_lucas(i, P = P, Q = Q)
        if res:
          primes.append(i)
    return primes

  # n is a range
  elif isinstance(n, range):
    
    primes = []
    for i in n:
      res = double_lucas(i, P = P, Q = Q)
      if res:
        primes.append(i)
    return primes


def generalized_lucas(n, P = 1, Q = '', R = ''):
  '''
  *** GENERALIZED LUCAS TEST ***
  INPUTS:
        / integer to be tested
    n = - path of a file with a number in each line
        \ range(min, max, step) containing the integers to be tested
        
    P, Q, R = non-zero integer parameters for Lucas sequence
      or
    only R  = adapted Selfridge method will be applied
          
  OUTPUTS:
    True or False if n is a single integer
      or
    list of primes if n is a range of integers or a file path

  Examples:
    generalized_lucas(range(3,8,2),R=2) = test 3, 5, 7 with generalized Lucas-Selfridge
    generalized_lucas(11,P=1,Q=2,R=3) = test 11 with generalized Lucas, for P = 1, Q = 2 and R = 3

    Test fermat pseudoprimes with base 2 up to 2**64, list taken from
    http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html
      generalized_lucas('psps-below-2-to-64.txt',R=2) = []
  '''

  # n is a single integer
  if isinstance(n, int):

    # exclude 0, 1, 2, 3, even and square numbers
    if n == 0 or n == 1: 
      return False
    if n == 2 or n == 3:
      return True
    if n % 2 == 0 or isqrt(n)**2 == n:
      return False

    # obtain P, Q, R, D and j=(D|n)
    if isinstance(Q, int):  # P, Q, R given
      D = P**2 - 4 * Q * R
      j = jacobi_symbol(D, n)
    else:                   # only R given, use adapted Selfridge method
      Q, D, j = selfridge(n, R) 

    # preliminary controls
    if P != 0 and Q != 0 and D == 0:    
      return 'D = 0, change P and/or Q.'
    if D % n == 0 and D != 0:
      return 'D is a multiple of ' + str(n) + ', so it does not pass the test.'
    if j == 0:
      return False
    if Q % n == 0:
      return 'Q is a multiple of ' + str(n) + ', so it does not pass the test.'
    
    # start the test
    U1, U0 = mod_mat_pow_1st_col(P, -Q, R, 0, n - j, n)
    if j == -1 and U0 == 0 and U1 == Q * R % n:
      return True
    if j == 1 and U0 == 0 and U1 == 1:
      return True
    return False

  # n is a txt file
  elif isinstance(n, str):                
    
    primes = []
    with open(n, 'r') as file:
      for line in file:
        i = int(line[ : -1])
        res = generalized_lucas(i, P = P, Q = Q, R = R)
        if res:
          primes.append(i)
    return primes

  # n is a range
  elif isinstance(n, range):
    
    primes = []
    for i in n:
      res = generalized_lucas(i, P = P, Q = Q, R = R)
      if res:
        primes.append(i)
    return primes


def generalized_pell(n, X = '', Y = '', D = ''):
  '''
  *** GENERALIZED PELL TEST ***
  INPUTS:
        / integer to be tested
    n = - path of a file with a number in each line
        \ range(min, max, step) containing the integers to be tested
        
    X, Y, D   = parameters for the test
      or
    only X, Y = Selfridge method will be applied
  
  OUTPUTS:
    True or False if n is a single integer
      or
    list of primes if n is a range of integers or a file path

  Examples:
    generalized_pell(range(3,8,2),X=1,Y=3) = test 3,5,7 with Pell-Selfridge, X = 1, Y = 3 
    generalized_pell(11,X=4,Y=2,D=13) = test 11 with Pell, for X = 4, Y = 2, D = 13
    
    Test fermat pseudoprimes with base 2 up to 2**64, list taken from
    http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html
      generalized_pell('psps-below-2-to-64.txt',X=3,Y=2) = []
  '''
  
  # n is a single integer
  if isinstance(n, int):

    # exclude 0, 1, 2, 3, even and square numbers
    if n == 0 or n == 1: 
      return False
    if n == 2 or n == 3:
      return True
    if n % 2 == 0 or isqrt(n)**2 == n:
      return False

    # obtain D and j=(D|n)
    if isinstance(D, int):  # D given
      j = jacobi_symbol(D, n)
    else:                   # D not given, use Selfridge method
      _, D, j = selfridge(n, 1)

    # preliminary controls
    if D % n == 0 and D != 0:
      return 'D is a multiple of ' + str(n) + ', so it does not pass the test.'
    if j == 0:
      return False

    # obtain Q and control its value
    Q = (X ** 2 - D * Y ** 2) % n
    if Q == 0:
      return 'Q is a multiple of ' + str(n) + ', so it does not pass the test.'

    # start the test
    Xk, Yk = mod_mat_pow_1st_col(X, D * Y, Y, X, n - j, n)
    if j == -1 and Xk == Q and Yk == 0:
      return True
    if j == 1 and Xk == 1 and Yk == 0:
      return True
    return False

  # n is a txt file
  elif isinstance(n, str):                
    
    primes = []
    with open(n, 'r') as file:
      for line in file:
        i = int(line[ : -1])
        res = generalized_pell(i, X = X, Y = Y, D = D)
        if res:
          primes.append(i)
    return primes

  # n is a range
  elif isinstance(n, range):
    
    primes = []
    for i in n:
      res = generalized_pell(i, X = X, Y = Y, D = D)
      if res:
          primes.append(i)
    return primes

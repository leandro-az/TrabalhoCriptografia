import operator

from sympy import *
from sympy.abc import x
import numpy as np
import scipy.linalg
from numpy.linalg import inv
from sympy.polys import rem

x1 = Symbol('x1')
x2 = Symbol('x2')
x3 = Symbol('x3')
x4 = Symbol('x4')
x5 = Symbol('x5')
x = Symbol('x')



def mdc(num1, num2):
    resto = None
    while resto is not 0:
        resto = num1 % num2
        num1  = num2
        num2  = resto

    return num1

def invmodp(a,p):

    resposta=0;
    for d in range(1, p):
        r = (int) ((d * a) % p)
        if r == 1:
            resposta=d


    return resposta


def modulo_lista(lista, modulo):
    return [item % modulo for item in lista]


def multiplica_matriz_por_vetor(matriz_a, vetor_b):
    if isinstance(matriz_a, list):
        matriz_a = np.matrix(matriz_a)
    if isinstance(vetor_b, list):
        vetor_b = np.matrix(vetor_b).T
    return (matriz_a * vetor_b).T.tolist()[0]


def transpor_lista(lista):
    if isinstance(lista, list):
        return np.matrix(lista).T.tolist()
    raise TypeError("Argumento precisa ser do tipo list!")


# y = Ax + b
def funcao_linear(matriz_a, x, b):
    return list(map(operator.add, multiplica_matriz_por_vetor(matriz_a, x), b))

############## primo escolhido #################

primo=31

############ Ordem do corpo finito  ############

q = 2

########### Tamanho do corpo finito ############

n = 5

########### Polinômio Sobre o corpo finito de Galois de grau 2 ############

f = (x**5 + x**4 + x**3 + x + 1 )

########### Definição do valor de t ########################

t= 3

############ Definição do valor de h ########################

h =(q**t) +1

############ Definição do valor de h' ########################

h_linha = invmodp(h,primo)

#####################################  verifica primalidade ####################################################

if(mdc(h,(q**n)-1)==1):
  pass

############# Definicao de Ma ###############################

Ma = [[1, 0, 1, 1, 0],
      [0, 1, 1, 0, 1],
      [1, 1, 0, 0, 1],
      [0, 1, 0, 1, 0],
      [0, 0, 0, 1, 1 ]]

############# Definicao de Mb ###############################

Mb = [[1, 0, 0, 1, 1],
      [0, 0, 1, 1, 0],
      [1, 1, 0, 0, 1],
      [1, 1, 0, 0, 0],
      [1, 0, 0, 0, 0]]

############# vetores  constantes #########################

c = [1, 0, 1, 1, 1]
d = [1, 0, 1, 0, 0]

vetX=[x1, x2, x3, x4, x5]

# Definição da matriz A inversa
Mainv = Matrix(Ma).inv_mod(2).tolist()

# Definição da matriz B inversa
Mbinv = Matrix(Mb).inv_mod(2).tolist()

############## definição do vetor intermediário u ############

u = funcao_linear(Ma, vetX, c)

ux=expand((u[0]+ expand(u[1]*(x)) + expand(u[2]*(x**2)) +  expand(u[3]*(x**3)) + expand(u[4]*(x**4)) )*
          (u[0] + expand(u[1] * (x**8)) + expand(u[2] * (x ** 16)) + expand(u[3] * (x ** 24)) + expand(u[4] * (x ** 32))))


q,r = div(ux,f,domain = GF(2))

f0= r + expand((-1)*(x**4)*(r.coeff(x**4))) + expand((-1)*(x**3)*(r.coeff(x**3))) + expand((-1)*(x**2)*(r.coeff(x**2))) + expand((-1)*(x)*(r.coeff(x)))
f4=r.coeff(x**4)
f3=r.coeff(x**3)
f2=r.coeff(x**2)
f1=r.coeff(x)

funcao = [f0, f1, f2, f3, f4]

# Cria o vetor v
v = [Poly(fi, domain=GF(2)) for fi in funcao]

# Subtrai o vetor d do vetor v
v_d = list(map(operator.sub, v, d))

y = multiplica_matriz_por_vetor(Mbinv, v_d)

# Reaplica o campo de Galois ao resultado (senão fica domínio ZZ)
for i, poli in enumerate(y):
    y[i] = poli.set_domain(GF(2))


print("PlainText >>>>" + str([1,1,1,1,1]))
print("\n")

cif=[y[0](1,1,1,1,1),y[1](1,1,1,1,1),y[2](1,1,1,1,1),y[3](1,1,1,1,1),y[4](1,1,1,1,1)]
print("Ciphertext >>>>> " + str(cif))



############################################ decriptacao  #################################################################


v_linha = funcao_linear(Mb, cif, d)

vx=(expand((v_linha[0] + expand(v_linha[1] * (x)) + expand(v_linha[2] * (x**2)) + expand(v_linha[3] * (x**3)) + expand(v_linha[4] * (x**4))) ** h_linha))


qr, vr = div(vx, f, domain = GF(2))

u_linha=[]

u_linha.append(int(vr + expand((-1)*(x**4)*(vr.coeff(x**4))) + expand((-1)*(x**3)*(vr.coeff(x**3))) + expand((-1)*(x**2)*(vr.coeff(x**2))) + expand((-1)*(x)*(vr.coeff(x)))))
u_linha.append(int(vr.coeff(x)))
u_linha.append( int(vr.coeff(x**2)))
u_linha.append( int(vr.coeff(x**3)))
u_linha.append(int(vr.coeff(x**4)))



u_linha_c_linha = list(map(operator.sub, u_linha, c))



x_linha = Mbinv = modulo_lista(multiplica_matriz_por_vetor(Mainv, u_linha_c_linha), 2)

print("\n")

print("Back to PlainText >>>> " + str(x_linha))

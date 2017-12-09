from  sympy import *
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

def multiplicaMatriz(A,B):
    linhaA=len(A)
    colunaA=len(A[0])
    linhaB = len(B)
    colunaB= len(B[0])
    resp = []
    if(colunaA != linhaB):
        print("Impossivel Multiplicar Matrizes colunaA dif linhaB")
        return resp
    for l in range(linhaA):
        resp.append([])
        for c in range(colunaB):
            resp[l].append(0)
            for k in range(colunaA):
                resp[l][c]+=A[l][k]*B[k][c]

    return resp


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

Ma = [[1, 0, 1, 1, 0  ], [0, 1, 1, 0, 1], [1, 1, 0, 0, 1 ], [0, 1, 0, 1, 0], [0, 0, 0, 1, 1 ]]

############# Definicao de Mainv ###############################

Mainv=[[[0, 0, 1, 1, 1  ], [1, 1, 1, 1, 0 ], [0, 1, 0, 1, 1 ], [1, 1, 1, 0, 0], [1, 1, 1, 0, 1 ]]]

############# Definicao de Mb ###############################

Mb = [[1, 0, 0, 1, 1 ], [0, 0, 1, 1, 0], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 0, 0, 0]]

############# Definicao de Mbinv ###############################

Mbinv = [[0, 0, 0, 0, 1 ], [0, 0, 0, 1, 1], [1, 1, 1, 1, 1], [1, 0, 1, 1, 1 ], [0, 0, 1, 1, 0]]

############# vetores  contantes #########################

c = [1, 0, 1, 1, 1]
d = [1, 0, 1, 0, 0]

vetX=[x1, x2, x3, x4, x5]

############## definição do vetor intermediário u ############


u=[]
u.append(((Ma[0][0]*vetX[0])+ (Ma[0][1]*vetX[1]) + (Ma[0][2]*vetX[2]) + (Ma[0][3]*vetX[3]) + (Ma[0][4]*vetX[4]) + c[0]))
u.append(((Ma[1][0]*vetX[0])+ (Ma[1][1]*vetX[1]) + (Ma[1][2]*vetX[2]) + (Ma[1][3]*vetX[3]) + (Ma[1][4]*vetX[4]) + c[1]))
u.append(((Ma[2][0]*vetX[0])+ (Ma[2][1]*vetX[1]) + (Ma[2][2]*vetX[2]) + (Ma[2][3]*vetX[3]) + (Ma[2][4]*vetX[4])+ c[2]))
u.append(((Ma[3][0]*vetX[0])+ (Ma[3][1]*vetX[1]) + (Ma[3][2]*vetX[2]) + (Ma[3][3]*vetX[3]) + (Ma[3][4]*vetX[4])+ c[3]))
u.append(((Ma[4][0]*vetX[0])+ (Ma[4][1]*vetX[1]) + (Ma[4][2]*vetX[2]) + (Ma[4][3]*vetX[3]) + (Ma[4][4]*vetX[4])+ c[4]))

ux=expand((u[0]+ expand(u[1]*(x)) + expand(u[2]*(x**2)) +  expand(u[3]*(x**3)) + expand(u[4]*(x**4)) )*
          (u[0] + expand(u[1] * (x**8)) + expand(u[2] * (x ** 16)) + expand(u[3] * (x ** 24)) + expand(u[4] * (x ** 32))))


q,r = div(ux,f,domain = GF(2))

f0= r + expand((-1)*(x**4)*(r.coeff(x**4))) + expand((-1)*(x**3)*(r.coeff(x**3))) + expand((-1)*(x**2)*(r.coeff(x**2))) + expand((-1)*(x)*(r.coeff(x)))
f4=r.coeff(x**4)
f3=r.coeff(x**3)
f2=r.coeff(x**2)
f1=r.coeff(x)


v=[]
v.append( Poly(f0,domain = GF(2)))
v.append( Poly(f1,domain = GF(2)))
v.append( Poly(f2,domain = GF(2)))
v.append( Poly(f3,domain = GF(2)))
v.append( Poly(f4,domain = GF(2)))






v_d=[]
v_d.append(v[0]-d[0])
v_d.append(v[1]-d[1])
v_d.append(v[2]-d[2])
v_d.append(v[3]-d[3])
v_d.append(v[4]-d[4])


y=[]

y.append(Poly((Mbinv[0][0]*v_d[0])+ (Mbinv[0][1]*v_d[1]) + (Mbinv[0][2]*v_d[2]) + (Mbinv[0][3]*v_d[3]) + (Mbinv[0][4]*v_d[4]) ,domain = GF(2)))
y.append(Poly((Mbinv[1][0]*v_d[0])+ (Mbinv[1][1]*v_d[1]) + (Mbinv[1][2]*v_d[2]) + (Mbinv[1][3]*v_d[3]) + (Mbinv[1][4]*v_d[4]) ,domain = GF(2)))
y.append(Poly((Mbinv[2][0]*v_d[0])+ (Mbinv[2][1]*v_d[1]) + (Mbinv[2][2]*v_d[2]) + (Mbinv[2][3]*v_d[3]) + (Mbinv[2][4]*v_d[4]) ,domain = GF(2)))
y.append(Poly((Mbinv[3][0]*v_d[0])+ (Mbinv[3][1]*v_d[1]) + (Mbinv[3][2]*v_d[2]) + (Mbinv[3][3]*v_d[3]) + (Mbinv[3][4]*v_d[4]) ,domain = GF(2)))
y.append(Poly((Mbinv[4][0]*v_d[0])+ (Mbinv[4][1]*v_d[1]) + (Mbinv[4][2]*v_d[2]) + (Mbinv[4][3]*v_d[3]) + (Mbinv[4][4]*v_d[4]) ,domain = GF(2)))



print("y1>>"+ str(y[1]))
print("v5>>" + str(v[4]))



cif=[y[0](1,1,1,1,1),y[1](1,1,1,1,1),y[2](1,1,1,1,1),y[3](1,1,1,1,1),y[4](1,1,1,1,1)]
print(">>>>> cifrado:" + str(cif))





dec=[]
dec.append(((Mb[0][0]*cif[0]) + (Mb[0][1]*cif[1]) + (Mb[0][2]*cif[2]) + (Mb[0][3]*cif[3]) + (Mb[0][4]*cif[4])+ d[0])%2)
dec.append(((Mb[1][0]*cif[0]) + (Mb[1][1]*cif[1]) + (Mb[1][2]*cif[2]) + (Mb[1][3]*cif[3]) + (Mb[1][4]*cif[4])+ d[1])%2)
dec.append(((Mb[2][0]*cif[0]) + (Mb[2][1]*cif[1]) + (Mb[2][2]*cif[2]) + (Mb[2][3]*cif[3]) + (Mb[2][4]*cif[4])+ d[2])%2)
dec.append(((Mb[3][0]*cif[0]) + (Mb[3][1]*cif[1]) + (Mb[3][2]*cif[2]) + (Mb[3][3]*cif[3]) + (Mb[3][4]*cif[4])+ d[3])%2)
dec.append(((Mb[4][0]*cif[0]) + (Mb[4][1]*cif[1]) + (Mb[4][2]*cif[2]) + (Mb[4][3]*cif[3]) + (Mb[4][4]*cif[4])+ d[4])%2)


Ma = np.array( [[1, 0, 1, 1, 0  ], [0, 1, 1, 0, 1], [1, 1, 0, 0, 1 ], [0, 1, 0, 1, 0], [0, 0, 0, 1, 1 ]])


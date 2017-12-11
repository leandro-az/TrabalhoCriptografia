import operator

import funcoes_utilitarias as utilitarias
from sympy import *
from sympy.abc import x

if __name__ == '__main__':
    # Definição dos pontos no campo como símbolos do sympy
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    x3 = Symbol('x3')
    x4 = Symbol('x4')
    x5 = Symbol('x5')

    # Primo usado
    primo = 31

    # Ordem do campo finito
    q = 2

    # Tamanho (dimensões) do campo finito
    n = 5

    # Polinômio a ser usado sobre o corpo finito de Galois de grau 2
    f = (x ** 5 + x ** 4 + x ** 3 + x + 1)

    # Definição do valor de t
    t = 3

    # Definição do valor de h
    h = (q ** t) + 1

    # Definição do valor de h'
    h_linha = utilitarias.inverso_modular(h, primo)

    # Verifica se h e (q^n) - 1 são coprimos
    if utilitarias.mdc(h, (q ** n) - 1) == 1:
        pass

    # Definição da matriz A e inversa
    matriz_a = [[1, 0, 1, 1, 0],
                [0, 1, 1, 0, 1],
                [1, 1, 0, 0, 1],
                [0, 1, 0, 1, 0],
                [0, 0, 0, 1, 1]]

    matriz_a_inversa = Matrix(matriz_a).inv_mod(2).tolist()

    # Definição da matriz B e inversa
    matriz_b = [[1, 0, 0, 1, 1],
                [0, 0, 1, 1, 0],
                [1, 1, 0, 0, 1],
                [1, 1, 0, 0, 0],
                [1, 0, 0, 0, 0]]

    matriz_b_inversa = Matrix(matriz_b).inv_mod(2).tolist()

    # Vetores constantes
    c = [1, 0, 1, 1, 1]
    d = [1, 0, 1, 0, 0]

    vetX = [x1, x2, x3, x4, x5]

    # Definição da matriz B inversa

    # Definição do vetor intermediário u
    u = utilitarias.funcao_linear(matriz_a, vetX, c)

    ux = expand((u[0] + expand(u[1] * x) +
                 expand(u[2] * (x ** 2)) +
                 expand(u[3] * (x ** 3)) +
                 expand(u[4] * (x ** 4))) *
                (u[0] + expand(u[1] * (x ** 8)) +
                 expand(u[2] * (x ** 16)) +
                 expand(u[3] * (x ** 24)) +
                 expand(u[4] * (x ** 32))))

    # Divide ux pelo polinômio f
    q, r = div(ux, f, domain=GF(2))

    f0 = r + expand((-1) * (x ** 4) * (r.coeff(x ** 4))) + expand((-1) * (x ** 3) * (r.coeff(x ** 3))) + expand(
        (-1) * (x ** 2) * (r.coeff(x ** 2))) + expand((-1) * x * (r.coeff(x)))
    f4 = r.coeff(x ** 4)
    f3 = r.coeff(x ** 3)
    f2 = r.coeff(x ** 2)
    f1 = r.coeff(x)

    funcao = [f0, f1, f2, f3, f4]

    # Cria o vetor v
    v = [Poly(fi, domain=GF(2)) for fi in funcao]

    # Subtrai o vetor d do vetor v
    v_d = list(map(operator.sub, v, d))

    # Multiplica
    y = utilitarias.multiplica_matriz_por_vetor(matriz_b_inversa, v_d)

    # Reaplica o campo de Galois ao resultado (senão fica domínio ZZ)
    for i, poli in enumerate(y):
        y[i] = poli.set_domain(GF(2))

    # Define a mensagem que desejamos cifrar
    mensagem_original = [1, 1, 1, 1, 1]

    print("PlainText >>>>" + str(mensagem_original))
    print("\n")

    # Cifra a mensagem (usando o polinômio em y)
    mensagem_cifrada = [i.eval(tuple(mensagem_original)) for i in y]

    print("Ciphertext >>>>> " + str(mensagem_cifrada))

    # Decifragem da mensagem
    v_linha = utilitarias.funcao_linear(matriz_b, mensagem_cifrada, d)

    vx = (expand((v_linha[0] + expand(v_linha[1] * x) + expand(v_linha[2] * (x ** 2)) + expand(
        v_linha[3] * (x ** 3)) + expand(v_linha[4] * (x ** 4))) ** h_linha))

    qr, vr = div(vx, f, domain=GF(2))

    u_linha = [int(
        vr +
        expand((-1) * (x ** 4) * (vr.coeff(x ** 4))) +
        expand((-1) * (x ** 3) * (vr.coeff(x ** 3))) +
        expand((-1) * (x ** 2) * (vr.coeff(x ** 2))) +
        expand((-1) * x * (vr.coeff(x)))),
        int(vr.coeff(x)), int(vr.coeff(x ** 2)), int(vr.coeff(x ** 3)), int(vr.coeff(x ** 4))]

    u_linha_c_linha = list(map(operator.sub, u_linha, c))

    mensagem_decifrada = utilitarias.modulo_lista(utilitarias.multiplica_matriz_por_vetor(matriz_a_inversa, u_linha_c_linha), 2)

    print("\n")

    print("Back to PlainText >>>> " + str(mensagem_decifrada))

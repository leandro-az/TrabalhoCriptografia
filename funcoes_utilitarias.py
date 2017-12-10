import numpy as np
import operator


def mdc(num1, num2):
    resto = None
    while resto is not 0:
        resto = num1 % num2
        num1 = num2
        num2 = resto
    return num1


def inverso_modular(valor, modulo):
    resposta = 0
    for d in range(1, modulo):
        r = int((d * valor) % modulo)
        if r == 1:
            resposta = d
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

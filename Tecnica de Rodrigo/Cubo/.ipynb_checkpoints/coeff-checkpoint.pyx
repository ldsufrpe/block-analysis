from functools import reduce
from sage.all import binomial
def dicionario_geral(poly):
    """INPUT - poly- Um polinomio como expressão simbolica na classe <class 'sage.symbolic.expression.Expression'>
    OUTPUT - Um dicionario com CHAVE do tipo tupla multigrau (d1,d2,...,dl) correpondendo ao monomio 
    variavel1^(d1)*variavel2^(d2)*...*variavell^(dl) e VALOR sendo o coeficiente deste monomio.
    OBSERVACAO- A ordem das variaveis imspota pelo SAGE manda na ordem das chaves do dicionario. 
    Exemplo-
    sage: z, y, x, u = var('z y x u')
    sage: p = (2100*u^2*x+2020*u^2*y^2)*z^0+(3031*u^3*y^3+1101*u*x)*z+(3022*u^3*y^2+202*x^2)*z^2
    sage: variaveis = p.free_variables(); variaveis
        (u,x,y,z)
    sage: dicionario_geral(p)
        {(0, 2, 0, 2): 202,
         (1, 1, 0, 1): 1101,
         (2, 0, 2, 0): 2020,
         (2, 1, 0, 0): 2100,
         (3, 0, 2, 2): 3022,
         (3, 0, 3, 1): 3031}
        
        
    sage: u,v = var('u v')
    sage: p3 = 23*u^2*v^3+21*u^2*v+11*u*v+10*u+4;
    sage: dicionario_geral(p3)
      {(0, 0): 4, (1, 0): 10, (1, 1): 11, (2, 1): 21, (2, 3): 23}
    """
    def nova_lista(lista, var):
        new_lista = [];
        for item in lista:
            sublista = item[0].coefficients(var)
            for subitem in sublista:
                subitem = [subitem[0],subitem[1:]+item[-1]]
                new_lista.append(subitem)
        return  new_lista

    #variaveis = (poly).free_variables(); Modifiquei aqui
    variaveis = (poly).variables()
#####################################################################
    variavel_1 = variaveis[0];
    l_variavel_1 = (poly).coefficients(variavel_1);
    l = [[item[0],item[1:]] for item in l_variavel_1]
    for var in variaveis[1:]:
        l = nova_lista(l, var)
    dicionario = {tuple(list(reversed(j[1]))) : j[0] for j in l}
    return dicionario


from functools import reduce
from libc.math cimport pow

def coeff_de_moebius(poly, graus, graus_total, chaves, extremos, novo_coefsxy):    
    
    
    soma_total = 0

    binomials_cache = {}
    powers_cache = {}

    for chave in chaves:
        somas = []
        for k in range(len(chave)):
            ck = chave[k]
            gt = graus_total[k]
            gks = graus[k]
            a_k = extremos[k][0]
            b_k = extremos[k][1]

            somak = 0
            for pp in range(ck + 1):
                binom_ck_pp = binomials_cache.get((ck, pp))
                if binom_ck_pp is None:
                    binom_ck_pp = binomial(ck, pp)
                    binomials_cache[(ck, pp)] = binom_ck_pp

                binom_gtck_gkspp = binomials_cache.get((gt - ck, gks - pp))
                if binom_gtck_gkspp is None:
                    binom_gtck_gkspp = binomial(gt - ck, gks - pp)
                    binomials_cache[(gt - ck, gks - pp)] = binom_gtck_gkspp

                a_k_pp = powers_cache.get((a_k, pp))
                if a_k_pp is None:
                    a_k_pp = pow(a_k, pp)
                    powers_cache[(a_k, pp)] = a_k_pp

                b_k_ck_pp = powers_cache.get((b_k, ck - pp))
                if b_k_ck_pp is None:
                    b_k_ck_pp = pow(b_k,(ck - pp))
                    powers_cache[(b_k, ck - pp)] = b_k_ck_pp

                somak += binom_ck_pp * binom_gtck_gkspp * a_k_pp * b_k_ck_pp

            somas.append(somak)
        somas_e_valor = somas + [novo_coefsxy[chave]]
        soma_total += reduce(lambda a, b: a * b, somas_e_valor)

    return soma_total

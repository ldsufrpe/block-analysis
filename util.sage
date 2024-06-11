y=var('y')

var('c', domain='real')
var('t', domain='real')
assume(t>0,t<1)
assume(c>0)

def simplify_expression(p, repeat=2):
    c = var('c')    
    i=0
    while i<repeat:
        p = p.collect(c).combine(deep=True).simplify_real()
        lst0 = p.coefficients(c)    
        lst = []    
        for l0 in lst0:
            expr0, k = l0
            lst1 = expr0.operands()
            if len(lst1)<2:
                lst.append(expr0*c^k)
            for expr1 in lst1:
                try:
                    #expr1 = expr1.canonicalize_radical().factor()
                    expr1 = expr1.expand().simplify_full()
                except AttributeError:
                    continue
                #if len(p.variables())==1:
                #    pass
                expr1 = expr1*c^k
                lst.append(expr1)
        p = sum(lst).collect(c).combine(deep=True)
        i+=1
    return p

def count_sign_changes(p):
    l = [c for c in p if not c.is_zero()]
    changes = [l[i]*l[i + 1] < 0 for i in range(len(l) - 1)]

    return changes.count(True)

def sturm(p, a, b):
    assert p.degree()> 2
    assert not (p(a) == 0)
    assert not (p(b) == 0)
    assert a <= b
    remains = [p, diff(p)]
    
    for i in range(p.degree() - 1):
        r = -(remains[i] % remains[i + 1])
        if r.is_zero():
            #print (r)
            break
            
        remains.append(r)
    evals = [[], []]
    for q in remains:
        evals[0].append(q(a))
        evals[1].append(q(b))
    return count_sign_changes(evals[0])-count_sign_changes(evals[1])


def expr_to_poly(p, ring=AA):    
    cc = p.coefficients()
    Q = 0
    y = polygen(ring, 'y')
    for c, i in cc:        
        Q+= ring(c)*y^i
    return Q


##Criando operação de substituição de Mobius
def subs_mobius_u_by_x(polinomio, a,b):
    phi=lambda a,b: (a*x+b)/(x+1);
    rational_result=polinomio.subs(u=phi(a,b))
    new_polinomial=rational_result.numerator()
    return new_polinomial

def subs_mobius_v_by_y(polinomio, a,b):
    phi=lambda a,b: (a*y+b)/(y+1);
    rational_result=polinomio.subs(v=phi(a,b))
    new_polinomial=rational_result.numerator()
    return new_polinomial

def is_complete_x(polinomial):
    """Verifica se o polinomio na variavel x eh completo"""
    coeffs=polinomial.coefficients(x);
    grau=polinomial.degree(x)
    if len(coeffs)<(grau+1):
        return False
    else:
        return True
    
def is_complete_y(polinomial):
    """Verifica se o polinomio na variavel y eh completo"""
    coeffs=polinomial.coefficients(y);
    grau=polinomial.degree(y)
    if len(coeffs)<(grau+1):
        return False
    else:
        return True
    
def has_signal_variations_x_y(polinomial):
    """Verifica se o polinomio nas variaveis x e y tem variacoes de sinais"""
    coeffs_y=polinomial.coefficients(y);
    #grau=polinomial.degree(y)
    for i in range(len(coeffs_y)):
        coeffs_x=coeffs_y[i][0].coefficients(x)
        #print(coeffs_x)
        for j in range(len(coeffs_x)):
            if bool(coeffs_x[j-1][0]*coeffs_x[j][0]<0):
                return True

    return False

def expand_poly(poly, var, a, b):
    coeffs = poly.coefficients()
    p = 0
    d = poly.degree(x)
    for i in range(len(coeffs)):
        c = coeffs[d-i]
        p+= c[0]*(a*var+b)^(d-i)*(var+1)^i       
        
    return p


def subs_mobius_2d_uv_errada(poly, a, b, c, d):
    """
    poly ---
    """    
    #rational_result =  expand_poly(poly, u, a, b) #poly(u=(a*u + b)/(u + 1))
    poly = expand_poly(poly, u, a, b)#rational_result.numerator()
    poly = expand_poly(poly, v, c, d) #poly(v=(c*v + d)/(v + 1))
    #poly = rational_result.numerator()
    return poly.numerator()
    
    
def subs_mobius_2d_uv(poly, a, b, c, d):
    """
    poly ---
    """    
    rational_result_u =  poly(u=(a*u + b)/(u + 1))   
    rational_result_uv = rational_result_u(v=(c*v + d)/(v + 1))    
    poly = rational_result_uv.numerator()
    return poly
            


def has_signal_variations_u_v(poly):
    """Verifica se o polinomio nas variaveis x e y tem variacoes de sinais"""
    coeffs_v = poly.coefficients(v);
    
    for i in range(len(coeffs_v)):
        coeffs_u=coeffs_v[i][0].coefficients(u)
        
        for j in range(len(coeffs_u)):
            if bool(coeffs_u[j-1][0]*coeffs_u[j][0]<0):
                return True

    return False


def has_signal_variations_u_v_correted_optimize(poly):
    """Verifica se o polinomio nas variaveis u e v tem variacoes de sinais"""
    
    coeffs_v = poly.coefficients(v);
    coeffs_u_and_grades = [x[0].coefficients(u) for x in coeffs_v];
    coeffs = [item[0] for items in coeffs_u_and_grades for item in items];
    for j in range(1, len(coeffs)):
        if sign(coeffs[j-1])*sign(coeffs[j])<0:
            return True

    return False


#t = var('t')
def extremos_uv_particionados(p , extremos_t):
    """ INPUT - UM POLINOMIO QUADRATICO NA VARIAVEL t DA FORMA
   p --> p(t)=a*t^2 + b*t+c 
   extremos_t ----> Os extremos para o t, (deve ser um subintervalo do INTERVALO [0,1])
   OUTPUT - OS VALORES MAXIMOS E MINIMOS DE p NO INTERVALO [t0,t1].
    """
    t0 =  extremos_t[0]
    t1 =  extremos_t[1]
    coefs = reversed(p.coefficients(t))
    a, b , c = map(lambda x: x[0], coefs )
    
    if t0 <=-b/(2*a) <=t1:
        
        lst = [p(t=t0), p(t=t1), p(t=-b/(2*a))]
        
    else:
        lst = [p(t=t0), p(t=t1)]
        
    return min(lst), max(lst)

###TESTAR SE f e g é um par u e v.
def teste_par_u_v(lista):
    """INPUT- Uma lista com o possivel par f e g para v e u(nesta ordem).
    OUTPUT- Um booleano Falso ou verdadeiro se f e g forma um par para a substituição"""
    f_lista = lista[0]
    
    g_lista = lista[1]
    
    b = (g_lista.coefficients()[1][0])
    
    return bool(((g_lista-f_lista)/(2*b)).full_simplify()== t)


from functools import reduce

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

#novo_coefsxy = dicionario_geral(poly_u_v) # virar variável global
def coeficientes_de_moebius_simbolico_varias_variaveis(poly, graus, extremos, novo_coefsxy ):
    """INPUTS -
    poly- Um polinomio como expressão simbolica na classe  <class 'sage.symbolic.expression.Expression'>;
    graus - Uma lista com naturais representando o multigrau do monomio que voce quer o coeficiente.
            Importante - A sua lista de graus deve estar na ordem que o Sage ordena as variaveis.
    extremos - Uma lista com pares (a_i,b_i) que dão o intervalo da substituição de Moebius para
    a variavel x_i, isto é,    x_i---> x_i = (a_i*u_i+b_i)/(u_i+1) para a variavel u_i.
    OUTPUT - O coeficiente do monomio com multigrau dado pela lista graus depois da substituicao de Moebius
    Exemplo - 1
    sage: poly_teste = 23*u^2*v^3 + 21*u^2*v + 11*u*v + 10*u + 4;
    sage: graus_teste = [1,2]; ## para acahar o coeficiente de (x_1)^1*(x_2)^2
    sage: extremos_teste = [(3, 5), (7, 11)];
    sage: novo_coefsxy_teste = dicionario_geral(poly_teste); # virar variável global;
    sage: c = coeficientes_de_moebius_simbolico_varias_variaveis(poly_teste, graus_teste, extremos_teste, novo_coefsxy_teste)
        1133944
    Exemplo - 2
    sage: poly = (2100*u^2*x+2020*u^2*y^2)*z^0+(3031*u^3*y^3+1101*u*x)*z+(3022*u^3*y^2+202*x^2)*z^2
    sage: variaveis = poly.free_variables(); variaveis
    sage: extremos_teste_2 = [(3,5),(7,11),(13,17),(19,23)];
    sage: ru, rx, ry, rz = var('ru, rx, ry, rz')
    sage: variaveis_r = [ru, rx, ry, rz]
    sage: novo_coefsxy_teste_2 = dicionario_geral(poly); 
    sage: c = coeficientes_de_moebius_simbolico_varias_variaveis(poly, [0, 0, 1, 0], extremos_teste_2, novo_coefsxy_teste_2) ## Procurando o coeficiente de ru^0 * rx^0 * ry^1 * rz^0
    sage: show(c)
     244374065874
       
       OBSERVACAO: Para mais detalhes veja o arquivo /Artigo_de_Aplicacoes/Reuniao-casa-marcelo-26-07/Icosaedro-bloco-2-reunião-casa-marcelo/Funcao_Para_Coeficientes_Com_Varias_Variaveis_Versao_Final.ipynb
        """
    #variaveis = (poly).free_variables(); Modifiquei aqui
    variaveis = (poly).variables()
    graus_total = [poly.degree(x) for x in variaveis];
    chaves = novo_coefsxy.keys()
    soma_total = 0
    for chave in chaves:
        somas = []
        for k in range(len(chave)):
            somak = 0
            ck = chave[k]
            gt = graus_total[k]
            gks = graus[k]
            a_k = extremos[k][0]
            b_k = extremos[k][1]
            for pp in range(ck+1):
                somak = somak+((binomial(ck,pp))*(binomial(gt-ck, gks-pp))*((a_k)^(pp))*((b_k)^(ck-pp)))
            somas.append(somak)
        somas_e_valor = somas + [novo_coefsxy[chave]]
        soma_total = soma_total + reduce(lambda a, b: a*b, somas_e_valor)
    return soma_total

















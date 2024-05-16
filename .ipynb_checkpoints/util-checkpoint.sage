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






















var('c', domain='real')
var('t', domain='real')
assume(t > 0, t < 1)
assume(c > 0)
var('u,v')
assume(u > 0, v > 0)


def simplify_fractional_radicals(expr):
    # Step 1: Separate each term of the expression into a list
    terms = expr.operands()
    # If the expression has no operands, then it is a single term
    if not terms:
        terms = [expr]

    # List to store simplified terms
    simplified_terms = []

    for term in terms:
        # Step 2: Identify which expressions contain a fractional power

        # Step 3: If it is a symbolic fractional expression, separate numerator and denominator
        numerator = term.numerator()
        denominator = term.denominator()

        # Step 4: Identify if the numerator or the denominator or both contain expressions with fractional powers
        if numerator.operator() == operator.pow and not numerator.operands()[1].is_integer():
            base, exp = numerator.operands()
            numerator_simplified = (base.factor() ** exp).simplify_real()
            # print(base, exp)
        else:
            numerator_simplified = numerator

        if denominator.operator() == operator.pow and not denominator.operands()[1].is_integer():
            base, exp = denominator.operands()
            denominator_simplified = (base.factor() ** exp).simplify_real()
            # print(base, exp)
            # print(denominator_simplified )
        else:
            denominator_simplified = denominator

        # Step 5: Reassemble each fractional expression after the previous simplification
        term_simplified = (numerator_simplified / denominator_simplified)

        # Add the simplified term to the list
        simplified_terms.append(term_simplified)

    # Step 7: Reassemble the complete expression after simplifications
    expr_simplified = sum(simplified_terms)

    return expr_simplified


def simplify_expression(p, repeat=2):
    c = var('c')
    i = 0
    while i < repeat:
        p = p.collect(c).combine(deep=True).simplify_real()
        lst0 = p.coefficients(c)
        lst = []
        for l0 in lst0:
            expr0, k = l0
            lst1 = expr0.operands()
            if len(lst1) < 2:
                lst.append(expr0 * c ^ k)
            for expr1 in lst1:
                try:
                    # expr1 = expr1.canonicalize_radical().factor()
                    expr1 = expr1.expand().simplify_full()
                except AttributeError:
                    continue
                # if len(p.variables())==1:
                #    pass
                expr1 = expr1 * c ^ k
                lst.append(expr1)
        p = sum(lst).collect(c).combine(deep=True)
        i += 1
    return p


def count_sign_changes(p):
    l = [c for c in p if not c.is_zero()]
    changes = [l[i] * l[i + 1] < 0 for i in range(len(l) - 1)]

    return changes.count(True)


def sturm(p, a, b):
    assert p.degree() > 2
    assert not (p(a) == 0)
    assert not (p(b) == 0)
    assert a <= b
    remains = [p, diff(p)]

    for i in range(p.degree() - 1):
        r = -(remains[i] % remains[i + 1])
        if r.is_zero():
            # print (r)
            break

        remains.append(r)
    evals = [[], []]
    for q in remains:
        evals[0].append(q(a))
        evals[1].append(q(b))
    return count_sign_changes(evals[0]) - count_sign_changes(evals[1])


def expr_to_poly(p, ring=AA):
    cc = p.coefficients()
    Q = 0
    y = polygen(ring, 'y')
    for c, i in cc:
        Q += ring(c) * y ^ i
    return Q


def polynomial_to_dict(polynomial):
    """
    Converts a polynomial into a dictionary representation.

    **Parameters:**
    - `polynomial`: A polynomial given as a symbolic expression (class `<class 'sage.symbolic.expression.Expression'>`).

    **Returns:**
    - A dictionary where the keys are tuples representing the multidegree (d1, d2, ..., dl) corresponding to the monomial
      `variable1^(d1) * variable2^(d2) * ... * variablel^(dl)`, and the values are the coefficients of these monomials.

    **Note:**
    The order of the variables as imposed by SAGE determines the order of the keys in the dictionary.

    **Examples:**

    sage: z, y, x, u = var('z y x u')
    sage: p = (2100*u^2*x + 2020*u^2*y^2)*z^0 + (3031*u^3*y^3 + 1101*u*x)*z + (3022*u^3*y^2 + 202*x^2)*z^2
    sage: variables = p.variables(); variables
    (u, x, y, z)
    sage: polynomial_to_dict(p)
    {
        (0, 2, 0, 2): 202,
        (1, 1, 0, 1): 1101,
        (2, 0, 2, 0): 2020,
        (2, 1, 0, 0): 2100,
        (3, 0, 2, 2): 3022,
        (3, 0, 3, 1): 3031
    }

    sage: u, v = var('u v')
    sage: p3 = 23*u^2*v^3 + 21*u^2*v + 11*u*v + 10*u + 4
    sage: polynomial_to_dict(p3)
    {
        (0, 0): 4,
        (1, 0): 10,
        (1, 1): 11,
        (2, 1): 21,
        (2, 3): 23
    }

    """

    def update_list(existing_list, variable):
        updated_list = []
        for item in existing_list:
            sublist = item[0].coefficients(variable)
            for subitem in sublist:
                new_item = [subitem[0], subitem[1:] + item[-1]]
                updated_list.append(new_item)
        return updated_list

    variables = polynomial.variables()
    first_variable = variables[0]

    initial_list = polynomial.coefficients(first_variable)
    formatted_list = [[item[0], item[1:]] for item in initial_list]

    for var in variables[1:]:
        formatted_list = update_list(formatted_list, var)

    result_dict = {tuple(list(reversed(entry[1]))): entry[0] for entry in formatted_list}
    return result_dict


from functools import reduce
from sage.all import binomial


def moebius_coefficient(polynomial, degree_per_variable, total_degree, variable_ranges, coefficient_dict):
    """
    Calculates the Möbius coefficient for a given polynomial expression.

    **Parameters:**
    - `poly`: A polynomial as a symbolic expression (class `<class 'sage.symbolic.expression.Expression'>`).
    - `degrees`: A list of natural numbers representing the multidegree of the monomial for which you want the coefficient.
                 **Important:** The order of degrees should follow the order in which Sage orders the variables.
    - `variable_ranges`: A list of pairs `(a_i, b_i)` that give the interval of the Moebius substitution for the variable `x_i`,
                i.e., `x_i --> x_i = (a_i*u_i + b_i)/(u_i + 1)` for the variable `u_i`.
    - `dic_coeffs`: A dictionary where the keys are tuples of degrees and the values are the coefficients of the monomials.

    **Returns:**
    - The Möbius coefficient for the given polynomial expression.

    **Examples:**
    sage: u, x, v, y, z = var('u x v y z')
    sage: p = (2100*u^2*x+2020*u^2*y^2)*z^0+(3031*u^3*y^3+1101*u*x)*z+(3022*u^3*y^2+202*x^2)*z^2
    sage: d = polynomial_to_dict(p)
    sage: vars = p.variables()
    sage: total_degs = [p.degree(x) for x in vars]
    sage: ranges = [(3, 5), (7, 11), (13, 17), (19, 23)]
    sage: calculate_moebius_coefficient(p, [0, 0, 1, 0], total_degs, ranges, d)
    244374065874

    sage: u, v = var('u v')
    sage: p = 23*u^2*v^3 + 21*u^2*v + 11*u*v + 10*u + 4
    sage: vars = p.variables()
    sage: total_degs = [p.degree(x) for x in vars]
    sage: degs = [1, 2]  # To find the coefficient of (x_1)^1 * (x_2)^2
    sage: bounds = [(3, 5), (7, 11)]
    sage: d = polynomial_to_dict(p)
    sage: calculate_moebius_coefficient(p, degs, total_degs, bounds, d)
    1133944
     """
    from sage.all import binomial
    from functools import reduce

    total_sum = 0

    binom_cache = {}
    pow_cache = {}

    for key in coefficient_dict.keys():

        partial_sums = []
        for k in range(len(key)):
            ck = key[k]
            tot_deg_k = total_degree[k]
            deg_k = degree_per_variable[k]

            a_k, b_k = variable_ranges[k]

            partial_sum = 0
            for pp in range(ck + 1):
                binom_ck_pp = binom_cache.get((ck, pp))
                if binom_ck_pp is None:
                    binom_ck_pp = binomial(ck, pp)
                    binom_cache[(ck, pp)] = binom_ck_pp

                binom_tot_deg_diff = binom_cache.get((tot_deg_k - ck, deg_k - pp))
                if binom_tot_deg_diff is None:
                    binom_tot_deg_diff = binomial(tot_deg_k - ck, deg_k - pp)
                    binom_cache[(tot_deg_k - ck, deg_k - pp)] = binom_tot_deg_diff

                a_k_pow = pow_cache.get((a_k, pp))
                if a_k_pow is None:
                    a_k_pow = a_k ^ pp
                    pow_cache[(a_k, pp)] = a_k_pow

                b_k_pow = pow_cache.get((b_k, ck - pp))
                if b_k_pow is None:
                    b_k_pow = b_k ^ (ck - pp)
                    pow_cache[(b_k, ck - pp)] = b_k_pow

                partial_sum += binom_ck_pp * binom_tot_deg_diff * a_k_pow * b_k_pow

            partial_sums.append(partial_sum)
        partial_sums_and_value = partial_sums + [coefficient_dict[key]]
        total_sum += reduce(lambda a, b: a * b, partial_sums_and_value)

    return total_sum


def quadratic_extrema_on_interval(p, t_interval):
    """
    INPUT - A QUADRATIC POLYNOMIAL IN THE VARIABLE t OF THE FORM:
    p --> p(t) = a*t^2 + b*t + c
    t_interval ----> The bounds for t (should be a subinterval of [0,1])

    OUTPUT - THE MINIMUM AND MAXIMUM VALUES OF p IN THE INTERVAL [t0, t1].
    """
    t0 = t_interval[0]
    t1 = t_interval[1]
    coefs = reversed(p.coefficients(t))  # Retrieve the coefficients [a, b, c]
    a, b, c = map(lambda x: x[0], coefs)

    # Check if the critical point -b/(2*a) lies within the interval [t0, t1]
    if t0 <= -b / (2 * a) <= t1:
        values = [p(t=t0), p(t=t1), p(t=-b / (2 * a))]
    else:
        values = [p(t=t0), p(t=t1)]

    return sqrt(min(values)), sqrt(max(values))








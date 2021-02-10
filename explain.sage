load('dynamic_f5.sage')
load('../gb_voodoo/f4_5.sage')

def parse_magma_poly_system(input_file_name, field, order='degrevlex'):
    with open(input_file_name) as input_file:
        content = input_file.read()
    content = "".join(content.split()) # remove all whitespaces
    content = content.replace("[", "").replace("]", "")
    string_system = content.split(",")

    max_var = 0
    decomposed_system = []
    for string_poly in string_system:
        decomposed_poly = []
        string_terms = string_poly.split("+")
        for string_term in string_terms:
            coeff, variables, exponents = 1, [], []
            rest = string_term.split("*")
            if rest[0][0] != "x":
                coeff, rest = int(rest[0]), rest[1:]
            for monom_string in rest:
                rest = monom_string.split("^")
                if len(rest) == 1:
                    rest += [1]
                variables += [rest[0]]
                exponents += [int(rest[1])]
                max_var = max(max_var, int(rest[0][1:]))
            decomposed_poly += [(coeff, variables, exponents)]
        decomposed_system += [decomposed_poly]
    ring = PolynomialRing(field, max_var, "x", order=order)
    exes = ring.gens()
    system = []
    for decomposed_poly in decomposed_system:
        poly = 0
        for coeff, string_variables, exponents in decomposed_poly:
            monom = 1
            variables = [exes[int(var[1:]) - 1] for var in string_variables]
            for var, exp in zip(variables, exponents):
                monom *= var**exp
            poly += coeff*monom
        system += [ring(poly)]
    return system

# —————————————————————————————————————————————————————————————————————————––– #

set_verbose(1)
f5 = F5() # or: F5R(), F5C(), F4F5()
f5_dynamic # this also exists
input_file_name = "./xlix_comps/grevlex_3-round_GF-101_system.txt"
system = parse_magma_poly_system(input_file_name, GF(101), "degrevlex")
G, voos = f5(system, homogenize=False)
[print(f"{v}\n") for v in voos]
involvement_buckets = [0] * 2**(len(voos[0]))
for voo in voos:
    idx = sum([2**i for i, x in enumerate(voo) if x])
    involvement_buckets[idx] += 1
print(f"––––––––––––\n Involvement: {involvement_buckets}")


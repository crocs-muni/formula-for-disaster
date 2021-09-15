from pyecsca.ec.mod import Mod, SymbolicMod
from pyecsca.ec.coordinates import AffineCoordinateModel
from pyecsca.ec.params import _create_params
from pyecsca.ec.model import ShortWeierstrassModel, MontgomeryModel, EdwardsModel, TwistedEdwardsModel
from pyecsca.ec.point import Point

from pyecsca.misc.cfg import TemporaryConfig

import sys
import click
import json

from sympy import symbols, fraction, sympify, FF
from sympy.parsing.sympy_parser import parse_expr
xP, xQ, xPQ, yP, yQ, yPQ = symbols('xP xQ xPQ yP yQ yPQ')


# ===========================================================================================================================
# json2curves read a list of curves in a given json file
def json2curves(filename: str):
    f = open(filename)
    curves = json.load(f)
    f.close()
    return curves, len(curves['curves'])


# ===========================================================================================================================
# get_curve return the curve structure of the i-th element of a given list of curves
def get_curve(curves : dict, i : int, coords: str):
    return _create_params(curves['curves'][i], coords, infty=False)


# ===========================================================================================================================
# converting from point structure to string (for a pretty print of the counterexamples)
def point2str(name : str, point, b : int):
    strpoint = name + '={'
    for c, v in point.coords.items():
        strpoint += c + ': 0x%0*x, ' % (b, int(v))
    return strpoint[:-2] + '}'


def algebraic_point(params, x, y):
    # --- Projective point representation of P
    P_expr = {var: Mod(1, params.curve.prime) for var in params.curve.coordinate_model.variables}
    if 'X' in params.curve.coordinate_model.variables:
        P_expr['X'] = SymbolicMod(x, params.curve.prime)
    else:
        if params.curve.model.shortname == "edwards":
            if params.curve.coordinate_model.name == "yz":
                P_expr['Y'] = SymbolicMod(y * params.curve.parameters['r'], params.curve.prime)
            elif params.curve.coordinate_model.name == "yzsquared":
                P_expr['Y'] = SymbolicMod((y ** 2) * params.curve.parameters['r'], params.curve.prime)
        else:
            P_expr['Y'] = SymbolicMod(y, params.curve.prime)
    return P_expr


def algebraic_point_xy(params, x, y):
    # --- Projective point representation of P
    P_expr = {var: Mod(1, params.curve.prime) for var in params.curve.coordinate_model.variables}
    P_expr['X'] = SymbolicMod(x, params.curve.prime)
    P_expr['Y'] = SymbolicMod(y, params.curve.prime)
    if 'T' in params.curve.coordinate_model.variables:
        # Extended coordinates
        P_expr['T'] = SymbolicMod(x * y, params.curve.prime)
    return P_expr


def print_magma(params, basis_str, ring_els):
    # --- Printing in Magma format, which can be copied/pasted and executed at http://magma.maths.usyd.edu.au/calc/
    click.echo(f'fp := GF({params.curve.prime});')
    click.echo(f'P<{", ".join(str(elem) for elem in ring_els)}> := PolynomialRing(fp, {len(ring_els)});')
    click.echo(f'I := ideal< P | {basis_str} >;')
    click.echo('Groebner(I);')
    #click.echo('print I;')
    click.echo('V := Variety(I);')
    click.echo('print V;\n')

def print_sagemath(params, basis_str, ring_els):
    # --- Printing in sagemath format
    click.echo(f'from sage.all import *')
    basis_str = basis_str.replace('^', '**')
    ring_els_tup = tuple(map(str, ring_els))
    click.echo(f'fp = GF({params.curve.prime})')
    click.echo(f'pr = PolynomialRing(fp, {ring_els_tup}, {len(ring_els)})')
    ring_els_str = ','.join(list(map(str, ring_els)))
    click.echo(f'{ring_els_str} = pr.gens()')
    click.echo(f'I = pr.ideal({basis_str})')
    click.echo(f'B = pr.ideal(I.groebner_basis())')
    click.echo('V = B.variety()')
    click.echo('print(V)\n')

def print_raw(params, basis_str, ring_els):
    click.echo(',\n'.join(basis_str.split(', ')))

print_code = {
    'sagemath': print_sagemath, 
    'magma':print_magma,
    'raw':print_raw
}

def test_add_algebraic(formula, params, format):
    p = params.curve.prime
    field = FF(params.curve.prime)
    curve_params = {key: SymbolicMod(field(int(val)), p) for key, val in params.curve.parameters.items()}

    # --- Projective point representation of P
    P_expr = algebraic_point_xy(params, xP, yP)

    # --- Obtaining curve equation in terms of the affine point representation of P
    P = Point(params.curve.coordinate_model, **P_expr)
    if params.curve.coordinate_model.name != "inverted":
        Ploc = {**curve_params, "x": SymbolicMod(xP, p)}
        yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
        Peq = parse_expr(str(yPsquared - yP**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))
    else:
        Ploc = {**curve_params, "x": SymbolicMod(1/xP, p)}
        yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
        yPsquared_quot = fraction(sympify(str(yPsquared)))
        assert( sympify(str(yPsquared)) == yPsquared_quot[0]/yPsquared_quot[1] )
        Peq = parse_expr(str(yPsquared_quot[1] - yPsquared_quot[0] *(yP ** 2) ).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point representation of Q
    Q_expr = algebraic_point_xy(params, xQ, yQ)

    # --- Obtaining curve equation in terms of the affine point representation of Q
    Q = Point(params.curve.coordinate_model, **Q_expr)
    if params.curve.coordinate_model.name != "inverted":
        Qloc = {**curve_params, "x": SymbolicMod(xQ, p)}
        yQsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Qloc)
        Qeq = parse_expr(str(yQsquared - yQ**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))
    else:
        Qloc = {**curve_params, "x": SymbolicMod(1/xQ, p)}
        yQsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Qloc)
        yQsquared_quot = fraction(sympify(str(yQsquared)))
        assert( sympify(str(yQsquared)) == yQsquared_quot[0]/yQsquared_quot[1] )
        Qeq = parse_expr(str(yQsquared_quot[1] - yQsquared_quot[0] *(yQ ** 2) ).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point addition P + Q
    result = formula(params.curve.prime, P, Q, **curve_params)[0]

    # --- Our variety will be determined by the equations of P, Q, and P + Q
    basis = [ parse_expr(str(coordinate).replace('Mod(', '').replace(', %d)' % params.curve.prime, '')) for coordinate in result.coords.values() ]
    basis_str = str([Peq, Qeq] + basis)[1:-1]
    basis_str = basis_str.replace('**', '^')

    #print_magma(params, basis_str, [xP, xQ, yP, yQ])
    print_code[format](params, basis_str, [xP, xQ, yP, yQ])


def test_dbl_algebraic(formula, params, format):
    p = params.curve.prime
    field = FF(params.curve.prime)
    curve_params = {key: SymbolicMod(field(int(val)), p) for key, val in params.curve.parameters.items()}

    # --- Projective point representation of P
    P_expr = algebraic_point_xy(params, xP, yP)

    # --- Obtaining curve equation in terms of the affine point representation of P
    P = Point(params.curve.coordinate_model, **P_expr)
    if params.curve.coordinate_model.name != "inverted":
        Ploc = {**curve_params, "x": SymbolicMod(xP, p)}
        yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
        Peq = parse_expr(str(yPsquared - yP**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))
    else:
        Ploc = {**curve_params, "x": SymbolicMod(1/xP, p)}
        yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
        yPsquared_quot = fraction(sympify(str(yPsquared)))
        assert( sympify(str(yPsquared)) == yPsquared_quot[0]/yPsquared_quot[1] )
        Peq = parse_expr(str(yPsquared_quot[1] - yPsquared_quot[0] *(yP ** 2) ).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point addition [2]P
    result = formula(p, P, **curve_params)[0]

    # --- Our variety will be determined by the equations of P, and [2]P
    basis = [ parse_expr(str(coordinate).replace('Mod(', '').replace(', %d)' % params.curve.prime, '')) for coordinate in result.coords.values() ]
    basis_str = str([Peq] + basis)[1:-1]
    basis_str = basis_str.replace('**', '^')

    #print_magma(params, basis_str, [xP, xQ, yP, yQ])
    print_code[format](params, basis_str, [xP, xQ, yP, yQ])


def test_dadd_algebraic(formula, params, format):
    p = params.curve.prime
    field = FF(params.curve.prime)
    curve_params = {key: SymbolicMod(field(int(val)), p) for key, val in params.curve.parameters.items()}

    P_expr = algebraic_point(params, xP, yP)

    # --- Obtaining curve equation in terms of the affine point representation of P
    P = Point(params.curve.coordinate_model, **P_expr)
    Ploc = {**curve_params, "x": SymbolicMod(xP, p)}
    yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
    Peq = parse_expr(str(yPsquared - yP**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point representation of Q
    Q_expr = algebraic_point(params, xQ, yQ)

    # --- Obtaining curve equation in terms of the affine point representation of Q
    Q = Point(params.curve.coordinate_model, **Q_expr)
    Qloc = {**curve_params, "x": SymbolicMod(xQ, p)}
    yQsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Qloc)
    Qeq = parse_expr(str(yQsquared - yQ**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point representation of P - Q
    PQ_expr = algebraic_point(params, xPQ, yPQ)

    # --- Obtaining curve equation in terms of the affine point representation of PQ
    PQ = Point(params.curve.coordinate_model, **PQ_expr)
    PQloc = {**curve_params, "x": SymbolicMod(xPQ, p)}
    yPQsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), PQloc)
    PQeq = parse_expr(str(yPQsquared - yPQ**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Affine point representation of P and Q
    affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
    P_expr = {'x': SymbolicMod(xP, p), 'y': SymbolicMod(yP, p)}
    affine_P = Point(affine_model, **P_expr)
    Q_expr = {'x': SymbolicMod(xQ, p), 'y': SymbolicMod(yQ, p)}
    affine_Q = Point(affine_model, **Q_expr)

    # --- Obtaining curve equation in terms of the affine point representation of P - Q
    DIFF = params.curve.affine_add(affine_P, params.curve.affine_negate(affine_Q))
    # x(P-Q) = xPQ
    xPQ_quot = fraction(sympify(str(DIFF.x)))
    assert( sympify(str(DIFF.x)) == xPQ_quot[0]/xPQ_quot[1] )
    xPQeq = parse_expr(str(xPQ_quot[1]*xPQ - xPQ_quot[0]).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))
    # y(P-Q) = yPQ
    yPQ_quot = fraction(sympify(str(DIFF.y)))
    assert( sympify(str(DIFF.y)) == yPQ_quot[0]/yPQ_quot[1] )
    yPQeq = parse_expr(str(yPQ_quot[1]*yPQ - yPQ_quot[0]).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point addition P + Q
    result = formula(p, PQ, P, Q, **curve_params)[0]

    # --- Our variety will be determined by the equations of P, Q, PQ, xPQeq, yPQeq, and P + Q
    basis = [ parse_expr(str(coordinate).replace('Mod(', '').replace(', %d)' % params.curve.prime, '')) for coordinate in result.coords.values() ]
    basis_str = str([Peq, Qeq, PQeq] + [xPQeq, yPQeq] + basis)[1:-1]
    basis_str = basis_str.replace('**', '^')

    #print_magma(params, basis_str, [xP, xQ, xPQ, yP, yQ, yPQ])
    print_code[format](params, basis_str, [xP, xQ, xPQ, yP, yQ, yPQ])


def test_ladd_algebraic(formula, params, format):
    p = params.curve.prime
    field = FF(params.curve.prime)
    curve_params = {key: SymbolicMod(field(int(val)), p) for key, val in params.curve.parameters.items()}

    # --- Projective point representation of P
    P_expr = algebraic_point(params, xP, yP)

    # --- Obtaining curve equation in terms of the affine point representation of P
    P = Point(params.curve.coordinate_model, **P_expr)
    Ploc = {**curve_params, "x": SymbolicMod(xP, p)}
    yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
    Peq = parse_expr(str(yPsquared - yP**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point representation of Q
    Q_expr = algebraic_point(params, xQ, yQ)

    # --- Obtaining curve equation in terms of the affine point representation of Q
    Q = Point(params.curve.coordinate_model, **Q_expr)
    Qloc = {**curve_params, "x": SymbolicMod(xQ, p)}
    yQsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Qloc)
    Qeq = parse_expr(str(yQsquared - yQ**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point representation of P - Q
    PQ_expr = algebraic_point(params, xPQ, yPQ)

    # --- Obtaining curve equation in terms of the affine point representation of PQ
    PQ = Point(params.curve.coordinate_model, **PQ_expr)
    PQloc = {**curve_params, "x": SymbolicMod(xPQ, p)}
    yPQsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), PQloc)
    PQeq = parse_expr(str(yPQsquared - yPQ**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Affine point representation of P and Q
    affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
    P_expr = {'x': SymbolicMod(xP, p), 'y': SymbolicMod(yP, p)}
    affine_P = Point(affine_model, **P_expr)
    Q_expr = {'x': SymbolicMod(xQ, p), 'y': SymbolicMod(yQ, p)}
    affine_Q = Point(affine_model, **Q_expr)

    # --- Obtaining curve equation in terms of the affine point representation of P - Q
    DIFF = params.curve.affine_add(affine_P, params.curve.affine_negate(affine_Q))
    # x(P-Q) = xPQ
    xPQ_quot = fraction(sympify(str(DIFF.x)))
    assert( sympify(str(DIFF.x)) == xPQ_quot[0]/xPQ_quot[1] )
    xPQeq = parse_expr(str(xPQ_quot[1]*xPQ - xPQ_quot[0]).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))
    # y(P-Q) = yPQ
    yPQ_quot = fraction(sympify(str(DIFF.y)))
    assert( sympify(str(DIFF.y)) == yPQ_quot[0]/yPQ_quot[1] )
    yPQeq = parse_expr(str(yPQ_quot[1]*yPQ - yPQ_quot[0]).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point addition P + Q
    result1, result2 = formula(p, PQ, P, Q, **curve_params)[:2]

    # --- Our variety will be determined by the equations of P, Q, PQ, xPQeq, yPQeq, P + Q, and [2]P
    basis1 = [ parse_expr(str(coordinate).replace('Mod(', '').replace(', %d)' % params.curve.prime, '')) for coordinate in result1.coords.values() ]
    basis2 = [ parse_expr(str(coordinate).replace('Mod(', '').replace(', %d)' % params.curve.prime, '')) for coordinate in result2.coords.values() ]
    basis_str = str([Peq, Qeq, PQeq] + [xPQeq, yPQeq] + basis1 + basis2)[1:-1]
    basis_str = basis_str.replace('**', '^')

    #print_magma(params, basis_str, [xP, xQ, xPQ, yP, yQ, yPQ])
    print_code[format](params, basis_str, [xP, xQ, xPQ, yP, yQ, yPQ])


def test_ddbl_algebraic(formula, params, format):
    p = params.curve.prime
    field = FF(params.curve.prime)
    curve_params = {key: SymbolicMod(field(int(val)), p) for key, val in params.curve.parameters.items()}

    # --- Projective point representation of P
    P_expr = algebraic_point(params, xP, yP)

    # --- Obtaining curve equation in terms of the affine point representation of P
    P = Point(params.curve.coordinate_model, **P_expr)
    Ploc = {**curve_params, "x": SymbolicMod(xP, p)}
    yPsquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), Ploc)
    Peq = parse_expr(str(yPsquared - yP**2).replace('Mod(', '').replace(', %d)' % params.curve.prime, ''))

    # --- Projective point addition [2]P
    result = formula(p, P, **curve_params)[0]

    # --- Our variety will be determined by the equations of P, and [2]P
    basis = [ parse_expr(str(coordinate).replace('Mod(', '').replace(', %d)' % params.curve.prime, '')) for coordinate in result.coords.values() ]
    basis_str = str([Peq] + basis)[1:-1]
    basis_str = basis_str.replace('**', '^')

    #print_magma(params, basis_str, [xP, xQ, xPQ, yP, yQ, yPQ])
    print_code[format](params, basis_str, [xP, xQ, xPQ, yP, yQ, yPQ])


@click.command()
@click.option('--std', default=None, type=str, help="File with extension .json containing the curves to be used (same format as in std-curve).")
@click.option('--efd', default=None, type=str, help="Formula to be used (same name as in the efd collection).")
@click.option('--format', default=None, type=click.Choice(("sagemath", "magma", "raw")), help="sagemath or magma or raw output format).")
def main(std, efd, format):
    original_stdout = sys.stdout # Save a reference to the original standard output
    # Model names
    model_names = {
            "Weierstrass": "Short Weierstrass curve models",
            "Montgomery": "Montgomery curve models",
            "TwistedEdwards": "Twisted Edwards curve models",
            "Edwards": "Edwards curve models",
            }

    # Curve Models
    model_list = {
            "Weierstrass": ShortWeierstrassModel(),
            "Montgomery": MontgomeryModel(),
            "TwistedEdwards": TwistedEdwardsModel(),
            "Edwards": EdwardsModel()
            }

    # Formulas to be used
    tests = {
        "add": test_add_algebraic,
        "dbl": test_dbl_algebraic
    }
    tests_kummer = {
        "dadd": test_dadd_algebraic,
        "ladd": test_ladd_algebraic,
        "dbl": test_ddbl_algebraic
    }

    curves, length = json2curves(std)  # Reading the file with extension .json that contains the list of curves to be used

    # ===========================================================================================================================
    for i in range(0, length, 1):

        FORMULA_NOT_FOUND = True    # This is order to ensure the dseired formula belongs to the EFD
        message = "PROCESSING: " + curves['curves'][i]['desc']
        print("-" * len(message) + "\n# " + message)
        with open('%s-%s.log' % (curves['curves'][i]['name'], efd), 'w') as f:
            sys.stdout = f  # Change the standard output to the file we created.
            curve_name = curves['curves'][i]['name']
            model = model_list[curves['curves'][i]['form']]
            assert( model in model_list.values() )
            print("# " + "=" * len(curves['curves'][i]['desc']) + "\n# " + curves['curves'][i]['desc'])
            for coord_name, coords in model.coordinates.items():

                try:
                    curve = get_curve(curves, i, coord_name)    # now infty = False (projective representation of the point at infinty)
                except ValueError as e:
                    click.echo(f"# {coord_name} {curve_name} {e}")
                    continue

                #print(efd, coords.formulas.keys())
                for formula_name, formula in coords.formulas.items():
                    if efd == formula_name:
                        FORMULA_NOT_FOUND = False
                        if (formula.shortname in tests) or (formula.shortname in tests_kummer):
                            print("# " + "-" * (len(coord_name + curve_name + formula_name) + 2))
                            click.echo(f"# {coord_name} {curve_name} {formula_name}\n")
                            try:
                                with TemporaryConfig() as cfg:
                                    cfg.ec.mod_implementation = "python"
                                    if len(coords.variables) > 2:
                                        # Rational points: point addition and doubling
                                        tests[formula.shortname](formula, curve, format)
                                    else:
                                        # Kummer line: differential point additions and doubling
                                        assert( len(coords.variables) == 2)
                                        tests_kummer[formula.shortname](formula, curve, format)
                            except NameError as e:
                                click.echo(f"# {coord_name} {curve_name} {formula_name} Bad formula! {e}")
                                continue

            sys.stdout = original_stdout  # Reset the standard output to its original value

        if FORMULA_NOT_FOUND:
            click.echo(f"# WARNING: {efd} is not a valid formula for {model_names[curves['curves'][i]['form']]}")


if __name__ == "__main__":
    main()

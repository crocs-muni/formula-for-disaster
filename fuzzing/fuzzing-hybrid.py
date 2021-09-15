import subprocess
import sys

import click

def to_sagemath_format(shortname: str, coords: str, operation: str, formula: str, local_pyecsca_path=None, verbose=False):
    if local_pyecsca_path is None:
        base_search_point = "$HOME"
    else:
        base_search_point = local_pyecsca_path
    path_str = f"find {base_search_point} -type f -name {formula}.op3 | grep pyecsca/pyecsca/ec/efd/{shortname}/{coords}/{operation}/{formula}.op3 "
    path_formula = subprocess.check_output("""%s""" % path_str, shell=True)
    path_formula = str(path_formula).replace('b\'', '').replace('\\n\'', '')

    click.echo('formula = {}')
    with open(path_formula) as formula:
        for line in formula:
            click.echo(f"{line.strip().replace('+', ' + ').replace('-', ' - ').replace('*', ' * ').replace('/', ' / ').replace('^', ' ** ')}")
            left_side, right_side = line.strip().replace('^', ' ** ').split(" = ")
            click.echo(f"formula['{left_side}'] = {left_side}")

    if verbose:
        click.echo('for key, value in formula.items():')
        click.echo("\tprint(f\'{key} = {value}\')")

    return 0


def export_formula(model, coords, efd, local_pyecsca_path=None):
    original_stdout = sys.stdout  # Save a reference to the original standard output
    name = 'formulae/addition/%s-%s-%s-addition.py' % (model, coords, efd)
    with open(name.replace('-', '_'), 'w') as f:

        sys.stdout = f  # Change the standard output to the file we created.
        click.echo(f'from sage.all import PolynomialRing, ZZ\n')

        if model == 'shortw':
            click.echo(f'pr = PolynomialRing(ZZ, (\'a\', \'b\', \'X1\', \'X2\', \'Y1\', \'Y2\', \'Z1\', \'Z2\'), 8)')
            click.echo(f'a, b, X1, X2, Y1, Y2, Z1, Z2 = pr.gens()')
            click.echo(f'ZZ1, ZZZ1, ZZ2, ZZZ2 = Z1**2, Z1**3, Z2**2, Z2**3')
            click.echo(f'b2 = 2 * b')
            click.echo(f'b4 = 4 * b')
            click.echo(f'half = 1 / 2')
        elif model == 'twisted':
            click.echo(f'pr = PolynomialRing(ZZ, (\'a\', \'d\', \'X1\', \'X2\', \'Y1\', \'Y2\', \'Z1\', \'Z2\'), 8)')
            click.echo(f'a, d, X1, X2, Y1, Y2, Z1, Z2 = pr.gens()')
            click.echo(f'k, d2 = 2 * d, 2 * d')
        elif model == 'edwards':
            click.echo(f'pr = PolynomialRing(ZZ, (\'c\', \'d\', \'i\', \'X1\', \'X2\', \'Y1\', \'Y2\', \'Z1\', \'Z2\'), 9)')
            click.echo(f'c, d, i, X1, X2, Y1, Y2, Z1, Z2 = pr.gens()')
            click.echo(f'ccd = c * c * d')
            click.echo(f'ccd2 = 2 * c * c * d')
            click.echo(f'c2 = 2 * c')
            click.echo(f'cc4 = 4 * c * c')
            click.echo(f'k = 1 / c')
        if 'extended' in coords:
            click.echo('T1 = X1 * Y1')
            click.echo('T2 = X2 * Y2')

        output = to_sagemath_format(model, coords, 'addition', efd, local_pyecsca_path)

    sys.stdout = original_stdout  # Reset the standard output to its original value
    return output


@click.command()
@click.option('--model', default=None, type=click.Choice(("shortw", "montgom", "twisted", "edwards")),
              help="Elliptic curve model.")
@click.option('--coords', default=None, type=str, help="Coordinate point representation.")
@click.option('--efd', default=None, type=str, help="Formula to be used (same name as in the efd collection).")
def main(model, coords, efd):
    return export_formula(model, coords, efd)


if __name__ == '__main__':
    assert (main() == 0)

# Fuzzing

Brute force search of exceptional points on elliptic curves in _Short Weierstrass_, _Montgomery_, _Edwards_, and _Twisted Edwards_ forms.

## Requirements

The toolkit [pyecsca](https://github.com/J08nY/pyecsca) (**Py**thon **E**lliptic **C**urve cryptography **S**ide-**C**hannel **A**nalysis), which
uses the [Explicit-Formulas Database](https://www.hyperelliptic.org/EFD/index.html) (EFD) by Daniel J. Bernstein and Tanja Lange.

### Brute force on all rational points
```bash
# Looking for all the possible exceptional points for each formula in the EFD on {5,8}-bits elliptic curves 
python3 fuzzing.py
```

The script `fuzzing.py` takes always as input the files `data/shortw.json`, data/motgom.json, data/twisted.json, and `data/edwards.json`.

### Brute force on a portion of the points

Both scripts `fuzzing.py` and `fuzzing-pruned.py` allows to work with any arbitrary elliptic curve instance.
However, `fuzzing-pruned.py` was designed to work for  elliptic curves of size `q`, `2q` and `4q` being `q` a _**large**_ prime number.

```bash
python3 fuzzing-pruned.py --help
Usage: fuzzing-pruned.py [OPTIONS]

Options:
  --std TEXT     File with extension .json containing the curves to be used
                 (same format as in std-curve).
  --efd TEXT     Formula to be used (same name as in the efd collection).
  --num INTEGER  Number of point addition to be performed.
  --help         Show this message and exit.

# Examples
python3 fuzzing-pruned.py --std data/shortw/shortw-8bits.json --efd add-2016-rcb --num 150
python3 fuzzing-pruned.py --std data/montgom/montgom-8bits.json --efd dadd-1987-m --num 75
python3 fuzzing-pruned.py --std data/twisted/twisted-8bits.json --efd madd-2008-hwcd --num 200
python3 fuzzing-pruned.py --std data/edwards/edwards-8bits.json --efd dadd-2006-g --num 50
python3 fuzzing-pruned.py --std data/edwards/edwards-8bits.json --efd madd-2007-bl-2 --num 250
python3 fuzzing-pruned.py --std data/edwards/edwards-8bits.json --efd add-2007-bl-4 --num 16
```

## Algebraic equations from a given formula
The script `fuzzing-algebraic.py` extends _pyecsca_ to handle the use of `sympy` polynomials (symbols) as coordinates.

```bash
# Output in Magma-code format
python3 fuzzing-algebraic.py --std data/shortw/shortw-8bits.json --efd add-2016-rcb --format magma
# Output in Sagemath-code format
python3 fuzzing-algebraic.py --std data/shortw/shortw-8bits.json --efd add-2016-rcb --format sagemath
# Next, sagemath can be runned in console mode as
sage -python shortw-prime-8bits-b-has-sqrt-add-2016-rcb.log

# Now, we can use another script (hybrid)
python3 fuzzing-hybrid.py --model twisted --coords extended-1 --efd add-2008-hwcd
python3 fuzzing-hybrid.py --model shortw --coords projective --efd add-2007-bl
python3 fuzzing-hybrid.py --model shortw --coords projective --efd add-2002-bj
python3 fuzzing-hybrid.py --model shortw --coords projective --efd add-2016-rcb
# For looking the list of equations, you can do it as in fuzzing-algebraic.py examples
sage -python twisted-extended-add-2008-hwcd-addition.log
sage -python shortw-projective-add-2007-bl-addition.log
sage -python shortw-projective-add-2002-bj-addition.log
sage -python shortw-projective-add-2016-rcb-addition.log

# ---
sage -python
> from formulae.twisted_extended_1_add_2008_hwcd_addition import *
> formula
{'A': X1*X2, 'B': Y1*Y2, 't0': d*X2*Y2, 'C': d*X1*X2*Y1*Y2, 'D': 1, 't1': X1 + Y1, 't2': X2 + Y2, 't3': X1*X2 + X2*Y1 + X1*Y2 + Y1*Y2, 't4': X2*Y1 + X1*Y2 + Y1*Y2, 'E': X2*Y1 + X1*Y2, 'F': -d*X1*X2*Y1*Y2 + 1, 'G': d*X1*X2*Y1*Y2 + 1, 't5': a*X1*X2, 'H': -a*X1*X2 + Y1*Y2, 'X3': -d*X1*X2^2*Y1^2*Y2 - d*X1^2*X2*Y1*Y2^2 + X2*Y1 + X1*Y2, 'Y3': -a*d*X1^2*X2^2*Y1*Y2 + d*X1*X2*Y1^2*Y2^2 - a*X1*X2 + Y1*Y2, 'T3': -a*X1*X2^2*Y1 - a*X1^2*X2*Y2 + X2*Y1^2*Y2 + X1*Y1*Y2^2, 'Z3': -d^2*X1^2*X2^2*Y1^2*Y2^2 + 1}
```

## Adding new curve instances
The curve instances must be stored as a file with extension `.json` and structured as follows

```json
# Short Weierstrass model
{
        "name": "shortw-test",
        "desc": "Elliptic curve instances in short Weierstrass form over small primes.",
        "curves": [
{
        "form": "Weierstrass",
        "name": "small-5bits-prime-curve-b-has-sqrt",
        "category": "shortw-test",
        "desc": "A small 5-bit prime curve with points of the form (0, y).",
        "field": {
                "type": "Prime",
                "p": "0x1F",
                "bits":  5
        },
        "params": {
                "a": { "raw": "0x1A" },
                "b": { "raw": "0x1C" }
        },
        "generator": {
                "x": { "raw": "0x3" },
                "y": { "raw": "0x3" }
        },
        "order": "0x1F",
        "cofactor": "0x1"
}
        ]
}

# Montgomery model
{
        "name": "montgom-test",
        "desc": "Elliptic curve instances in Montgomery form over small primes.",
        "curves": [
{
        "form": "Montgomery",
        "name": "small-5bits-montgomery-curve-b1",
        "category": "montgom-test",
        "desc": "A small 5-bit Montgomery curve with with b = 1 and a different from 0 and 6.",
        "field": {
                "type": "Prime",
                "p": "0x13",
                "bits":  5
        },
        "params": {
                "a": { "raw": "0x8" },
                "b": { "raw": "0x1" }
        },
        "generator": {
                "x": { "raw": "0xC" },
                "y": { "raw": "0x2" }
        },
        "order": "0x1C",
        "cofactor": "0x1"
}
        ]
}

# Twisted Edwards model
{
        "name": "twisted-test",
        "desc": "Elliptic curve instances in Twisted Edwards form over small primes.",
        "curves": [
{
        "form": "TwistedEdwards",
        "name": "small-5bits-twisted-edwards-curve-a1",
        "category": "twisted-test",
        "desc": "A small 5-bit Twisted Edwards curve with with a = 1.",
        "field": {
                "type": "Prime",
                "p": "0x13",
                "bits":  5
        },
        "params": {
                "a": { "raw": "0x1" },
                "d": { "raw": "0x8" }
        },
        "generator": {
                "x": { "raw": "0x10" },
                "y": { "raw": "0xE" }
        },
        "order": "0x70",
        "cofactor": "0x1"
}
        ]
}

# Edwards
{
        "name": "edwards-test",
        "desc": "Elliptic curve instances in Edwards form over small primes.",
        "curves": [
{
        "form": "Edwards",
        "name": "small-5bits-edwards-curve-c1-d-has-sqrt",
        "category": "edwards-test",
        "desc": "A small 5-bit Edwards curve with with c = 1 and d has sqrt.",
        "field": {
                "type": "Prime",
                "p": "0x1F",
                "bits":  5
        },
        "params": {
                "c": { "raw": "0x1" },
                "d": { "raw": "0x8" }
        },
        "generator": {
                "x": { "raw": "0xB" },
                "y": { "raw": "0x13" }
        },
        "order": "0x50",
        "cofactor": "0x1"
}
        ]
}
```
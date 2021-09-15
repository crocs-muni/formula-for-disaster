# A formula for disaster: a unified approach to elliptic curve special-point-based attacks

This repository contains data and scripts used in the [A formula for disaster: a unified approach to elliptic curve special-point-based attacks](https://crocs.fi.muni.cz/public/papers/formulas_asiacrypt21) paper at ASIACRYPT 2021.

## Abstract

> The Refined Power Analysis, Zero-Value Point, and Exceptional Procedure 
> attacks introduced side-channel attack techniques against
> specific cases of elliptic curve cryptography. The three attacks recover
> bits of a static ECDH key adaptively, collecting information on whether
> a certain multiple of the input point was computed. We unify and generalize
> these attacks in a common framework, and solve the corresponding
> problem for a broader class of inputs. We also introduce a version of
> the attack against windowed scalar multiplication methods, recovering
> the full scalar instead of just a part of it. Finally, we systematically
> analyze elliptic curve point addition formulas from the Explicit-Formulas
> Database, classify all non-trivial exceptional points, and find them
> in new formulas. These results indicate the usefulness of our tooling for
> unrolling formulas and finding special points, potentially of independent
> research interest.


## Contents

 * `unrolling/` -> Scripts/notebooks and data of unrolled formulas from the Explicit-Formulas Database.
 * `epa/` -> Scripts and data related to the Exceptional Procedure Attack.
 * `rpa/` -> Scripts and data related to the Refined Power Analysis attack.
 * `zvp/` -> Scripts and data related to the Zero-Value Point attack.
 * `fuzzing/` -> Scripts and data related to the [fuzzing](./fuzzing/README.md) search (brute force search).

## Requirements

The notebooks are Jupyter notebooks and as such require Jupyter to run.
Some notebooks also require a SageMath kernel.
The scripts and notebooks use the [pyecsca](https://neuromancer.sk/pyecsca/)
toolkit and an export of the Explicit-Formulas Database available at [efd](https://github.com/J08nY/efd),
which is also a part of pyecsca.


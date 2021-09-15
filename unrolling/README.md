# Formula unrolling

This directory contains code for unrolling formulas from the EFD as well as the unrolled
formulas themselves.

## Contents

 * `add_unroll.ipynb` -> Addition formula unrolling script.
 * `dbl_unroll.ipynb` -> Doubling formula unrolling script.
 * `ladd_unroll.ipynb` -> Ladder and DifferentialAddition formula unrolling script.
 * `unroll/` -> Unrolled formulas with all intermediate variables. The value assigned to the
   variables is the one that was left in the variable at the end of the execution of the formula.
   Each variable might have been reused and these intermediates are thus not captured in
   the output (but pyecsca tracks these and they can be extracted).
 * `unroll_out/` -> Unrolled formulas with only the outputs.

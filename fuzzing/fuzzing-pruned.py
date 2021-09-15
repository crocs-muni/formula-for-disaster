from pyecsca.ec.params import load_params, _create_params
from pyecsca.ec.model import ShortWeierstrassModel, MontgomeryModel, EdwardsModel, TwistedEdwardsModel
from pyecsca.ec.formula import AdditionFormula, DoublingFormula
from pyecsca.ec.point import InfinityPoint, Point

from pyecsca.ec.coordinates import AffineCoordinateModel
from pyecsca.ec.mod import Mod

from pkg_resources import resource_filename
import click

import json
import yaml

# ===========================================================================================================================
# This function maps from projective to affine coordinates (even if the point at infinity doesn't have affine representation)
def proj_to_affine(point : Point, params):

    # ---------------------------------------------------------------------------------------------------------
    # Scaling the point at infinity must be handled in a different way (at least for shortw and montgom models)
    tmp = { k : point.coords[k] * params.curve.neutral.coords[k] for k in point.coords.keys() }
    infty = Point(params.curve.coordinate_model, **tmp)

    # --------------------------------------------------------------
    # We are mapping the degenerate zero-vector state into (0,0) ... 
    # This should be modified (a new class, posibly named DegeneratedPoint, must be included in pyecsca)
    affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
    degstate = { 'x':Mod(0, params.curve.prime), 'y':Mod(0, params.curve.prime)}
    degstate = Point(affine_model, **degstate)

    try:
        # If to_affine() arises an error, then there is a "division" by 0
        affine_point = point.to_affine()
        return affine_point

    except:

        if (max([ int(v) for v in point.coords.values()]) > 0):
            # Point at infinity
            return params.curve.neutral
        else:
            # Degenerated zero-vector
            return degstate

# ===========================================================================================================================
# json2curves read a list of curves in a given json file
def json2curves(filename: str):
    f = open(filename,)
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
    for c,v in point.coords.items():
        strpoint += c + ': 0x%0*x, ' % (b, int(v)) 
    return strpoint[:-2] + '}' 

# ===========================================================================================================================
# Next function computes a 4-order point
def fullorder_generator(P : Point, params):

    # Recall, For edwards and twisted models, the generator is assumed to have prime order, thus we require a 4-order point
    if params.curve.model.shortname == 'twisted' or params.curve.model.shortname == 'edwards':

        # If the input curve is assumed to have h-order with h > 4, then the next lines must be modified
        affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
        if params.curve.model.shortname == 'twisted':
            sqrt_a = params.curve.parameters['a'].sqrt()
            invsqrt_a = Mod(1, params.curve.prime) / sqrt_a
            P4_raw = { 'x':int(invsqrt_a), 'y':Mod(0, params.curve.prime)}                        # Affine representation of one 4-order point (its negative is (-1/sqrt(a),0))

        else:
            P4_raw = { 'x':int(params.curve.parameters['c']), 'y':Mod(0, params.curve.prime)}     # Affine representation of one 4-order point (its negative is (-c,0))

        P4_affine = Point(affine_model, **P4_raw)
        assert(params.curve.is_on_curve(P4_affine))
        P_full = params.curve.affine_add(P, P4_affine)
        assert(params.curve.is_on_curve(P_full))
        return P_full
    else:

        # For shortw and montgom, we have a full order generator. No extra point is required
        return P

# ===========================================================================================================================
# For inverted coordinates, the mapping of the torsion-4 points is seperately handled
def mapping_special_points(point : Point, params):

    if params.curve.coordinate_model.name == "inverted":

        if point.x == 0:
            # torsion-2 point case
            result = { 'X':point.y, 'Y':Mod(0, params.curve.prime), 'Z':Mod(0, params.curve.prime)}
            return Point(params.curve.coordinate_model, **result)

        elif point.y == 0:
            # 4-order point case
            result = { 'X':Mod(0, params.curve.prime), 'Y':-point.x, 'Z':Mod(0, params.curve.prime)}
            return Point(params.curve.coordinate_model, **result)

        else: 
            # General case
            return point.to_model(params.curve.coordinate_model, params.curve)
    else:
        # General case
        return point.to_model(params.curve.coordinate_model, params.curve)

# ===========================================================================================================================
# Testing point addition (P + Q, P + Q, P + O, and 0 + P)
def test_add_formula(formula, params, num):
    # Kummer-line has elements in the projective space P^2
    assert( len(params.curve.coordinate_model.variables) > 2 )
    affine_gen = params.generator.to_affine()

    ok = True
    ok_counterexample = []
    unified = True
    unified_counterexample = []
    complete = True #not isinstance(params.curve.neutral, InfinityPoint)
    complete_counterexample = []

    # Let's precompute the cyclic subgroup
    points = []
    assert(params.curve.affine_multiply(affine_gen, params.order - 1) == params.curve.affine_negate(affine_gen))
    points.append(affine_gen)
    assert(params.curve.is_on_curve(points[0]))
    points.append(params.curve.affine_double(points[0]))                    # [2]g
    assert(params.curve.is_on_curve(points[1]))
    for i in range(2, min(num, params.order) - 1, 1):
        points.append(params.curve.affine_add(points[0], points[i - 1]))    # [i + 1]g
        assert(params.curve.is_on_curve(points[i]))

    digits = len(str(int(min(num, params.order))))
    for i in range(1, min(num, params.order)):
        #affine_P = params.curve.affine_multiply(affine_gen, i)
        affine_P = points[i - 1]
        P = mapping_special_points(affine_P, params)

        # Next IF should be modified for catching possible errors in the inverted coordinates (catching special points in inverted)
        if params.curve.coordinate_model.name != "inverted":
            assert(params.curve.is_on_curve(P))
        
        # Do P + Q for all P,Q on the curve, s.t. P != Q and P + Q != infty then compare with affine add
        for j in range(1, min(num, params.order)):
            # TODO: Test also the case when P + Q = infty
            if i == j or i + j == params.order:
                continue
            #affine_Q = params.curve.affine_multiply(affine_gen, j)
            affine_Q = points[j - 1]
            Q = mapping_special_points(affine_Q, params)

            # Next IF should be modified for catching possible errors in the inverted coordinates (catching special points in inverted)
            if params.curve.coordinate_model.name != "inverted":
                assert(params.curve.is_on_curve(Q))

            try:
                # If the formulae are not exception-free, then the next to point additions will raise an Unsatisfied Assumption Error
                result = formula(params.curve.prime, P, Q, **params.curve.parameters)[0]
                affine_result = proj_to_affine(result, params)
                affine_expected = params.curve.affine_add(affine_P, affine_Q)
                if affine_result != affine_expected:
                    ok = False
                    ok_counterexample.append(\
                            point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.order)[2:])) \
                            + ', and ' + \
                            point2str(f'[{j:{digits}d}]P', affine_Q, len(hex(params.order)[2:])) \
                            )
            except:
                # This branch corresponds with an Unsatisfied Assumption Error (and thus, it is not an exception-free formula)
                ok = False
                ok_counterexample.append(\
                        point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.order)[2:])) \
                        + ', and ' + \
                        point2str(f'[{j:{digits}d}]P', affine_Q, len(hex(params.order)[2:])) \
                        )

        try:
            # Do P + P and compare with affine doubling of P. If this works, the formulas are unified
            # If the formulae are not exception-free, then the next to point additions will raise an Unsatisfied Assumption Error
            result = formula(params.curve.prime, P, P, **params.curve.parameters)[0]
            #affine_result = result.to_affine()
            affine_result = proj_to_affine(result, params)

            # =======================================================================================
            # THE AFFINE IMPLEMENTATION OF THE AFFINE POINT DOUBLING FAILS WHEN P IS AN 2-ORDER POINT
            if params.curve.affine_negate(affine_P) == affine_P:
                if params.curve.model.shortname == 'twisted' or params.curve.model.shortname == 'edwards':
                    affine_expected = params.curve.affine_double(affine_P)
                else:
                    affine_expected = params.curve.neutral
            else:
                affine_expected = params.curve.affine_double(affine_P)
            # ========================================================================================

            assert(params.curve.is_on_curve(affine_result) or\
                    (affine_result == params.curve.neutral) or (max([ int(v) for v in result.coords.values()]) == 0) or\
                    params.curve.model.shortname == "twisted" or params.curve.model.shortname == "edwards")

            assert(params.curve.is_on_curve(affine_expected) or (affine_expected == params.curve.neutral))

            if affine_result != affine_expected:
                unified = False
                unified_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )

        except:
            # This branch corresponds with an Unsatisfied Assumption Error (and thus, it is not an exception-free formula)
            unified = False
            unified_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )

        try:
            # Do P + infty, and infty + P, if infty is representable in coordinate system. If this works, the formulas are complete
            # If the formulae are not exception-free, then the next to point additions will raise an Unsatisfied Assumption Error
            result2 = formula(params.curve.prime, params.curve.neutral, P, **params.curve.parameters)[0]
            result1 = formula(params.curve.prime, P, params.curve.neutral, **params.curve.parameters)[0]
            #affine_result1 = result1.to_affine()
            #affine_result2 = result2.to_affine()
            affine_result1 = proj_to_affine(result1, params)
            affine_result2 = proj_to_affine(result2, params)
            if affine_result1 != affine_P or affine_result2 != affine_P:
                complete = False
                complete_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )
        except:
            # This branch corresponds with an Unsatisfied Assumption Error (and thus, it is not an exception-free formula)
            complete = False
            complete_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )


    return {\
            'P + Q' : { 'Are there failed computations?' : not ok,
                'Failures' : ok_counterexample, 
                'Number of failures' : len(ok_counterexample)},\
            'P + P' : { 'Are there failed computations?' : not unified, 
                'Failures' : unified_counterexample if len(unified_counterexample) < (min(num, params.order) - 1) else "All the points (not including the point at infinity)", 
                'Number of failures' : len(unified_counterexample) },\
            'P + infty and infty + P' : { 'Are there failed computations?' : not complete, 
                'Failures' : complete_counterexample if len(complete_counterexample) < (min(num, params.order) - 1) else "All the points (not including the point at infinity)", 
                'Number of failures' : len(complete_counterexample) }\
            }


# ===========================================================================================================================
# Testing point doubling
def test_dbl_formula(formula, params, num):
    # Kummer-line has elements in the projective space P^2
    assert( len(params.curve.coordinate_model.variables) > 2 )
    affine_gen = params.generator.to_affine()

    ok = True
    ok_counterexample = []
    # Let's precompute the cyclic subgroup
    points = []
    assert(params.curve.affine_multiply(affine_gen, params.order - 1) == params.curve.affine_negate(affine_gen))
    points.append(affine_gen)
    assert(params.curve.is_on_curve(points[0]))
    points.append(params.curve.affine_double(points[0]))                    # [2]g
    assert(params.curve.is_on_curve(points[1]))
    for i in range(2, min(num, params.order) - 1, 1):
        points.append(params.curve.affine_add(points[0], points[i - 1]))    # [i + 1]g
        assert(params.curve.is_on_curve(points[i]))

    digits = len(str(int(min(num, params.order))))
    for i in range(1, min(num, params.order)):
        #affine_P = params.curve.affine_multiply(affine_gen, i)
        affine_P = points[i - 1]
        P = mapping_special_points(affine_P, params)

        try:
            result = formula(params.curve.prime, P, **params.curve.parameters)[0]
            #affine_result = result.to_affine()
            affine_result = proj_to_affine(result, params)
            #affine_expected = params.curve.affine_double(affine_P)
            if params.curve.affine_negate(affine_P) == affine_P:
                if params.curve.model.shortname == 'twisted' or params.curve.model.shortname == 'edwards':
                    affine_expected = params.curve.affine_double(affine_P)
                else:
                    affine_expected = params.curve.neutral
            else:
                affine_expected = params.curve.affine_double(affine_P)

            if affine_result != affine_expected:
                ok = False
                ok_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )

        except:
            # This branch corresponds with an Unsatisfied Assumption Error (and thus, it is not an exception-free formula)
            ok = False
            ok_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )


    return { '[2]P' : { 'Are there failed computations?' : not ok, 
        'Failures' : ok_counterexample if len(ok_counterexample) < (min(num, params.order) - 1) else "All the points (not including the point at infinity)", 
        'Number of failures' : len(ok_counterexample) } }

# ===========================================================================================================================
# Mapping from the Kummer line to affine points
# Remark: the kummer line corresponds with E/{+1,-1}, which implies the affine points P and -P have the same Kummer representation
# The case of yzsquared coordinates is more special (the mapping is 1-to-4), thus the tests must take it in account
def kummer_to_affine(P, params):

    assert( params.curve.model.shortname != "twisted" ) # The current op3 doesn't have Kummer line arithmetic for the twisted edwards
    values_P = [ int(v) for v in P.coords.values() ]
    assert( len(values_P) == 2)
    result = {}
    if max(values_P) > 0:
        # Recall, the point at infinity is (1 : 0)
        if values_P[0] != 0 and values_P[1] == 0:
            return params.curve.neutral

        # Next, we must to compute X/Z or Y/Z  (no matter which model curve is used, this "division" is required)
        Z = Mod(1, params.curve.prime) / values_P[1]
        result['x'] = Z * values_P[0]
        if params.curve.model.shortname == "edwards":
            assert(params.curve.coordinate_model.name == 'yz' or params.curve.coordinate_model.name == 'yzsquared')
            # Recall, d = r^2
            r_inv = params.curve.parameters['r'].inverse()                          # 1/r
            if params.curve.coordinate_model.name == "yz":
                result['x'] = result['x'] * r_inv                                   # Recall, r*y = Y/Z
            if params.curve.coordinate_model.name == "yzsquared":
                result['x'] = result['x'] * r_inv                                   # Recall, r*y^2 = Y/Z
                result['x'] = result['x'].sqrt()

        # Evaluating x(P) in the curve equation and computing y(P) by using sqrt computations inf GF(p)
        loc = {**params.curve.parameters, "x": result['x']}
        ysquared = eval(compile(params.curve.model.ysquared, "", mode="eval"), loc)
        result['y'] = ysquared.sqrt()

        if params.curve.model.shortname == "edwards":
            # Edwards curves have particularities, we just need to swap x <-> y
            result = { 'x':result['y'], 'y':result['x'] }

        affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
        return Point(affine_model, **result)
    else:
        assert(values_P == [0,0])
        result['x'] = Mod(0, params.curve.prime)
        result['y'] = Mod(0, params.curve.prime)
        # This case implies P is the degenerated zero-vector
        affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
        return Point(affine_model, **result)

# ===========================================================================================================================
def compare_kummerline_points(affine_expected, affine_result, params):

    # Handling special point negation
    if affine_expected == params.curve.neutral:
        affine_expected_neg = params.curve.neutral
    else:
        affine_expected_neg = params.curve.affine_negate(affine_expected)

    # Handling special coordinate
    if params.curve.coordinate_model.name == "yzsquared":
        if affine_expected == params.curve.neutral:
            affine_expected_tmp = params.curve.neutral
            affine_expected_aux = params.curve.neutral
        else:
            affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
            tmp = {"x":affine_expected.x, "y":-affine_expected.y}
            affine_expected_tmp = Point(affine_model, **tmp)
            assert(params.curve.is_on_curve(affine_expected_tmp))
            affine_expected_aux = params.curve.affine_negate(affine_expected_tmp)
            assert(params.curve.is_on_curve(affine_expected_aux))
    else:
        affine_expected_tmp = affine_expected
        affine_expected_aux = affine_expected_neg

    return affine_result == affine_expected or affine_result == affine_expected_neg or affine_result == affine_expected_tmp or affine_result == affine_expected_aux

# ===========================================================================================================================
# Just for testing the differential addition on the Kummer line
def test_dadd_formula(formula, params, num):
    # Kummer-line has elements in the projective space P^2
    assert( len(params.curve.coordinate_model.variables) == 2)
    affine_gen = kummer_to_affine(params.generator, params)

    # Special case for the YZsquared-coordinates
    if params.curve.coordinate_model.name == "yzsquared":
        if params.curve.affine_multiply(affine_gen, params.order - 1) != params.curve.affine_negate(affine_gen):
            affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
            affine_tmp = {"x":affine_gen.x, "y":-affine_gen.y}
            affine_gen = Point(affine_model, **affine_tmp)
            assert(params.curve.is_on_curve(affine_gen))
        
    ok = True
    ok_counterexample = []
    unified = True
    unified_counterexample = []
    complete = True #not isinstance(params.curve.neutral, InfinityPoint)
    complete_counterexample = []

    # Let's precompute the cyclic subgroup
    points = []
    assert(params.curve.affine_multiply(affine_gen, params.order - 1) == params.curve.affine_negate(affine_gen))
    points.append(affine_gen)
    assert(params.curve.is_on_curve(points[0]))
    points.append(params.curve.affine_double(points[0]))                    # [2]g
    assert(params.curve.is_on_curve(points[1]))
    for i in range(2, min(num, params.order) - 1, 1):
        points.append(params.curve.affine_add(points[0], points[i - 1]))    # [i + 1]g
        assert(params.curve.is_on_curve(points[i]))

    digits = len(str(int(min(num, params.order))))
    for i in range(1, min(num, params.order)):
        #affine_P = params.curve.affine_multiply(affine_gen, i)
        affine_P = points[i - 1]
        P = affine_P.to_model(params.curve.coordinate_model, params.curve)

        # Do P + Q for all P,Q on the curve, s.t. P != Q and P + Q != infty then compare with affine add
        for j in range(1, min(num, params.order)):
            # TODO: Test also the case when P + Q = infty
            if i == j or i + j == params.order:
                continue
            #affine_Q = params.curve.affine_multiply(affine_gen, j)
            affine_Q = points[j - 1]
            Q = affine_Q.to_model(params.curve.coordinate_model, params.curve)
            affine_PQ= params.curve.affine_add(affine_P, params.curve.affine_negate(affine_Q))
            PQ = affine_PQ.to_model(params.curve.coordinate_model, params.curve)
            
            result = formula(params.curve.prime, PQ, P, Q, **params.curve.parameters)[0]
            affine_result = kummer_to_affine(result, params)
            affine_expected = params.curve.affine_add(affine_P, affine_Q)
            if not compare_kummerline_points(affine_expected, affine_result, params):
                ok = False
                ok_counterexample.append(\
                        point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.order)[2:])) \
                        + '; ' + \
                        point2str(f'[{j:{digits}d}]P', affine_Q, len(hex(params.order)[2:])) \
                        + '; ' + \
                        point2str(f'[{i:{digits}d}]P - [{j:{digits}d}]P', affine_PQ, len(hex(params.order)[2:]))
                        )

        # Do P + P and compare with affine doubling of P. If this works, the formulas are unified (but this will never occurs, I guess ...)
        try:
            result = formula(params.curve.prime, params.curve.neutral, P, P, **params.curve.parameters)[0]
            affine_result = kummer_to_affine(result, params)
            # Handling correct point doubling computation of the point at infinity
            if params.curve.affine_negate(affine_P) == affine_P:
                if params.curve.model.shortname == 'twisted' or params.curve.model.shortname == 'edwards':
                    affine_expected = params.curve.affine_double(affine_P)
                else:
                    affine_expected = params.curve.neutral
            else:
                affine_expected = params.curve.affine_double(affine_P)

            if not compare_kummerline_points(affine_expected, affine_result, params):
                unified = False
                unified_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )
        except:
            unified = False
            unified_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )

        # Do P + infty, and infty + P, if infty is representable in coordinate system. If this works, the formulas are complete
        result1 = formula(params.curve.prime, P, P, params.curve.neutral, **params.curve.parameters)[0]
        result2 = formula(params.curve.prime, P, params.curve.neutral, P, **params.curve.parameters)[0]
        affine_result1 = kummer_to_affine(result1, params)
        affine_result2 = kummer_to_affine(result2, params)
        if not compare_kummerline_points(affine_P, affine_result1, params) or not compare_kummerline_points(affine_P, affine_result2, params):
            complete = False
            complete_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )

    return {\
            'P + Q' : { 'Are there failed computations?' : not ok,
                'Failures' : ok_counterexample, 
                'Number of failures' : len(ok_counterexample)},\
            'P + P' : { 'Are there failed computations?' : not unified, 
                'Failures' : unified_counterexample if len(unified_counterexample) < (min(num, params.order) - 1) else "All the points (not including the point at infinity)", 
                'Number of failures' : len(unified_counterexample) },\
            'P + infty and infty + P' : { 'Are there failed computations?' : not complete, 
                'Failures' : complete_counterexample if len(complete_counterexample) < (min(num, params.order) - 1) else "All the points (not including the point at infinity)", 
                'Number of failures' : len(complete_counterexample) }\
            }

# ===========================================================================================================================
# Just for testing the simultaneous computation of dadd and dbl on the Kummer line
well_defined = lambda x: {True:'wrong', False:'okay'}[x]
def test_ladd_formula(formula, params, num):
    # Kummer-line has elements in the projective space P^2
    assert( len(params.curve.coordinate_model.variables) == 2 )
    affine_gen = kummer_to_affine(params.generator, params)
    
    # Special case for the YZsquared-coordinates
    if params.curve.coordinate_model.name == "yzsquared":
        if params.curve.affine_multiply(affine_gen, params.order - 1) != params.curve.affine_negate(affine_gen):
            affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
            affine_tmp = {"x":affine_gen.x, "y":-affine_gen.y}
            affine_gen = Point(affine_model, **affine_tmp)
            assert(params.curve.is_on_curve(affine_gen))

    ok = True
    ok_counterexample = []

    # Let's precompute the cyclic subgroup
    points = []
    assert(params.curve.affine_multiply(affine_gen, params.order - 1) == params.curve.affine_negate(affine_gen))
    points.append(affine_gen)
    assert(params.curve.is_on_curve(points[0]))
    points.append(params.curve.affine_double(points[0]))                    # [2]g
    assert(params.curve.is_on_curve(points[1]))
    for i in range(2, min(num, params.order) - 1, 1):
        points.append(params.curve.affine_add(points[0], points[i - 1]))    # [i + 1]g
        assert(params.curve.is_on_curve(points[i]))

    digits = len(str(int(min(num, params.order))))
    for i in range(1, min(num, params.order)):
        #affine_P = params.curve.affine_multiply(affine_gen, i)
        affine_P = points[i - 1]
        P = affine_P.to_model(params.curve.coordinate_model, params.curve)
        
        # Do P + Q for all P,Q on the curve, s.t. P != Q and P + Q != infty then compare with affine add
        for j in range(1, min(num, params.order)):
            # TODO: Test also the case when P + Q = infty
            if i == j or i + j == params.order:
                continue
            #affine_Q = params.curve.affine_multiply(affine_gen, j)
            affine_Q = points[j - 1]
            Q = affine_Q.to_model(params.curve.coordinate_model, params.curve)
            affine_PQ= params.curve.affine_add(affine_P, params.curve.affine_negate(affine_Q))
            PQ = affine_PQ.to_model(params.curve.coordinate_model, params.curve)
            
            result1, result2 = formula(params.curve.prime, PQ, P, Q, **params.curve.parameters)[:2]
            affine_result1 = kummer_to_affine(result1, params)
            affine_result2 = kummer_to_affine(result2, params)

            # Catching affine point doubling with input 2-order point
            if params.curve.affine_negate(affine_P) == affine_P:
                if params.curve.model.shortname == 'twisted' or params.curve.model.shortname == 'edwards':
                    affine_expected1 = params.curve.affine_double(affine_P)
                else:
                    affine_expected1 = params.curve.neutral
            else:
                affine_expected1 = params.curve.affine_double(affine_P)
            affine_expected2 = params.curve.affine_add(affine_P, affine_Q)

            if not compare_kummerline_points(affine_expected1, affine_result1, params) or not compare_kummerline_points(affine_expected2, affine_result2, params):
                ok = False
                ok_counterexample.append(\
                        point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.order)[2:])) \
                        + '; ' + \
                        point2str(f'[{j:{digits}d}]P', affine_Q, len(hex(params.order)[2:])) \
                        + '; ' + \
                        point2str(f'[{i:{digits}d}]P - [{j:{digits}d}]P', affine_PQ, len(hex(params.order)[2:])) \
                        + '; {' + \
                        'dbl: ' + well_defined(affine_result1 != affine_expected1 and affine_result1 != params.curve.affine_negate(affine_expected1)) \
                        + ', ' + \
                        'add: ' + well_defined((affine_result2 != affine_expected2 and affine_result2 != params.curve.affine_negate(affine_expected2))) \
                        + '}'
                        )

    return {\
            'P + Q' : { 'Are there failed computations?' : not ok,
                'Failures' : ok_counterexample, 
                'Number of failures' : len(ok_counterexample)},\
            }

# ===========================================================================================================================
# Testing the differential doubling on the Kummer line
def test_ddbl_formula(formula, params, num):
    # Kummer-line has elements in the projective space P^2
    assert( len(params.curve.coordinate_model.variables) == 2 )
    affine_gen = kummer_to_affine(params.generator, params)
    
    # Special case for the YZsquared-coordinates
    if params.curve.coordinate_model.name == "yzsquared":
        if params.curve.affine_multiply(affine_gen, params.order - 1) != params.curve.affine_negate(affine_gen):
            affine_model = AffineCoordinateModel(params.curve.coordinate_model.curve_model)
            affine_tmp = {"x":affine_gen.x, "y":-affine_gen.y}
            affine_gen = Point(affine_model, **affine_tmp)
            assert(params.curve.is_on_curve(affine_gen))

    ok = True
    ok_counterexample = []
    # Let's precompute the cyclic subgroup
    points = []
    assert(params.curve.affine_multiply(affine_gen, params.order - 1) == params.curve.affine_negate(affine_gen))
    points.append(affine_gen)
    assert(params.curve.is_on_curve(points[0]))
    points.append(params.curve.affine_double(points[0]))                    # [2]g
    assert(params.curve.is_on_curve(points[1]))
    for i in range(2, min(num, params.order) - 1, 1):
        points.append(params.curve.affine_add(points[0], points[i - 1]))    # [i + 1]g
        assert(params.curve.is_on_curve(points[i]))

    digits = len(str(int(min(num, params.order))))
    for i in range(1, min(num, params.order)):
        #affine_P = params.curve.affine_multiply(affine_gen, i)
        affine_P = points[i - 1]
        P = affine_P.to_model(params.curve.coordinate_model, params.curve)
        result = formula(params.curve.prime, P, **params.curve.parameters)[0]
        affine_result = kummer_to_affine(result, params)

        # Handling correct point doubling computation of the point at infinity
        if params.curve.affine_negate(affine_P) == affine_P:
            if params.curve.model.shortname == 'twisted' or params.curve.model.shortname == 'edwards':
                affine_expected = params.curve.affine_double(affine_P)
            else:
                affine_expected = params.curve.neutral
        else:
            affine_expected = params.curve.affine_double(affine_P)

        if not compare_kummerline_points(affine_expected, affine_result, params):
            ok = False
            ok_counterexample.append( point2str(f'[{i:{digits}d}]P', affine_P, len(hex(params.curve.prime)[2:])) )

    return { '[2]P' : { 'Are there failed computations?' : not ok, 
        'Failures' : ok_counterexample if len(ok_counterexample) < (min(num, params.order) - 1) else "All the points (not including the point at infinity)", 
        'Number of failures' : len(ok_counterexample) } }

# ===========================================================================================================================
# ===========================================================================================================================
import sys

@click.command()
@click.option('--std', default=None, type=str, help="File with extension .json containing the curves to be used (same format as in std-curve).")
@click.option('--efd', default=None, type=str, help="Formula to be used (same name as in the efd collection).")
@click.option('--num', default=None, type=int, help="Number of point addition to be performed.")
def main(std, efd, num):
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
        "add": test_add_formula,
        "dbl": test_dbl_formula
    }
    tests_kummer = {
        "dadd": test_dadd_formula,
        "ladd": test_ladd_formula,
        "dbl": test_ddbl_formula,
    }

    curves, length = json2curves(std) # Reading the file with extension .json that contains the list of curves to be used

    # ===========================================================================================================================
    for i in range(0, length, 1):

        FORMULA_NOT_FOUND = True    # This is order to ensure the dseired formula belongs to the EFD
        message = "PROCESSING: " + curves['curves'][i]['desc']
        print("-" * len(message)  + "\n" + message)
        with open('%s-%s.log' % (curves['curves'][i]['name'], efd), 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            curve_name = curves['curves'][i]['name']
            model = model_list[curves['curves'][i]['form']]
            assert( model in model_list.values() )
            print("=" * len(curves['curves'][i]['desc'])  + "\n" + curves['curves'][i]['desc'])
            for coord_name, coords in model.coordinates.items():

                try:
                    curve = get_curve(curves, i, coord_name)    # now infty = False (projective representation of the point at infinty)
                except ValueError as e:
                    click.echo(f"{coord_name} {curve_name} {e}")
                    continue

                #print(efd, coords.formulas.keys())
                for formula_name, formula in coords.formulas.items():
                    if efd == formula_name:
                        FORMULA_NOT_FOUND = False
                        if (formula.shortname in tests) or (formula.shortname in tests_kummer):
                            try:
                                if len(coords.variables) > 2:
                                    # Rational points: point addition and doubling
                                    res = tests[formula.shortname](formula, curve, num)
                                else:
                                    # Kummer line: differential point additions and doubling
                                    assert( len(coords.variables) == 2)
                                    res = tests_kummer[formula.shortname](formula, curve, num)
                            except NameError as e:
                                click.echo(f"{coord_name} {curve_name} {formula_name} Bad formula! {e}")
                                continue
                            print( "-" * (len(coord_name + curve_name + formula_name) + 2))
                            click.echo(f"{coord_name} {curve_name} {formula_name}\n" + yaml.dump(res, default_flow_style=False, width=float("inf")))

            sys.stdout = original_stdout # Reset the standard output to its original value

        if FORMULA_NOT_FOUND:
            click.echo(f"WARNING: {efd} is not a valid formula for {model_names[curves['curves'][i]['form']]}")

if __name__ == "__main__":
    main()

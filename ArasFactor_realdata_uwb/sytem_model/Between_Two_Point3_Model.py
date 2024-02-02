# Setup
import numpy as np
import symforce
import sympy as sp
import os 

symforce.set_symbolic_api("symengine")
symforce.set_log_level("warning")

symforce.set_epsilon_to_symbol()

from symforce import codegen
from symforce.codegen import codegen_util
from symforce import ops
import symforce.symbolic as sf
from symforce.values import Values
from symforce.notebook_util import display, display_code, display_code_file
from itertools import permutations

out_put_save_directory = os.getcwd()

def calculate_norm(vector):
    return sf.sqrt(sum(component**2 for component in vector))

a = sf.V3.symbolic("a")
b = sf.V3.symbolic("b")
ofset, distance_measure = sf.symbols("ofset distance_measure")

cost = sf.M11((a - b).norm() + ofset - distance_measure)

def error_model_between_tow_point(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), ofset: sf.Symbol('ofset'), distance_measure: sf.Symbol('distance_measure'), epsilon: sf.Scalar = 0
                ) -> sf.Vector1:
    return sf.V1(cost) 
resedual_func_codegen = codegen.Codegen.function(func=error_model_between_tow_point, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

def error_model_between_tow_point_wrt_point1(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), ofset: sf.Symbol('ofset'), distance_measure: sf.Symbol('distance_measure'), epsilon: sf.Scalar = 0
                ) -> sf.Vector3:
    vector_cost_wrt_point1 = cost.jacobian(a)
    return sf.M13(vector_cost_wrt_point1) 
resedual_func_codegen = codegen.Codegen.function(func=error_model_between_tow_point_wrt_point1, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

def error_model_between_tow_point_wrt_point2(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), ofset: sf.Symbol('ofset'), distance_measure: sf.Symbol('distance_measure'), epsilon: sf.Scalar = 0
                ) -> sf.Vector3:
    vector_cost_wrt_point2 = cost.jacobian(b)
    return sf.M13(vector_cost_wrt_point2) 
resedual_func_codegen = codegen.Codegen.function(func=error_model_between_tow_point_wrt_point2, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

def error_model_between_tow_point_wrt_ofset(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), ofset: sf.Symbol('ofset'), distance_measure: sf.Symbol('distance_measure'), epsilon: sf.Scalar = 0
                ) -> sf.Vector3:
    vector_cost_wrt_ofset = cost.diff(ofset)
    return sf.M11(vector_cost_wrt_ofset) 
resedual_func_codegen = codegen.Codegen.function(func=error_model_between_tow_point_wrt_ofset, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

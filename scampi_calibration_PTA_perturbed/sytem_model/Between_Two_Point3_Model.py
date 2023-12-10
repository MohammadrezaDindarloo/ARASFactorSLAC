# Setup
import numpy as np
import symforce
import sympy as sp

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

out_put_save_directory = "/home/mohammad/scampi_calibration/scampi_factor_graph"

def calculate_norm(vector):
    return sf.sqrt(sum(component**2 for component in vector))

a = sf.V3.symbolic("a")
b = sf.V3.symbolic("b")
vector_cost_wrt_point1 = sf.V3.symbolic("")
vector_cost_wrt_point2 = sf.V3.symbolic("")

cost = (a - b).norm()

vector_cost_wrt_point1[0] = cost.diff(a[0])
vector_cost_wrt_point1[1] = cost.diff(a[1])
vector_cost_wrt_point1[2] = cost.diff(a[2])

vector_cost_wrt_point2[0] = cost.diff(b[0])
vector_cost_wrt_point2[1] = cost.diff(b[1])
vector_cost_wrt_point2[2] = cost.diff(b[2])

print("**********************")
print((cost))
print("**********************")
print((vector_cost_wrt_point1))
print("**********************")
print((vector_cost_wrt_point2))

def error_model(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), epsilon: sf.Scalar = 0
                ) -> sf.Vector1:
    return sf.V1(cost) 
resedual_func_codegen = codegen.Codegen.function(func=error_model, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


def error_model_wrt_point1(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), epsilon: sf.Scalar = 0
                ) -> sf.Vector3:
    return sf.V3(vector_cost_wrt_point1) 
resedual_func_codegen = codegen.Codegen.function(func=error_model_wrt_point1, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

def error_model_wrt_point2(a: sf.V3.symbolic('a'), b: sf.V3.symbolic('b'), epsilon: sf.Scalar = 0
                ) -> sf.Vector3:
    return sf.V3(vector_cost_wrt_point2) 
resedual_func_codegen = codegen.Codegen.function(func=error_model_wrt_point2, config=codegen.CppConfig(),)
resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


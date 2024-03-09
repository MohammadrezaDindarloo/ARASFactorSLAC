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

# data structure
class GeometricVariables:
    def __init__(self):
        self.p_in_w = []   # Pulley locations in world frame
        self.b_in_w = []   # Body attachments in world frame
        self.b_rot = []    # Rotated Body attachment points to W
        self.sx = []       # Body to pulley directions projected on xy plane
        self.s=[]           #unit vector along cable direction (assuming no sag)
        self.r_to_cog = sf.Vector3(0, 0, 0)


class CatenaryVariables:
    def __init__(self):
        self.length = []  # The true distance between the anchor locations and body attachment points
        self.c1 = []      # Catenary parameter 1
        self.c2 = []      # Catenary parameter 2
        self.y0_cat = []  # Z component of the pulley location in world coordinate
        self.yl_cat = []  # Z component of the body attachment point in world coordinate
        self.lc_cat = []  # The true label length

class RobotState:
    def __init__(self):
        self.p_platform = sf.Vector3.symbolic("P")  # Pulley locations in world frame
        self.rot_platform = sf.Rot3.symbolic("R")
        self.cable_forces = []    # Horizontal and vertical cable forces at body attachment point
        self.cable_forces_compact = sf.Matrix([0] * 8)
        self.wrench = sf.Matrix([0] * 6)
        self.static_constrain = sf.Matrix21.zero()  # The static equilibrium constrain

class RobotParameters:
    def __init__(self):
        self.pulleys = []   # Pulley locations in world frame
        self.ef_points = [] # Cable attachment points in ef frame
        self.r_to_cog = sf.Matrix([0, 0, 0])
        self.g_c = 0.0
        self.f_g = 0.0

class CatDataOut:
    def __init__(self):
        self.b_in_w = []    # Body attachments in world frame
        self.cable_forces = []    # Horizontal and vertical cable forces at body attachment point
        self.c1 = []      # Catenary parameter 1
        self.c2 = []      # Catenary parameter 2
        self.lc_cat = []  # The true label length
        self.rot_platform = sf.Rot3.symbolic("R_cat_out")

class IKDataOut:
    def __init__(self):
        self.b_in_w = []    # Body attachments in world frame
        self.cable_forces = []    # Horizontal and vertical cable forces at body attachment point
        self.c1 = []      # Catenary parameter 1
        self.c2 = []      # Catenary parameter 2
        self.lc_cat = []  # The true label length
        self.rot_platform = sf.Rot3.symbolic("R_IK_out")

class FKDataOut:
    def __init__(self):
        self.b_in_w = []    # Body attachments in world frame
        self.cable_forces = []    # Horizontal and vertical cable forces at body attachment point
        self.c1 = []      # Catenary parameter 1
        self.c2 = []      # Catenary parameter 2
        self.lc_cat = []  # The true label length
        self.p_platform = sf.Matrix([0, 0, 0])  # The optimized end-effector position
        self.rot_platform = sf.Rot3.symbolic("R_FK_out")  # The optimized end-effector orientation



numbers = "0123"
permutations_list = list(permutations(numbers))

name_list = [f"l{''.join(order)}" for order in permutations_list]

from itertools import permutations

numbers = "0123"
permutations_list = list(permutations(numbers))

list_names_without_l = [f"{''.join(order)}" for order in permutations_list]



# utilities functions
def calculate_norm(vector):
    return sf.sqrt(sum(component**2 for component in vector))

    # Compute the geometric variables based on the robot's state and parameters
def getGeometricVariables(robotstate_, robotparams_, geometric_vars: GeometricVariables):
    geometric_vars.b_rot.clear()
    geometric_vars.b_in_w.clear()
    geometric_vars.sx.clear()
    geometric_vars.s.clear()

    for i in range(len(robotparams_.pulleys)):
        # Rotate the body attachment points into the world frame
        rot_b = robotstate_.rot_platform.to_rotation_matrix() * sf.Matrix(robotparams_.ef_points[i])
        geometric_vars.b_rot.append(rot_b)
        # translate the rotated body point to the end-effector location
        b_in_w = rot_b + sf.Matrix(robotstate_.p_platform)
        geometric_vars.b_in_w.append(b_in_w)
        # compute the x-y projection of the direction from end-effector attachment point to the pulley
        sx = sf.Matrix(b_in_w) - sf.Matrix(robotparams_.pulleys[i])
        s = sf.Matrix(b_in_w) - sf.Matrix(
            robotparams_.pulleys[i])  # unit vector along cable direction (assuming no sag)
        sx[2] = 0.0
        geometric_vars.sx.append(sx / calculate_norm(sx))
        geometric_vars.s.append(s / calculate_norm(s))
        geometric_vars.p_in_w.append(sf.Matrix(robotparams_.pulleys[i]))

    geometric_vars.r_to_cog = robotstate_.rot_platform.to_rotation_matrix() * robotparams_.r_to_cog

    # Compute the catenary variables from the robot's state and geometric parameters
def getCatenaryVariables(robot_state, robot_params, geometric_vars, cat_vars):
    gc = robot_params.g_c  # for the sake of clarity

    cat_vars.length.clear()
    cat_vars.c1.clear()
    cat_vars.c2.clear()
    cat_vars.lc_cat.clear()
    cat_vars.y0_cat.clear()
    cat_vars.yl_cat.clear()

    for i in range(len(robot_params.pulleys)):
        # Compute projection length of the vector from the body to the pulley on the x-y plane
        p = geometric_vars.p_in_w[i][:2]
        b = geometric_vars.b_in_w[i][:2]

        L = calculate_norm(p-b)
        cat_vars.length.append(L)

        # Compute the C1 parameter
        fh = robot_state.cable_forces[i][0]
        fv = robot_state.cable_forces[i][1]
        C1 = fh / gc * sf.asinh(-fv / fh) - L

        cat_vars.c1.append(C1)

        # Extract the Z components of the pulleys and body attachment points
        y0_cat = geometric_vars.p_in_w[i][2]
        yl_cat = geometric_vars.b_in_w[i][2]
        cat_vars.y0_cat.append(y0_cat)
        cat_vars.yl_cat.append(yl_cat)

        # Compute the C2 parameter
        C2 = sf.cosh(C1 * gc / fh) - gc / fh * y0_cat

        cat_vars.c2.append(C2)

        # Compute the true cable length
        lc_cat = fh / gc * (sf.sinh((gc / fh) * (L + C1)) - sf.sinh((gc / fh) * C1))
        cat_vars.lc_cat.append(lc_cat)
    
# A function to make the A matrix for the getCableForces function
def constructStructuralMat(sx, b_rot):
    const_vec = sf.Matrix([0, 0, 1])
    A = sf.Matrix68.zero()

    for i in range(len(sx)):
        A[0:3, 2*i:2*i+1] = -sx[i]
        A[3:6, 2*i:2*i+1] = -b_rot[i].cross(sx[i])
        A[0:3, 2*i+1:2*i+2] = const_vec
        A[3:6, 2*i+1:2*i+2] = b_rot[i].cross(const_vec)

    return A

def constructStructuralMat_Non_Flatten(s, b_rot):
    const_vec = sf.Matrix([0, 0, 1])
    A = sf.Matrix64.zero()

    for i in range(4):
        A[0:3, i] = s[i]
        A[3:6, i] = b_rot[i].cross(s[i])
    return A

# Compute the catenary variables from the robot's state and geometric parameters
def getCableForces(fh, fv, robotstate_:RobotState, robotparams_, geometric_vars):
    structure_matrix = sf.Matrix68.zero()
    const_vec = sf.Matrix([0, 0, -1])
    # The wrench applied to the end-effector
    platform_wrench = sf.Matrix([0, 0, robotparams_.f_g, 0, 0, 0])

    # The cog modification recently added to the matlab implementation for the real data experiment
    
    platform_wrench[3:6, 0] = - geometric_vars.r_to_cog.cross(const_vec) * (robotparams_.f_g)
    # End of cog modification

    f_v = sf.Matrix([fh, fv])

    structure_matrix = constructStructuralMat(geometric_vars.sx, geometric_vars.b_rot)

    structure_matrix_Non_Flatten=constructStructuralMat_Non_Flatten(geometric_vars.s, geometric_vars.b_rot)

    J_aT = structure_matrix_Non_Flatten[0:4, :]
    J_uT = structure_matrix_Non_Flatten[4:6, :]
    W_a = structure_matrix_Non_Flatten[0:4, 0]
    W_u = structure_matrix_Non_Flatten[4:6, 0]

    # a subset of the structure matrix that defines the force in other cables as a function of forces of the first cable
    sub_structure_matrix = structure_matrix[0:6, 2:8]
    other_forces = sub_structure_matrix.inv() * (platform_wrench - structure_matrix[0:6, 0:2] * f_v[0:2, 0])
    cable_forces = sf.Matrix81.zero()
    cable_forces[0:2, 0] = f_v
    cable_forces[2:8, 0] = other_forces

    static_constrain=J_uT*J_aT.inv()*W_a-W_u
    robotstate_.static_constrain=static_constrain


    F = structure_matrix[0:3, 0:8] * cable_forces - sf.Matrix([0, 0, robotparams_.f_g])
    Torque = structure_matrix[3:6, 0:8] * cable_forces
    robotstate_.wrench[0:3, 0] = F
    robotstate_.wrench[3:6, 0] = Torque

    robotstate_.cable_forces_compact = cable_forces
    for i in range(len(geometric_vars.sx)):
        robotstate_.cable_forces.append(cable_forces[2*i:2*i+2, 0])


def get_cable_order(largest_cable):
    if set(largest_cable) != set('0123'):
        raise ValueError(f'Invalid input for cable length. '
                 f'The cable number is incorrect, largest cable should be one of the permutations of 0, 1, 2, 3. Input: {largest_cable}')
    order = [int(c) for c in largest_cable]
    name = f'l{largest_cable}'
    display_order(order)
    return order, name

def display_order(order):
    print("Cable order:", order)


# A function for the IK solver that puts the cable with the larges length as the first cable for numerical stability
def changeOrderForSolver(state, params, largest_cable):
    N = len(params.pulleys)
    order = sf.Matrix([i for i in range(N)])

    if largest_cable == '0123':
        order[0] = 0
        order[1] = 1
        order[2] = 2
        order[3] = 3
        name = 'l0123'

    elif largest_cable == '0132':
        order[0] = 0
        order[1] = 1
        order[2] = 3
        order[3] = 2
        name = 'l0132'

    elif largest_cable == '0213':
        order[0] = 0
        order[1] = 2
        order[2] = 1
        order[3] = 3
        name = 'l0213'
        
    elif largest_cable == '0231':
        order[0] = 0
        order[1] = 2
        order[2] = 3
        order[3] = 1
        name = 'l0231'

    elif largest_cable == '0312':
        order[0] = 0
        order[1] = 3
        order[2] = 1
        order[3] = 2
        name = 'l0312'

    elif largest_cable == '0321':
        order[0] = 0
        order[1] = 3
        order[2] = 2
        order[3] = 1
        name = 'l0321'

    elif largest_cable == '1023':
        order[0] = 1
        order[1] = 0
        order[2] = 2
        order[3] = 3
        name = 'l1023'

    elif largest_cable == '1032':
        order[0] = 1
        order[1] = 0
        order[2] = 3
        order[3] = 2
        name = 'l1032'

    elif largest_cable == '1203':
        order[0] = 1
        order[1] = 2
        order[2] = 0
        order[3] = 3
        name = 'l1203'

    elif largest_cable == '1230':
        order[0] = 1
        order[1] = 2
        order[2] = 3
        order[3] = 0
        name = 'l1230'

    elif largest_cable == '1302':
        order[0] = 1
        order[1] = 3
        order[2] = 0
        order[3] = 2
        name = 'l1302'

    elif largest_cable == '1320':
        order[0] = 1
        order[1] = 3
        order[2] = 2
        order[3] = 0
        name = 'l1320'

    elif largest_cable == '2013':
        order[0] = 2
        order[1] = 0
        order[2] = 1
        order[3] = 3
        name = 'l2013'

    elif largest_cable == '2031':
        order[0] = 2
        order[1] = 0
        order[2] = 3
        order[3] = 1
        name = 'l2031'

    elif largest_cable == '2103':
        order[0] = 2
        order[1] = 1
        order[2] = 0
        order[3] = 3
        name = 'l2103'

    elif largest_cable == '2130':
        order[0] = 2
        order[1] = 1
        order[2] = 3
        order[3] = 0
        name = 'l2130'

    elif largest_cable == '2301':
        order[0] = 2
        order[1] = 3
        order[2] = 0
        order[3] = 1
        name = 'l2301'

    elif largest_cable == '2310':
        order[0] = 2
        order[1] = 3
        order[2] = 1
        order[3] = 0
        name = 'l2310'

    elif largest_cable == '3012':
        order[0] = 3
        order[1] = 0
        order[2] = 1
        order[3] = 2
        name = 'l3012'

    elif largest_cable == '3021':
        order[0] = 3
        order[1] = 0
        order[2] = 2
        order[3] = 1
        name = 'l3021'

    elif largest_cable == '3102':
        order[0] = 3
        order[1] = 1
        order[2] = 0
        order[3] = 2
        name = 'l3102'

    elif largest_cable == '3120':
        order[0] = 3
        order[1] = 1
        order[2] = 2
        order[3] = 0
        name = 'l3120'

    elif largest_cable == '3201':
        order[0] = 3
        order[1] = 2
        order[2] = 0
        order[3] = 1
        name = 'l3201'

    elif largest_cable == '3210':
        order[0] = 3
        order[1] = 2
        order[2] = 1
        order[3] = 0
        name = 'l3210'
        
    else:
        raise NameError('invalid input for cable lenght')
        print("The cable number is incorrect, largest cable should be one of the 1, 2, 3, 0")
    display(order)
    
    # Create a reordered params data structure
    params_reordered = RobotParameters
    params_reordered.r_to_cog = params.r_to_cog
    params_reordered.g_c = params.g_c
    params_reordered.f_g = params.f_g
    
    params_reordered.pulleys = [params.pulleys[int(i)] for i in order]
    params_reordered.ef_points = [params.ef_points[int(i)] for i in order]

    new_order = order
    return params_reordered, new_order, name

# Return the order of the cables back to the normal configuration
def reverseOrderForSolver(robotstate_, geometric_vars, cat_vars, order):
    N = len(geometric_vars.b_in_w)
    reorder_idx = sf.Matrix([0 for _ in range(N)])
    for i in range(N):
        reorder_idx[order[i]] = i

    fixed_cat_out = IKDataOut()

    for i in range(N):
        fixed_cat_out.b_in_w.append(geometric_vars.b_in_w[reorder_idx[i]])
        fixed_cat_out.c1.append(cat_vars.c1[reorder_idx[i]])
        fixed_cat_out.c2.append(cat_vars.c2[reorder_idx[i]])
        fixed_cat_out.cable_forces.append(robotstate_.cable_forces[reorder_idx[i]])
        fixed_cat_out.lc_cat.append(cat_vars.lc_cat[reorder_idx[i]])

    return fixed_cat_out

# Compute an initial value for the cable forces for the IK solver
def computeInitCableForces(p_platform, rot_platform, robotparams_):
    b_w = p_platform + rot_platform * robotparams_.ef_points[0]
    sx = b_w - robotparams_.pulleys[0]
    sx[2] = 0
    sx = sx / calculate_norm(sx)
    alpha = sf.acos((sx.transpose() * (b_w - robotparams_.pulleys[0]))[0] / calculate_norm((b_w - robotparams_.pulleys[0])))
    fv0 = -robotparams_.f_g / 4.0
    fh0 = abs(fv0) / sf.tan(alpha)
    return fh0, fv0



class CableRobotParams:
    def __init__(self, g_c, f_g):
        self.f_g_ = f_g
        self.g_c_ = g_c
        self.p1_ = sf.Matrix([0, 0, 0])
        self.p2_ = sf.Matrix([0, 0, 0])
        self.p3_ = sf.Matrix([0, 0, 0])
        self.p4_ = sf.Matrix([0, 0, 0])
        self.b1_ = sf.Matrix([0, 0, 0])
        self.b2_ = sf.Matrix([0, 0, 0])
        self.b3_ = sf.Matrix([0, 0, 0])
        self.b4_ = sf.Matrix([0, 0, 0])
        self.r_to_cog_ = sf.Matrix([0, 0, 0])

    def setPulleyPoses(self, p1, p2, p3, p4):
        self.p1_ = sf.Matrix(p1)
        self.p2_ = sf.Matrix(p2)
        self.p3_ = sf.Matrix(p3)
        self.p4_ = sf.Matrix(p4)

    def setEEAnchors(self, b1, b2, b3, b4):
        self.b1_ = sf.Matrix(b1)
        self.b2_ = sf.Matrix(b2)
        self.b3_ = sf.Matrix(b3)
        self.b4_ = sf.Matrix(b4)

    def setCog(self, r_to_cog):
        self.r_to_cog_ = sf.Matrix(r_to_cog)




def ikSolver(p_platform, rot_init, params, largest_cable, result):
    # reorder the cable forces and choose the cable with largest length as the the first cable (for numerical stability)
    reorder_idx = np.zeros(len(params.pulleys), dtype=int)
    params_reord = RobotParameters()
    state = RobotState()
    state.rot_platform = rot_init
    state.p_platform = p_platform
#     this name going to use for nameing the header files
    params_reord, reorder_idx, name = changeOrderForSolver(state, params, largest_cable)
#     Compute initil cable forces as starting points for the solver
    fh0, fv0 = computeInitCableForces(p_platform, rot_init, params_reord)
    
    # here we should create the cost function (state, params_reord, rot_init)
    state_ = state
    params_ = params_reord
    rot_init_ = rot_init

    fh1 = sf.Symbol('fh1')
    fv1 = sf.Symbol('fv1')
    DeltaRot = sf.Rot3.symbolic("DeltaRot")
    encoder = sf.Vector4.symbolic("encoder")

    gc = params_.g_c
    init_R = rot_init_

    geom_vars = GeometricVariables()
    cat_vars = CatenaryVariables()
    state = RobotState()
    params = RobotParameters()

    for i in range(len(params_.ef_points)):
        params.ef_points.append(params_.ef_points[i])
        params.pulleys.append(params_.pulleys[i])
    
    params.f_g = params_.f_g
    params.g_c = params_.g_c
    params.r_to_cog = params_.r_to_cog
    state.p_platform = state_.p_platform

    delta_rot = DeltaRot 
    state.rot_platform  = init_R * delta_rot
    
    getGeometricVariables(state, params, geom_vars)
    getCableForces(fh1, fv1, state, params, geom_vars)
    getCatenaryVariables(state, params, geom_vars, cat_vars)
    

    residuals = [0] * (12)
    for i in range(len(params.ef_points)):
        fh = state.cable_forces[i][0]
        residuals[i] = cat_vars.yl_cat[i] - fh/gc * ( sf.cosh(gc/fh * (cat_vars.length[i] + cat_vars.c1[i])) - cat_vars.c2[i] )
        residuals[i+4] = cat_vars.lc_cat[i] - calculate_norm(geom_vars.p_in_w[i]-geom_vars.b_in_w[i])
        residuals[i+8] = sf.exp(-fh)
    
    cost_z = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_z[i] = residuals[i] 
        
    cost_cat = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_cat[i] = residuals[i+4] 

    cost_force = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_force[i] = residuals[i+8] 
        
    if name == name_list[0]:
        print("--------------------generating header files name_list[0] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),                        
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_fh1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        
        def IK_residual_func_cost1_wrt_fv1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[0] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[0] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl0(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl0, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[0] IK----------")    

        
    if name == name_list[1]:
        print("--------------------generating header files name_list[1] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[1] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[1] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl1, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[1] IK----------") 
   

    if name == name_list[2]:
        print("--------------------generating header files name_list[2] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[2] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[2] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl2, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[2] IK----------") 


    if name == name_list[3]:
        print("--------------------generating header files name_list[3] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[3] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[3] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl3, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[3] IK----------") 
        
    if name == name_list[4]:
        print("--------------------generating header files name_list[4] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[4] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[4] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl4, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[4] IK----------") 


    if name == name_list[5]:
        print("--------------------generating header files name_list[5] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[5] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[5] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl5, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[5] IK----------") 
        
    if name == name_list[6]:
        print("--------------------generating header files name_list[6] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[6] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[6] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl6, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[6] IK----------") 
        
    if name == name_list[7]:
        print("--------------------generating header files name_list[7] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[7] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[7] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl7(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl7, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[7] IK----------") 
        
    if name == name_list[8]:
        print("--------------------generating header files name_list[8] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[8] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[8] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl8(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl8, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[8] IK----------") 
        
    if name == name_list[9]:
        print("--------------------generating header files name_list[9] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[9] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[9] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl9(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl9, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[9] IK----------") 
        
    if name == name_list[10]:
        print("--------------------generating header files name_list[10] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[10] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[10] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl10(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl10, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[10] IK----------") 
        
    if name == name_list[11]:
        print("--------------------generating header files name_list[11] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[11] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[11] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl11(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl11, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[11] IK----------") 
        
    if name == name_list[12]:
        print("--------------------generating header files name_list[12] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[12] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)



        def IK_residual_func_cost2_wrt_DeltaRot_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[12] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl12(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl12, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[12] IK----------") 
        
    if name == name_list[13]:
        print("--------------------generating header files name_list[13] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[13] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[13] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl13(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl13, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[13] IK----------") 

    if name == name_list[14]:
        print("--------------------generating header files name_list[14] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[14] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[14] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl14(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl14, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[14] IK----------") 
        
    if name == name_list[15]:
        print("--------------------generating header files name_list[15] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[15] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[15] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl15(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl15, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[15] IK----------") 
        
    if name == name_list[16]:
        print("--------------------generating header files name_list[16] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[16] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[16] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl16(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl16, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[16] IK----------") 
        
    if name == name_list[17]:
        print("--------------------generating header files name_list[17] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[17] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[17] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl17(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl17, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[17] IK----------") 
        
    if name == name_list[18]:
        print("--------------------generating header files name_list[18] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[18] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[18] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl18(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl18, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[18] IK----------") 
        
    if name == name_list[19]:
        print("--------------------generating header files name_list[19] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[19] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[19] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl19(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl19, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_forces header files generated name_list[19] IK----------") 
        
    if name == name_list[20]:
        print("--------------------generating header files name_list[20] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[20] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[20] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl20(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl20, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[20] IK----------") 

    if name == name_list[21]:
        print("--------------------generating header files name_list[21] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[21] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[21] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl21(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl21, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[21] IK----------") 
        
    if name == name_list[22]:
        print("--------------------generating header files name_list[22] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[22] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[22] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl22(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl22, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[22] IK----------") 
    
    if name == name_list[23]:
        print("--------------------generating header files name_list[23] IK---------------------")
        cost = cost_z
        def IK_residual_func_cost1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_DeltaRot_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_DeltaRot_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fh1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fh1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost1_wrt_fv1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_fv1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pa_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pa_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pb_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pb_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pc_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pc_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost1_wrt_pd_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost1_wrt_pd_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_z header files generated name_list[23] IK---------------")

        cost = cost_cat
        def IK_residual_func_cost2_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_DeltaRot_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_DeltaRot_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fh1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fh1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost2_wrt_fv1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_fv1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pa_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pa_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pb_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pb_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pc_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pc_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost2_wrt_pd_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost2_wrt_pd_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
        print("--------------------cost_catenary header files generated name_list[23] IK--------")

        cost = cost_force
        def IK_residual_func_cost3_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            return sf.V4(cost) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_DeltaRot_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_DeltaRot = cost.jacobian(DeltaRot)
            return sf.Matrix43(diff_DeltaRot)  
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_DeltaRot_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fh1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fh1 = cost.diff(fh1)
            return sf.V4(diff_fh1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fh1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


        def IK_residual_func_cost3_wrt_fv1_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                        position_vector: sf.Vector3.symbolic("position_vector"),
                        encoder: sf.V4.symbolic('encoder'),
                        p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                        Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                        ) -> sf.Vector4:
            diff_fv1 = cost.diff(fv1)
            return sf.V4(diff_fv1) 
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_fv1_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pa_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pa = cost.jacobian(p_a)
            return sf.Matrix43(diff_pa)

        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pa_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pb_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pb = cost.jacobian(p_b)
            return sf.Matrix43(diff_pb)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pb_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pc_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pc = cost.jacobian(p_c)
            return sf.Matrix43(diff_pc)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pc_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

        def IK_residual_func_cost3_wrt_pd_Nl23(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                            position_vector: sf.Vector3.symbolic("position_vector"),
                            encoder: sf.V4.symbolic('encoder'),
                            p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                            Rot_init: sf.Rot3.symbolic("Rot_init"), epsilon: sf.Scalar = 0
                            ) -> sf.Vector4:
            diff_pd = cost.jacobian(p_d)
            return sf.Matrix43(diff_pd)
        resedual_func_codegen = codegen.Codegen.function(func=IK_residual_func_cost3_wrt_pd_Nl23, config=codegen.CppConfig(),)
        resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory) 
        print("--------------------cost_forces header files generated name_list[23] IK----------") 
        print("********************** Generating IK cost function finished **********************") 



def fkSolver(lc_cat_measure, rtation_init, params_, fk_result):
    
    fh1 = sf.Symbol('fh1')
    fv1 = sf.Symbol('fv1')
    DeltaRot = sf.Rot3.symbolic("DeltaRot")
    TransformationMatrix = sf.Pose3.symbolic("TransformationMatrix")
    offset = sf.Vector4.symbolic("offset")

    gc = params_.g_c
    init_R = TransformationMatrix.rotation()
    geom_vars = GeometricVariables()
    cat_vars = CatenaryVariables()
    state = RobotState()
    params = RobotParameters()

    for i in range(len(params_.ef_points)):
        params.ef_points.append(params_.ef_points[i])
        params.pulleys.append(params_.pulleys[i])
    
    params.f_g = params_.f_g
    params.g_c = params_.g_c
    params.r_to_cog = params_.r_to_cog

    delta_rot = DeltaRot
    state.rot_platform  = init_R * delta_rot
    state.p_platform = sf.Vector3(TransformationMatrix.position()[0], TransformationMatrix.position()[1], TransformationMatrix.position()[2])
    
    getGeometricVariables(state, params, geom_vars)
    getCableForces(fh1, fv1, state, params, geom_vars)
    getCatenaryVariables(state, params, geom_vars, cat_vars)
    

    residuals = [0] * (22)
    for i in range(len(params.ef_points)):
        fh = state.cable_forces[i][0]
        residuals[i] = cat_vars.yl_cat[i] - fh/gc * ( sf.cosh(gc/fh * (cat_vars.length[i] + cat_vars.c1[i])) - cat_vars.c2[i] )
        residuals[i+4] = cat_vars.lc_cat[i] - calculate_norm(geom_vars.p_in_w[i]-geom_vars.b_in_w[i])
        residuals[i+8] = cat_vars.lc_cat[i] + offset[i] - lc_cat_measure[i]
        residuals[i+12] = calculate_norm(geom_vars.p_in_w[i]-geom_vars.b_in_w[i]) + offset[i] - lc_cat_measure[i]
        residuals[i+16] = cat_vars.y0_cat[i] - fh/gc * ( sf.cosh(gc/fh * (0.0 + cat_vars.c1[i])) - cat_vars.c2[i] )
    residuals[20]=state.static_constrain[0]
    residuals[21]=state.static_constrain[1]
    # calculate_norm(geom_vars.p_in_w[i]-geom_vars.b_in_w[i]) - lc_cat_measure[i]     # rigid cabel

    cost_z = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_z[i] = residuals[i] 
        
    cost_cat = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_cat[i] = residuals[i+4] 

    cost_encoder = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_encoder[i] = residuals[i+8]

    cost_straigth_encoder = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_straigth_encoder[i]= residuals[i+12]

    cost_z_0 = sf.Vector4.symbolic("P")
    for i in range(4):
        cost_z_0[i]= residuals[i+16]

    cost_static = sf.Vector2.symbolic("P")
    for i in range(2):
        cost_static[i]= residuals[i+20]

    print("\n--------------------generating header files FK---------------------")
    cost = cost_z
    def FK_residual_func_cost1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        return sf.V4(cost) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_DeltaRot(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_DeltaRot = cost.jacobian(DeltaRot)
        return sf.Matrix43(diff_DeltaRot)  
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_DeltaRot, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_fh1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fh1 = cost.diff(fh1)
        return sf.V4(diff_fh1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_fh1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_fv1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fv1 = cost.diff(fv1)
        return sf.V4(diff_fv1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_fv1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    
    def FK_residual_func_cost1_wrt_transformation(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_transformation = cost.jacobian(TransformationMatrix)
        return sf.Matrix46(diff_transformation)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_transformation, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_pa(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pa = cost.jacobian(p_a)
        return sf.Matrix43(diff_pa)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_pa, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_pb(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pb = cost.jacobian(p_b)
        return sf.Matrix43(diff_pb)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_pb, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_pc(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pc = cost.jacobian(p_c)
        return sf.Matrix43(diff_pc)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_pc, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost1_wrt_pd(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pd = cost.jacobian(p_d)
        return sf.Matrix43(diff_pd)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost1_wrt_pd, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    print("--------------------cost_z header files generated FK---------------")
    
    cost = cost_cat
    def FK_residual_func_cost2(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        return sf.V4(cost) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


    def FK_residual_func_cost2_wrt_DeltaRot(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_DeltaRot = cost.jacobian(DeltaRot)
        return sf.Matrix43(diff_DeltaRot)  
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_DeltaRot, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_fh1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fh1 = cost.diff(fh1)
        return sf.V4(diff_fh1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_fh1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_fv1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fv1 = cost.diff(fv1)
        return sf.V4(diff_fv1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_fv1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_transformation(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_transformation = cost.jacobian(TransformationMatrix)
        return sf.Matrix46(diff_transformation)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_transformation, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_pa(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pa = cost.jacobian(p_a)
        return sf.Matrix43(diff_pa)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_pa, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_pb(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pb = cost.jacobian(p_b)
        return sf.Matrix43(diff_pb)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_pb, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_pc(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pc = cost.jacobian(p_c)
        return sf.Matrix43(diff_pc)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_pc, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost2_wrt_pd(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pd = cost.jacobian(p_d)
        return sf.Matrix43(diff_pd)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost2_wrt_pd, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    print("--------------------cost_catenary header files generated FK--------")

    cost = cost_encoder
    def FK_residual_func_cost3(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        return sf.V4(cost) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


    def FK_residual_func_cost3_wrt_DeltaRot(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_DeltaRot = cost.jacobian(DeltaRot)
        return sf.Matrix43(diff_DeltaRot)  
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_DeltaRot, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost3_wrt_fh1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fh1 = cost.diff(fh1)
        return sf.V4(diff_fh1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_fh1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


    def FK_residual_func_cost3_wrt_fv1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fv1 = cost.diff(fv1)
        return sf.V4(diff_fv1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_fv1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    
    def FK_residual_func_cost3_wrt_transformation(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_transformation = cost.jacobian(TransformationMatrix)
        return sf.Matrix46(diff_transformation)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_transformation, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost3_wrt_pa(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pa = cost.jacobian(p_a)
        return sf.Matrix43(diff_pa)

    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_pa, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost3_wrt_pb(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pb = cost.jacobian(p_b)
        return sf.Matrix43(diff_pb)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_pb, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost3_wrt_pc(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pc = cost.jacobian(p_c)
        return sf.Matrix43(diff_pc)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_pc, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost3_wrt_pd(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pd = cost.jacobian(p_d)
        return sf.Matrix43(diff_pd)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_pd, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost3_wrt_offset(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_offset = cost.jacobian(offset)
        return sf.Matrix44(diff_offset)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost3_wrt_offset, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    print("--------------------cost_encoder header files generated FK----------")

    cost = cost_z_0
    def FK_residual_func_cost4(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        return sf.V4(cost) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_DeltaRot(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_DeltaRot = cost.jacobian(DeltaRot)
        return sf.Matrix43(diff_DeltaRot)  
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_DeltaRot, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_fh1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fh1 = cost.diff(fh1)
        return sf.V4(diff_fh1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_fh1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_fv1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fv1 = cost.diff(fv1)
        return sf.V4(diff_fv1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_fv1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    
    def FK_residual_func_cost4_wrt_transformation(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_transformation = cost.jacobian(TransformationMatrix)
        return sf.Matrix46(diff_transformation)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_transformation, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_pa(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pa = cost.jacobian(p_a)
        return sf.Matrix43(diff_pa)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_pa, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_pb(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pb = cost.jacobian(p_b)
        return sf.Matrix43(diff_pb)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_pb, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_pc(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pc = cost.jacobian(p_c)
        return sf.Matrix43(diff_pc)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_pc, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost4_wrt_pd(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pd = cost.jacobian(p_d)
        return sf.Matrix43(diff_pd)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost4_wrt_pd, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    print("-----------cost_z_0 header files generated FK----------")

    cost = cost_straigth_encoder
    def FK_residual_func_cost5(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        return sf.V4(cost) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


    def FK_residual_func_cost5_wrt_DeltaRot(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_DeltaRot = cost.jacobian(DeltaRot)
        return sf.Matrix43(diff_DeltaRot)  
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_DeltaRot, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost5_wrt_fh1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fh1 = cost.diff(fh1)
        return sf.V4(diff_fh1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_fh1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)


    def FK_residual_func_cost5_wrt_fv1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_fv1 = cost.diff(fv1)
        return sf.V4(diff_fv1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_fv1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    
    def FK_residual_func_cost5_wrt_transformation(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_transformation = cost.jacobian(TransformationMatrix)
        return sf.Matrix46(diff_transformation)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_transformation, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost5_wrt_pa(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pa = cost.jacobian(p_a)
        return sf.Matrix43(diff_pa)

    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_pa, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost5_wrt_pb(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pb = cost.jacobian(p_b)
        return sf.Matrix43(diff_pb)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_pb, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost5_wrt_pc(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pc = cost.jacobian(p_c)
        return sf.Matrix43(diff_pc)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_pc, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost5_wrt_pd(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_pd = cost.jacobian(p_d)
        return sf.Matrix43(diff_pd)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_pd, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost5_wrt_offset(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       offset: sf.V4.symbolic('offset'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector4:
        diff_offset = cost.jacobian(offset)
        return sf.Matrix44(diff_offset)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost5_wrt_offset, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    print("------------cost_encoder_straight header files generated FK--------")

    cost = cost_static
    def FK_residual_func_cost6(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        return sf.V2(cost) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_DeltaRot(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_DeltaRot = cost.jacobian(DeltaRot)
        return sf.Matrix23(diff_DeltaRot)  
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_DeltaRot, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_fh1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_fh1 = cost.diff(fh1)
        return sf.V2(diff_fh1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_fh1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_fv1(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_fv1 = cost.diff(fv1)
        return sf.V2(diff_fv1) 
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_fv1, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    
    def FK_residual_func_cost6_wrt_transformation(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"),
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_transformation = cost.jacobian(TransformationMatrix)
        return sf.Matrix26(diff_transformation)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_transformation, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_pa(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_pa = cost.jacobian(p_a)
        return sf.Matrix23(diff_pa)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_pa, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_pb(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_pb = cost.jacobian(p_b)
        return sf.Matrix23(diff_pb)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_pb, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_pc(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_pc = cost.jacobian(p_c)
        return sf.Matrix23(diff_pc)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_pc, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)

    def FK_residual_func_cost6_wrt_pd(fh1: sf.Symbol('fh1'), fv1: sf.Symbol('fv1'), 
                       DeltaRot: sf.Rot3.symbolic("DeltaRot"),
                       TransformationMatrix: sf.Pose3.symbolic("TransformationMatrix"),
                       encoder: sf.V4.symbolic('encoder'),
                       p_a: sf.V3.symbolic('p_a'), p_b: sf.V3.symbolic('p_b'), p_c: sf.V3.symbolic("p_c"),  p_d: sf.V3.symbolic("p_d"), 
                       epsilon: sf.Scalar = 0
                      ) -> sf.Vector2:
        diff_pd = cost.jacobian(p_d)
        return sf.Matrix23(diff_pd)
    resedual_func_codegen = codegen.Codegen.function(func=FK_residual_func_cost6_wrt_pd, config=codegen.CppConfig(),)
    resedual_func_codegen_data = resedual_func_codegen.generate_function(output_dir=out_put_save_directory)
    print("--------------------cost_static header files generated FK---------------")
    print("************************** FK & IK files **************************")


def inverseKinematicsSolver(cale_robo_param, p_platform, rot_init, largest_cable):
    params_ = RobotParameters()

    params_.pulleys.append(cale_robo_param.p1_)
    params_.pulleys.append(cale_robo_param.p2_)
    params_.pulleys.append(cale_robo_param.p3_)
    params_.pulleys.append(cale_robo_param.p4_)

    params_.ef_points.append(cale_robo_param.b1_)
    params_.ef_points.append(cale_robo_param.b2_)
    params_.ef_points.append(cale_robo_param.b3_)
    params_.ef_points.append(cale_robo_param.b4_)

    params_.r_to_cog = cale_robo_param.r_to_cog_
    params_.f_g = cale_robo_param.f_g_
    params_.g_c = cale_robo_param.g_c_

    ik_result = IKDataOut()
    ikSolver(p_platform, rot_init, params_, largest_cable, ik_result)    



def forwardKinematicsSolver(cale_robo_param, lc_cat, rtation_init):
    params_ = RobotParameters()

    params_.pulleys.append(cale_robo_param.p1_)
    params_.pulleys.append(cale_robo_param.p2_)
    params_.pulleys.append(cale_robo_param.p3_)
    params_.pulleys.append(cale_robo_param.p4_)

    params_.ef_points.append(cale_robo_param.b1_)
    params_.ef_points.append(cale_robo_param.b2_)
    params_.ef_points.append(cale_robo_param.b3_)
    params_.ef_points.append(cale_robo_param.b4_)

    params_.r_to_cog = cale_robo_param.r_to_cog_
    params_.f_g = cale_robo_param.f_g_
    params_.g_c = cale_robo_param.g_c_

    fk_result = FKDataOut()
    fkSolver(lc_cat, rtation_init, params_, fk_result)

# Define variables used in the main function
Pulley_a_GT = sf.Vector3(-125.0, -110.0, 48.0)
Pulley_b_GT = sf.Vector3( 125.0, -110.0, 48.0)
Pulley_c_GT = sf.Vector3( 125.0,  110.0, 48.0)
Pulley_d_GT = sf.Vector3(-125.0,  110.0, 48.0)

Pulley_a = sf.V3.symbolic("p_a")
Pulley_b = sf.V3.symbolic("p_b")
Pulley_c = sf.V3.symbolic("p_c")
Pulley_d = sf.V3.symbolic("p_d")

Ee_a = sf.Vector3(-0.21, -0.21, -0.011)
Ee_b = sf.Vector3(0.21, -0.21, -0.011)
Ee_c = sf.Vector3(0.21, 0.21, -0.011)
Ee_d = sf.Vector3(-0.21, 0.21, -0.011)

r_to_cog = sf.Matrix([0, 0, -0.12])

g_c = 0.7100703113867337
f_g = 333.54

cale_robo_param_forward = CableRobotParams(g_c, f_g)
cale_robo_param_forward .setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d)
cale_robo_param_forward .setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d)
cale_robo_param_forward .setCog(r_to_cog)

cale_robo_param_inv = CableRobotParams(g_c, f_g)
cale_robo_param_inv.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d)
cale_robo_param_inv.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d)
cale_robo_param_inv.setCog(r_to_cog)

# -----------------------------------------------------------------------#
for largest_cable in list_names_without_l:
    rot_init= sf.Rot3.symbolic("Rot_init")  
    p_platform = sf.Vector3.symbolic("position_vector")

    # inverseKinematicsSolver(cale_robo_param_inv, p_platform, rot_init, largest_cable)

# -----------------------------------------------------------------------#
lc_cat = sf.Vector4.symbolic("encoder")
rtation_init= sf.Rot3.symbolic("Rot_init")  

forwardKinematicsSolver(cale_robo_param_forward, lc_cat, rtation_init)
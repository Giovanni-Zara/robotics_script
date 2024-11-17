import sympy as sp
import numpy as np

#use .evalf() to compute calculations
# use result.applyfunc(sp.trigsimp) to simplify trigonometric symbolic matrixes

"""
psi = sp.symbols('psi')


nm = sp.Matrix([[1, 0, 0],
            [0,1/3,1/3],
            [0, sp.cos(psi), 2/3]])

rx = SymbolicElementaryRotationMatrix.Rx('theta')

sp.pprint(SymbolicElementaryRotationMatrix.multiply_matrices(nm, rx))

"""

class SymbolicElementaryRotationMatrix:

    """
    ANGLE CAN BE : "theta" "psi" "phi"
    """


    """
    # example of usage 

    Rx = SymbolicElementaryRotationMatrix.Rx('theta')
    Rz = SymbolicElementaryRotationMatrix.Rz('phi')
    Rz1 = SymbolicElementaryRotationMatrix.Rz('psi')


    # Multiply Ry and Rz
    result = SymbolicElementaryRotationMatrix.multiply_matrices(Rz, Rx, Rz1)

    # Display the result
    print(f"\n\n")
    sp.pprint(result)
    print(f"\n\n")


    """

    @staticmethod
    def Rx(angle):  
        angle = sp.symbols(angle)  # Symbolic variable for theta
        return sp.Matrix([
            [1, 0, 0],
            [0, sp.cos(angle), -sp.sin(angle)],
            [0, sp.sin(angle), sp.cos(angle)]
        ])
    

    
    @staticmethod
    def Ry(angle):  
        phi = sp.symbols(angle)  # Symbolic variable for phi
        return sp.Matrix([
            [sp.cos(angle), 0, sp.sin(angle)],
            [0, 1, 0],
            [-sp.sin(angle), 0, sp.cos(angle)]
        ])
    
    @staticmethod
    def Rz(angle):  
        angle = sp.symbols(angle)  # Symbolic variable for psi
        return sp.Matrix([
            [sp.cos(angle), -sp.sin(angle), 0],
            [sp.sin(angle), sp.cos(angle), 0],
            [0, 0, 1]
        ])
    @staticmethod
    def multiply_matrices(*matrices):
        # Start with the identity matrix of the appropriate size
        result = sp.eye(3)
        # Multiply all matrices in the order they are provided
        for matrix in matrices:
            result = result * matrix
        return result
    

class PositionAndOrientation:

    @staticmethod
    def solve_direct_problem(r, theta):

        """
        Computes the rotation matrix given a vector r and an angle theta.
        Parameters:
        r (sympy Matrix or numpy array): a 3-dimensional vector.
        theta (float or symbolic): the rotation angle in radians, can be numeric or symbolic.
        """

        # NUMERIC EXAMPLE

        """
        # Define numeric vector and angle
        r_numeric = [1, 0, 0]  # A 3-dimensional vector
        theta_numeric = sp.pi / 4  # 45 degrees in radians

        # Compute the rotation matrix
        rotation_matrix_numeric = PositionAndOrientation.solve_direct_problem(r_numeric, theta_numeric)
        print("Rotation Matrix (numeric theta):")
        print(f"\n\n")
        sp.pprint(rotation_matrix_numeric)
        print(f"\n\n")

        """

        # SYMBOLIC EXAMPLE 

        """
        # Define numeric vector and symbolic angle
        r_numeric = [0, 0, 1]  # A 3-dimensional vector
        theta_symbolic = sp.symbols('theta')  # Symbolic angle

        # Compute the rotation matrix with symbolic theta
        rotation_matrix_symbolic = PositionAndOrientation.solve_direct_problem(r_numeric, theta_symbolic)
        print("Rotation Matrix (symbolic theta):")
        print(f"\n\n")
        sp.pprint(rotation_matrix_symbolic)
        print(f"\n\n")

        """

        # Ensure r is a column vector and normalize it if needed
        r = sp.Matrix(r).reshape(3, 1)  # Convert to a column matrix if not already

        # Calculate the norm of r and normalize if necessary
        norm_r = sp.sqrt(r.dot(r))
        if norm_r != 1:
            print("Normalizing vector r")
            r = r / norm_r  # Normalize the vector

        # Define symbolic variables for theta if not numeric
        theta = sp.symbols('theta') if isinstance(theta, str) else theta

        # Identity matrix
        identity_matrix = sp.eye(3)

        # r * r^T (outer product)
        r_transposed = r * r.T

        # Extract elements for the cross-product matrix
        r_x, r_y, r_z = r[0], r[1], r[2]

        # Skew-symmetric cross-product matrix
        s_of_r_matrix = sp.Matrix([
            [0, -r_z, r_y],
            [r_z, 0, -r_x],
            [-r_y, r_x, 0]
        ])

        # Compute the rotation matrix using Rodrigues' rotation formula
        rotation_matrix = (r_transposed +
                           (identity_matrix - r_transposed) * sp.cos(theta) +
                           s_of_r_matrix * sp.sin(theta))

        return rotation_matrix


    @staticmethod
    def solve_inverse_problem(R):
        """
        Computes the rotation axis vector r and the angle theta from a given rotation matrix R.
        
        Parameters:
        R (sympy Matrix or numpy array): A 3x3 rotation matrix, can contain numeric or symbolic values.
        
        Returns:
        tuple: (r, theta) where r is the rotation axis (as a symbolic/numeric array) and theta is the rotation angle in radians.
        """

        # EXAMPLE FOR NUMERICAL PROBLEM

        """
                
        # Define a numeric rotation matrix
        R_numeric = np.array([
            [-1, 0, 0],
            [0, -1 / np.sqrt(2), -1 / np.sqrt(2)],
            [0, -1 / np.sqrt(2), 1 / np.sqrt(2)]
        ])

        # Solve for r and theta
        r_numeric, theta_numeric = PositionAndOrientation.solve_inverse_problem(R_numeric)

        """

        # EXAMPLE FOR SYMBOLIC PROBLEM 

        """
        
        # Define a symbolic rotation matrix
        phi = sp.symbols('phi')
        R_symbolic = sp.Matrix([
            [sp.cos(phi), -sp.sin(phi), 0],
            [sp.sin(phi), sp.cos(phi), 0],
            [0, 0, 1]
        ])

        # Solve for r and theta
        r_symbolic, theta_symbolic = PositionAndOrientation.solve_inverse_problem(R_symbolic)
        
        """

        # Convert numpy array to sympy Matrix if necessary
        if isinstance(R, np.ndarray):
            R = sp.Matrix(R)

        # Extract elements from R for symbolic computation
        R_1_2, R_2_1 = R[0, 1], R[1, 0]
        R_1_3, R_3_1 = R[0, 2], R[2, 0]
        R_2_3, R_3_2 = R[1, 2], R[2, 1]
        R_1_1, R_2_2, R_3_3 = R[0, 0], R[1, 1], R[2, 2]

        # Compute symbolic expressions for the components
        positive_arg_inside_atan2 = sp.sqrt((R_1_2 - R_2_1)**2 + (R_1_3 - R_3_1)**2 + (R_2_3 - R_3_2)**2)
        trace_minus_one = R_1_1 + R_2_2 + R_3_3 - 1

        # Compute theta options symbolically
        theta_of_pos = sp.atan2(positive_arg_inside_atan2, trace_minus_one)
        theta_of_neg = sp.atan2(-positive_arg_inside_atan2, trace_minus_one)

        
        # Choose theta
        theta = (theta_of_neg, theta_of_pos)
        print("This is theta before:", theta)

        if theta_of_neg == 0 and theta_of_pos == 0:
            print("Singular case, no unique solution")
            return None
        
        elif 0 in theta:
            # Choose the non-zero theta value
            theta = theta_of_pos if theta_of_pos != 0 else theta_of_neg
            array_to_compute_r = sp.Matrix([R_3_2 - R_2_3, R_1_3 - R_3_1, R_2_1 - R_1_2])
            r = (1 / (2 * sp.sin(theta))) * array_to_compute_r
       
        elif sp.pi in theta:

            print("special case, theta is pi, set sin 𝜃 = 0, cos 𝜃 = −1 and solve for...")
            theta = sp.pi  # Set theta to pi explicitly
                # Compute magnitudes of r components
            rx = sp.sqrt((R_1_1 + 1) / 2)
            ry = sp.sqrt((R_2_2 + 1) / 2)
            rz = sp.sqrt((R_3_3 + 1) / 2)
            
            # Determine signs based on the off-diagonal elements of R
            # Using rx * ry = R[1, 2] / 2
            if (R_1_2 / 2).is_real:
                sign_rx_ry = sp.sign(R_1_2)
                rx = sign_rx_ry * rx if rx != 0 else rx
                ry = sign_rx_ry * ry if ry != 0 else ry
            
            # Using rx * rz = R[2, 1] / 2
            if (R_2_1 / 2).is_real:
                sign_rx_rz = sp.sign(R_1_3)
                rx = sign_rx_rz * rx if rx != 0 else rx
                rz = sign_rx_rz * rz if rz != 0 else rz
            
            # Using ry * rz = R[1, 3] / 2
            if (R_1_3 / 2).is_real:
                sign_ry_rz = sp.sign(R_1_3)
                ry = sign_ry_rz * ry if ry != 0 else ry
                rz = sign_ry_rz * rz if rz != 0 else rz

            # Construct the final vector

            print(f"\n this is the result after using the relationship to solve sign ambiguities")
            r1 = sp.Matrix([rx, ry, rz])
            r2 =  sp.Matrix([-rx, -ry, -rz])

            r= (r1,r2)
                
        
        else:
            theta = theta_of_pos  # Choose the positive theta as the primary solution
            array_to_compute_r = sp.Matrix([R_3_2 - R_2_3, R_1_3 - R_3_1, R_2_1 - R_1_2])
            
            # Two possible solutions for r
            rx = (1 / (2 * sp.sin(theta_of_neg))) * array_to_compute_r
            ry = (1 / (2 * sp.sin(theta_of_pos))) * array_to_compute_r
            r = (rx, ry)  # Return both possible solutions for r

        print("r is:")
        print(f"\n")
        sp.pprint(r)
        print(f"\n")
        print("THETA IS IN RAD")
        print(f"\n")
        sp.pprint(theta)
        print(f"\n")

        return r, theta
    


class EulerAngles:
    @staticmethod
    def solve_direct_problem(rotation_1, rotation_2, rotation_3, first_angle, second_angle, third_angle):
        """
        Computes the rotation matrix given three consecutive rotations of different angles.
        The rotation sequence must be valid, e.g., [Rx, Rz, Ry] or [Ry, Rz, Rx], but not [Rx, Rx, Ry].
        
        Supports both symbolic and numeric angles for each rotation.
        """

        # EXAMPLE USAGE SYMBOLIC FOR (ALL)

        """
        
        result = EulerAngles.solve_direct_problem('Rz', 'Rx', 'Rz', 'phi', 'theta', 'psi')
        print(f"\n")
        sp.pprint(result)  # This will print the matrix, showing mixed symbolic and numeric values
        print(f"\n")

        """

        # EXAMPLE USAGE PARTIAL NUMERIC FORM

        """
        
        result = EulerAngles.solve_direct_problem('Rz', 'Rx', 'Rz', sp.pi/2, 'theta', 'psi')
        print(f"\n")
        sp.pprint(result)  # This will print the matrix, showing mixed symbolic and numeric values
        print(f"\n")
        
        
        """
        
        # Define rotation mapping to symbolic or numeric matrices as needed
        def get_rotation_matrix(axis, angle):
            if isinstance(angle, str):  # If angle is a symbolic string
                return getattr(SymbolicElementaryRotationMatrix, axis)(angle)
            else:  # Otherwise, assume angle is numeric
                return getattr(ElementaryRotationMatrix, axis)(angle)

     
        
        # Get the rotation matrices based on the specified sequence and angles
        matrix_1 = get_rotation_matrix(rotation_1, first_angle)
        matrix_2 = get_rotation_matrix(rotation_2, second_angle)
        matrix_3 = get_rotation_matrix(rotation_3, third_angle)
        
        # Compute the final rotation matrix by multiplying in the given sequence
        final_rotation_matrix = SymbolicElementaryRotationMatrix.multiply_matrices(matrix_1, matrix_2, matrix_3)
        
        return final_rotation_matrix
    
class RollPitchYawAngles:
    @staticmethod
    def solve_direct_problem(rotation_1, rotation_2, rotation_3, first_angle, second_angle, third_angle):
        """
        Computes the rotation matrix given three consecutive rotations of different angles.
        The rotation sequence must be valid, e.g., [Rx, Rz, Ry] or [Ry, Rz, Rx], but not [Rx, Rx, Ry].
        
        Supports both symbolic and numeric angles for each rotation.
        """

        # EXAMPLE USAGE SYMBOLIC FOR (ALL)

        """
        
        result = RollPitchYawAngles.solve_direct_problem('Rx', 'Ry', 'Rz', 'psi', 'theta', 'phi')
        print(f"\n")
        sp.pprint(result)  # This will print the matrix, showing mixed symbolic and numeric values
        print(f"\n")

        """

        # EXAMPLE USAGE PARTIAL NUMERIC FORM

        """
        
        result = EulerAngles.solve_direct_problem('Rz', 'Rx', 'Rz', sp.pi/2, 'theta', 'psi')
        print(f"\n")
        sp.pprint(result)  # This will print the matrix, showing mixed symbolic and numeric values
        print(f"\n")
        
        
        """
        
        # Define rotation mapping to symbolic or numeric matrices as needed
        def get_rotation_matrix(axis, angle):
            if isinstance(angle, str):  # If angle is a symbolic string
                return getattr(SymbolicElementaryRotationMatrix, axis)(angle)
            else:  # Otherwise, assume angle is numeric
                return getattr(ElementaryRotationMatrix, axis)(angle)

     
        
        # Get the rotation matrices based on the specified sequence and angles
        matrix_1 = get_rotation_matrix(rotation_3, third_angle)
        matrix_2 = get_rotation_matrix(rotation_2, second_angle)
        matrix_3 = get_rotation_matrix(rotation_1, first_angle)
        
        # Compute the final rotation matrix by multiplying in the given sequence
        final_rotation_matrix = SymbolicElementaryRotationMatrix.multiply_matrices(matrix_1, matrix_2, matrix_3)
        
        return final_rotation_matrix



class SerialRobotDH:
    def __init__(self, dhtable):
        """
        Initialize the robot with a given DH parameter table.
        The table should be of the form:
        [
            [alpha1, a1, d1, theta1],
            [alpha2, a2, d2, theta2],
            ...
        ]
        """

        # Example usage:
        """
        # Define symbolic variables

        q1, q2, q3, q4 = sp.symbols('q1 q2 q3 q4')
        a1, a2, d1, d4 = sp.symbols('a1 a2 d1 d4')

        # Define DH table for SCARA robot
        DHTABLE = [
            [0,     a1, d1, q1],
            [0,     a2, 0,  q2],
            [0,     0,  q3, 0],
            [sp.pi, 0,  d4, q4]
        ]

        # Create an instance of the SerialRobotDH class
        scara_robot = SerialRobotDH(DHTABLE)

        # Compute symbolic forward kinematics
        T_symbolic = scara_robot.forward_kinematics()
        print("Symbolic Transformation Matrix (T0N):")
        sp.pprint(T_symbolic)

        # Extract position and orientation
        position, x_axis, y_axis, z_axis = scara_robot.get_position_orientation(T_symbolic)
        print("\nPosition (End-Effector):")
        sp.pprint(position)
        print("\nxN Axis:")
        sp.pprint(x_axis)
        print("\nyN Axis:")
        sp.pprint(y_axis)
        print("\nzN Axis:")
        sp.pprint(z_axis)

        # Compute numeric forward kinematics with specific values
        numeric_result = scara_robot.numeric_forward_kinematics(q1=0.5, q2=0.3, q3=0.2, q4=0.1, a1=1, a2=1, d1=0.5, d4=0.1)
        print("\nNumeric Transformation Matrix (T0N):")
        sp.pprint(numeric_result)
        """
        self.dhtable = dhtable
        self.N = len(dhtable)  # Number of joints
        
    def dh_transform(self, alpha, a, d, theta):
        """
        Compute the DH transformation matrix using the provided parameters.
        """
        return sp.Matrix([
            [sp.cos(theta), -sp.sin(theta)*sp.cos(alpha),  sp.sin(theta)*sp.sin(alpha), a*sp.cos(theta)],
            [sp.sin(theta),  sp.cos(theta)*sp.cos(alpha), -sp.cos(theta)*sp.sin(alpha), a*sp.sin(theta)],
            [0,              sp.sin(alpha),               sp.cos(alpha),              d],
            [0,              0,                           0,                          1]
        ])
    
    def forward_kinematics(self):
        """
        Compute the forward kinematics (end-effector position and orientation) by multiplying the
        individual transformation matrices.
        """
        T = sp.eye(4)  # Initialize with identity matrix
        
        # Compute the transformation matrix for each link
        for i in range(self.N):
            alpha, a, d, theta = self.dhtable[i]
            A_i = self.dh_transform(alpha, a, d, theta)
            T = T * A_i
            T = sp.simplify(T)
        
        return T
    
    def get_position_orientation(self, T):
        """
        Extract position (x, y, z) and orientation vectors (xN, yN, zN axes) from the transformation matrix.
        """
        # Extract position vector
        p = T[0:3, 3]
        
        # Extract orientation vectors (rotation axes)
        n = T[0:3, 0]  # xN axis
        s = T[0:3, 1]  # yN axis
        a = T[0:3, 2]  # zN axis
        
        return p, n, s, a
    
    def numeric_forward_kinematics(self, **kwargs):
        """
        Compute numeric forward kinematics using provided values for the DH parameters.
        kwargs should include values for symbols used in the DH table.
        """
        T = self.forward_kinematics()
        return T.evalf(subs=kwargs)


import math as m

q1, q2, q3, q4 = sp.symbols('q1 q2 q3 q4')
a1, a2, d1, d4 = sp.symbols('a1 a2 d1 d4')

# Define DH table for SCARA robot
DHTABLE = [
[sp.pi,   -0.5, 0, q1],
[-sp.pi/2,  0.6, 0,  q2]
]
scara_robot = SerialRobotDH(DHTABLE)
numeric_result = scara_robot.numeric_forward_kinematics(q1=sp.pi/2, q2=-sp.pi/2)
print("\nNumeric Transformation Matrix (T0N):")
sp.pprint(numeric_result)
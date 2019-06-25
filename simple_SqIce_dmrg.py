#!/usr/bin/env python
#
# Simple DMRG for square ice
#  - Infinite system algorithm
#  - Finite system algorithm
#
# This code is based on the DMRG Tutorial of:
# Copyright 2013 James R. Garrison and Ryan V. Mishmash.
# Open source under the MIT license.  Source code at
# <https://github.com/simple-dmrg/simple-dmrg/>

# This code will run under any version of Python >= 2.6.  The following line
# provides consistency between python2 and python3.
from __future__ import print_function, division  # requires Python >= 2.6

# numpy and scipy imports
import numpy as np
from scipy.sparse import kron, identity
from scipy.sparse.linalg import eigsh  # Lanczos routine from ARPACK

# We will use python's "namedtuple" to represent the Block and EnlargedBlock
# objects
from collections import namedtuple

Block = namedtuple("Block", ["length", "basis_size", "operator_dict"])
EnlargedBlock = namedtuple("EnlargedBlock", ["length", "basis_size", "operator_dict"])

def is_valid_block(block):
    for op in block.operator_dict.values():
        # here we have to differentiate between an operator list and a
        # single operator
        if type(op) is list:
            for op_element in op:
                if op_element.shape[0] != block.basis_size or op_element.shape[1] != block.basis_size:
                    return False
        else:
            if op.shape[0] != block.basis_size or op.shape[1] != block.basis_size:
                return False
    return True

def mkron( *kron_list):
    A = kron_list[0]
    for i,val in  enumerate(kron_list[1:],start=1):
        A = np.kron(A, val)  
    
    return A

# This function should test the same exact things, so there is no need to
# repeat its definition.
is_valid_enlarged_block = is_valid_block

# Model-specific code for 2 x Lx chain
model_d = 2**6  # single-site basis size

I2 = np.eye(2)

Sp1 = np.array([[0, 1], [0, 0]], dtype='d')  # single-site S^+
Sm1 = Sp1.conjugate().transpose()
Sz1 = np.array([[0.5, 0], [0, -0.5]], dtype='d')  # single-site S^z

Plq_right = [ mkron(I2,I2,Sp1,I2,Sm1,Sp1), mkron(I2,I2,I2,Sp1,Sp1,Sm1)] #right half plaquette
Plq_right_deg = [mkron(I2,I2,Sm1,I2,Sp1,Sm1),mkron(I2,I2,I2,Sm1,Sm1,Sp1)] #right half plaquette

Plq_left = [-mkron(Sp1,Sm1,Sp1,I2,I2,I2),-mkron(Sm1,Sp1,I2,Sp1,I2,I2)]
Plq_left_deg = [-mkron(Sm1,Sp1,Sm1,I2,I2,I2), -mkron(Sp1,Sm1,I2,Sm1,I2,I2)]

H1 = np.zeros([2**6, 2**6] , dtype='d')  # single-site portion of H is zero

def H2(hleft_list, hright_list):  # two-site part of H
    """Given the operators S^z and S^+ on two sites in different Hilbert spaces
    (e.g. two blocks), returns a Kronecker product representing the
    corresponding two-site term in the Hamiltonian that joins the two sites.
    """
    Hinteract = 0
    for hleft,hright in zip(hleft_list,hright_list):    
        Hinteract += kron(hleft,hright) 
    
    J = Jz = 1.
    Hinteract = J * Hinteract
    return Hinteract

# conn refers to the connection operator, that is, the operator on the edge of
# the block, on the interior of the chain.  We need to be able to represent S^z
# and S^+ on that site in the current basis in order to grow the chain.
initial_block_left = Block(length=1, basis_size=model_d, operator_dict={
    "H": H1,
    "h1loc": Plq_right, #h1n1 h1n2
    "h2loc": Plq_right_deg, #h2n1 h2n2
})

initial_block_right = Block(length=1, basis_size=model_d, operator_dict={
    "H": H1,
    "h1loc": Plq_left,
    "h2loc": Plq_left_deg,
})


def enlarge_block(block, side="left"):
    """This function enlarges the provided Block by a single site, returning an
    EnlargedBlock.
    """
    mblock = block.basis_size
    o = block.operator_dict

    # Create the new operators for the enlarged block.  Our basis becomes a
    # Kronecker product of the Block basis and the single-site basis.  NOTE:
    # `kron` uses the tensor product convention making blocks of the second
    # array scaled by the first.  As such, we adopt this convention for
    # Kronecker products throughout the code.
    
    if side == "left":
        left_block = o["h1loc"] + o["h1loc"]
        right_block = Plq_right + Plq_right_deg
        
        enlarged_operator_dict = {
                "H": kron(o["H"], identity(model_d)) + kron(identity(mblock), H1) + H2(right_block, left_block),
                "h1loc": [kron(identity(mblock), op) for op in Plq_left],
                "h2loc": [kron(identity(mblock), op) for op in Plq_left_deg],
                }
    else:
        right_block = o["h1loc"] + o["h1loc"]
        left_block = Plq_left + Plq_left_deg
        
        enlarged_operator_dict = {
                "H": kron(o["H"], identity(model_d)) + kron(identity(mblock), H1) + H2(right_block, left_block),
                "h1loc": [kron(op, identity(mblock)) for op in Plq_right],
                "h2loc": [kron(op, identity(mblock)) for op in Plq_right_deg],

                }
        

    return EnlargedBlock(length=(block.length + 1),
                         basis_size=(block.basis_size * model_d),
                         operator_dict=enlarged_operator_dict)

def rotate_and_truncate(operator, transformation_matrix):
    """Transforms the operator to the new (possibly truncated) basis given by
    `transformation_matrix`.
    """
    return transformation_matrix.conjugate().transpose().dot(operator.dot(transformation_matrix))

def single_dmrg_step(sys, env, m):
    """Performs a single DMRG step using `sys` as the system and `env` as the
    environment, keeping a maximum of `m` states in the new basis.
    """
    assert is_valid_block(sys)
    assert is_valid_block(env)

    # Enlarge each block by a single site.
    sys_enl = enlarge_block(sys,'left')
    if sys is env:  # no need to recalculate a second time
        env_enl = sys_enl
    else:
        env_enl = enlarge_block(env,'right')

    assert is_valid_enlarged_block(sys_enl)
    assert is_valid_enlarged_block(env_enl)

    # Construct the full superblock Hamiltonian.
    m_sys_enl = sys_enl.basis_size
    m_env_enl = env_enl.basis_size
    sys_enl_op = sys_enl.operator_dict
    env_enl_op = env_enl.operator_dict
    superblock_hamiltonian = kron(sys_enl_op["H"], identity(m_env_enl)) + kron(identity(m_sys_enl), env_enl_op["H"]) + \
                             H2(sys_enl_op["h1loc"] + sys_enl_op["h2loc"], env_enl_op["h1loc"] + env_enl_op["h2loc"])

    # Call ARPACK to find the superblock ground state.  ("SA" means find the
    # "smallest in amplitude" eigenvalue.)
    (energy,), psi0 = eigsh(superblock_hamiltonian, k=1, which="SA")

    # Construct the reduced density matrix of the system by tracing out the
    # environment
    #
    # We want to make the (sys, env) indices correspond to (row, column) of a
    # matrix, respectively.  Since the environment (column) index updates most
    # quickly in our Kronecker product structure, psi0 is thus row-major ("C
    # style").
    psi0 = psi0.reshape([sys_enl.basis_size, -1], order="C")
    rho = np.dot(psi0, psi0.conjugate().transpose())

    # Diagonalize the reduced density matrix and sort the eigenvectors by
    # eigenvalue.
    evals, evecs = np.linalg.eigh(rho)
    possible_eigenstates = []
    for eval, evec in zip(evals, evecs.transpose()):
        possible_eigenstates.append((eval, evec))
    possible_eigenstates.sort(reverse=True, key=lambda x: x[0])  # largest eigenvalue first

    # Build the transformation matrix from the `m` overall most significant
    # eigenvectors.
    my_m = min(len(possible_eigenstates), m)
    transformation_matrix = np.zeros((sys_enl.basis_size, my_m), dtype='d', order='F')
    for i, (eval, evec) in enumerate(possible_eigenstates[:my_m]):
        transformation_matrix[:, i] = evec

    truncation_error = 1 - sum([x[0] for x in possible_eigenstates[:my_m]])
    print("truncation error:", truncation_error)

    # Rotate and truncate each operator.
    new_operator_dict = {}
    for name, op in sys_enl.operator_dict.items():
        if type(op) is list:
            new_operator_dict[name] = [rotate_and_truncate(o, transformation_matrix) for o in op]
        else:
            new_operator_dict[name] = rotate_and_truncate(op, transformation_matrix)

    new_left_block = Block(length=sys_enl.length,
                     basis_size=my_m,
                     operator_dict=new_operator_dict)
    
        # Rotate and truncate each operator.
    new_operator_dict = {}
    for name, op in sys_enl.operator_dict.items():
        if type(op) is list:
            new_operator_dict[name] = [rotate_and_truncate(o, transformation_matrix) for o in op]
        else:
            new_operator_dict[name] = rotate_and_truncate(op, transformation_matrix)

    new_right_block = Block(length=sys_enl.length,
                     basis_size=my_m,
                     operator_dict=new_operator_dict)

    return new_left_block, new_right_block, energy

def graphic(sys_block, env_block, sys_label="l"):
    """Returns a graphical representation of the DMRG step we are about to
    perform, using '=' to represent the system sites, '-' to represent the
    environment sites, and '**' to represent the two intermediate sites.
    """
    assert sys_label in ("l", "r")
    graphic = ("=" * sys_block.length) + "**" + ("-" * env_block.length)
    if sys_label == "r":
        # The system should be on the right and the environment should be on
        # the left, so reverse the graphic.
        graphic = graphic[::-1]
    return graphic

def infinite_system_algorithm(L, m):
    block = initial_block
    # Repeatedly enlarge the system by performing a single DMRG step, using a
    # reflection of the current block as the environment.
    while 2 * block.length < L:
        print("L =", block.length * 2 + 2)
        block, energy = single_dmrg_step(block, block, m=m)
        print("E/L =", energy / (block.length * 2))

def finite_system_algorithm(L, m_warmup, m_sweep_list):
    assert L % 2 == 0  # require that L is an even number

    # To keep things simple, this dictionary is not actually saved to disk, but
    # we use it to represent persistent storage.
    block_disk = {}  # "disk" storage for Block objects

    # Use the infinite system algorithm to build up to desired size.  Each time
    # we construct a block, we save it for future reference as both a left
    # ("l") and right ("r") block, as the infinite system algorithm assumes the
    # environment is a mirror image of the system.
    right_block = initial_block_right
    left_block = initial_block_left
    block_disk["l", right_block.length] = right_block
    block_disk["r", left_block.length] = left_block
    while left_block.length + right_block.length < L:
        # Perform a single DMRG step and save the new Block to "disk"
        print(graphic(left_block, right_block))
        left_block, right_block, energy = single_dmrg_step(left_block, right_block, m=m_warmup)
        print("E/L =", energy / (left_block.length * 2))
        block_disk["l", left_block.length] = left_block
        block_disk["r", right_block.length] = right_block

    # Now that the system is built up to its full size, we perform sweeps using
    # the finite system algorithm.  At first the left block will act as the
    # system, growing at the expense of the right block (the environment), but
    # once we come to the end of the chain these roles will be reversed.
    sys_label, env_label = "l", "r"
    sys_block = left_block; del left_block; del right_block  # rename the variable
    for m in m_sweep_list:
        while True:
            # Load the appropriate environment block from "disk"
            env_block = block_disk[env_label, L - sys_block.length - 2]
            if env_block.length == 1:
                # We've come to the end of the chain, so we reverse course.
                sys_block, env_block = env_block, sys_block
                sys_label, env_label = env_label, sys_label

            # Perform a single DMRG step.
            print(graphic(sys_block, env_block, sys_label))
            sys_block, env_block, energy = single_dmrg_step(sys_block, env_block, m=m)

            print("E/L =", energy / L)

            # Save the block from this step to disk.
            block_disk[sys_label, sys_block.length] = sys_block

            # Check whether we just completed a full sweep.
            if sys_label == "l" and 2 * sys_block.length == L:
                break  # escape from the "while True" loop

if __name__ == "__main__":
    np.set_printoptions(precision=10, suppress=True, threshold=10000, linewidth=300)

    #infinite_system_algorithm(L=100, m=20)
    finite_system_algorithm(L=10, m_warmup=2, m_sweep_list=[10, 20, 30, 40, 40])
# -*- coding: utf-8 -*-


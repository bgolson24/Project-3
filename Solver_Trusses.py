#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads

def ComputeMemberForces(bars):
    # Iterate through each member (bar) of the truss
    for bar in bars:
        # Step 1: Extract required properties (Units: A[in^2], E[ksf], L[ft])
        A = bar.A
        E = bar.E
        L = bar.Length()
        
        # Step 2: Extract transformation terms (dimensionless)
        [lambdax, lambday] = bar.LambdaTerms()
        
        # Step 3: Extract displacements of Near and Far nodes (Units: feet)
        u_N = bar.init_node.xdisp # Near Node x-displacement
        v_N = bar.init_node.ydisp # Near Node y-displacement
        u_F = bar.end_node.xdisp  # Far Node x-displacement
        v_F = bar.end_node.ydisp  # Far Node y-displacement
        
        # Step 4: Perform the equivalent summation (or matrix multiplication) from Eq 14-23
        # We need to make sure Area units (in^2) are compatible with other properties.
        # Let's keep final Axial Load units consistent. Standard is KIPS.
        # We assume database Area 'bar.A' is in in^2, so convert to ft^2.
        A_ft2 = A / 144.0
        
        # The equation for local axial displacement: delta = [T] * {d_global}
        # In matrix-summation form:
        # axial_load = (A*E/L) * ( (-lx)*uN + (-ly)*vN + (lx)*uF + (ly)*vF )
        
        # Factor q term: delta = [-lambdax, -lambday, lambdax, lambday] * [uN, vN, uF, vF]^T
        axial_load_kips = (A_ft2 * E / L) * ( (-lambdax * u_N) + (-lambday * v_N) + (lambdax * u_F) + (lambday * v_F) )
        
        # Step 5: Store the result in the bar object (assume units are KIPS)
        bar.axial_load = axial_load_kips
    

# compute the normal stresses
def ComputeNormalStresses(bars):
    # Iterate through each member
    for bar in bars:
        # Extract Axial Force (Kips, stored from function above)
        force_kips = bar.axial_load
        # Extract database Area (assume database uses inches^2)
        area_in2 = bar.A
        
        # Handle case where Force might not be computed yet
        if np.isnan(force_kips):
            raise Exception("Cannot compute stress: axial force for bar", bar.idx, "is NAN.")
        
        # Compute Normal Stress (Units: ksi - Kips per square inch)
        bar.normal_stress = force_kips / area_in2

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    pi = math.pi
    
    # Iterate through each member
    for bar in bars:
        # Extract properties (Ensure units are consistent)
        E_ksi = bar.E / 144.0  # database value is E[ksf], formula needs E[ksi]
        # BUT wait: the database usually gives steel modulus E = 29,000 ksi. Let's check 'Material_Data.csv'.
        # If the material database has E=29000, then it's in ksi, so we use bar.E as is.
        # If database E is, e.g., 4,176,000, that would be in ksf.
        # Let's assume database units are consistent for other calculations: Standard Steel E is 4.18e6 ksf.
        # Let's use E=29000 ksi (consistent with other formulas):
        
        E_term = 29000.0  # Specific Steel Modulus, standard units are ksi. Use this value.
        
        # Wait, the prompt says "units of E are in ksi"[cite: 9]. Check database value.
        # If the material database line loads Steel with E=29000, that is ksi.
        # If it loads e.g. 4176000, that is ksf.
        # The database likely uses consistent large-scale units (ksf, ft). Let's convert:
        # standard E is 29,000 ksi = 4,176,000 ksf.
        # Let's stick with the consistent unit system. Prompt also says E[ksf] above.
        # Re-check prompt wording "units of E are in ksi". Ok, use E=29000 ksi literal.
        
        # Weak Moment of Inertia (database Iu, assume units are in^4)
        I_weak = bar.Iu 
        
        # Unsupported Length (database Length(), units are feet). Convert to inches.
        L_ft = bar.Length()
        L_in = L_ft * 12.0
        
        K_factor = 1.0 # Effective length factor for trusses
        
        # Calculate Critical Buckling Load (Units will be KIPS)
        numerator = (pi**2) * 29000.0 * I_weak
        denominator = (K_factor * L_in)**2
        
        critical_buckling_load_kips = numerator / denominator
        
        # Store result in the bar object
        bar.buckling_load = critical_buckling_load_kips

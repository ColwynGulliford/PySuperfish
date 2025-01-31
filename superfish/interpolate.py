from superfish.parsers import parse_fish_t7, parse_poisson_t7

from pmd_beamphysics import FieldMesh

import numpy as np

from glob import glob
import os


def get_t7(path):
    """
    Returns a list of T7 files in path
    """
    return glob(os.path.join(path, '*T7'))


def interpolate2d(sf,
                  x1min=-1000, x1max=1000, nx1=100,
                  x2min=-1000, x2max=1000, nx2=100,
                  return_fieldmesh=False):
    """
    
    Runs SF7.EXE on a Superfish object sf, requesting interpolating data, and reads the Parmela T7 file. 
    
    Units are in the input units of the program. 
    
    Labels will label the output columns. The user must know what these are from the problem.
    
    TODO: detect these from SF object
    
    SF7 automatically adjusts the bounds if they are requested outside of the computational domain.
    
    # TODO write up the return
    """

    if sf.geometry == 'cylindrical':
        return interpolate_cylindrical(sf,
                                       zmin=x1min, zmax=x1max, nz=nx1,
                                       rmin=x2min, rmax=x2max, nr=nx2, 
                                       return_fieldmesh=return_fieldmesh)
        
    elif sf.geometry == 'rectangular':
        interpolate_xy(sf, 
                   xmin=x1min, xmax=x1max, nx=nx1,
                   ymin=x2max, ymax=x2max, ny=nx2,
                   return_fieldmesh=False)

        
    

def interpolate_xy(sf, 
                   xmin=-1000, xmax=1000, nx=100,
                   ymin=-1000, ymax=1000, ny=200,
                   return_fieldmesh=False):

    problem = sf.problem 
    
    # fish and poisson have the opposite conventions:
    if problem == 'fish':
    # The interpolation grid
        F=f"""
{xmin} {ymin} {xmax} {ymax}
{nx-1} {ny-1}
End"""
        
    elif problem == 'poisson':
        F=f"""
{ymin} {xmin} {ymax} {xmax}
{ny-1} {nx-1}
End"""
    else:
        raise ValueError(f'unknowm problem: {problem}')
    

    # Clear old T7
    for f in get_t7(sf.path):
        os.remove(f)
    
    # Write 
    ifile = sf.basename+'.IN7'
    with open(os.path.join(sf.path, ifile), 'w') as f:
        f.write(F)
    
    # Needed so that the fields aren't normalized to 1 MV/m average
    inifile = os.path.join(sf.path, 'SF.INI')
    with open(inifile, 'w') as f:
        f.write("""[global]
Force1MVperMeter=No""")

    # Needed on WSL, otherwise optional
    t35file = sf.basename+'.T35'
    
        
    # Run
    sf.run_cmd('sf7', ifile, t35file)

    return # Right now just run
    
    # Get the filename
    t7file = get_t7(sf.path)
    assert len(t7file) == 1, f'T7 file is missing.'
    t7file = t7file[0]
    
    
    # Optional fieldmesh parsing
    if return_fieldmesh:
        # Parsing is different for each:
        if problem == 'fish':
            type = 'electric'
            return FieldMesh.from_superfish(t7file)
        else:
            if sf.output['sfo']['header']['variable']['XJFACT'] == 0:
                type = 'electric'
            else:
                type = 'magnetic'
            return FieldMesh.from_superfish(t7file, type=type)
        
        
    
    # Parsing is different for each:
    if problem == 'fish':
        type = 'electric'
        d =  parse_fish_t7(t7file)
    else:
        if sf.output['sfo']['header']['variable']['XJFACT'] == 0:
            type = 'electric'
        else:
            type = 'magnetic'
        d =  parse_poisson_t7(t7file, 'cylindrical', type=type)
    
    return d


def interpolate_cylindrical(sf,
                  zmin=-1000, zmax=1000, nz=100,
                  rmin=0, rmax=0, nr=1, return_fieldmesh=False):
    """
    
    Runs SF7.EXE on a Superfish object sf, requesting interpolating data, and reads the Parmela T7 file. 
    
    Units are in the input units of the program. 
    
    Labels will label the output columns. The user must know what these are from the problem.
    
    TODO: detect these from SF object
    
    SF7 automatically adjusts the bounds if they are requested outside of the computational domain.
    
    Returns a dict with:
        rmin: minimum radius in cm
        rmax: maximum radius in cm
        nr:   number of radius points
        zmin: minimum z in cm
        zmax: maximum z in cm
        nz:   number of z points
        freq: frequency in MHz
        data: 2D array of shape (nr, nz)
    
    
    """
    
    problem = sf.problem 
    
    # fish and poisson have the opposite conventions:
    if problem == 'fish':
    # The interpolation grid
        F=f"""Parmela
{zmin} {rmin} {zmax} {rmax}
{nz-1} {nr-1}
End"""
        
    elif problem == 'poisson':
        F=f"""Parmela
{rmin} {zmin} {rmax} {zmax}
{nr-1} {nz-1}
End"""
    else:
        raise ValueError(f'unknowm problem: {problem}')
    

    # Clear old T7
    for f in get_t7(sf.path):
        os.remove(f)
    
    # Write 
    ifile = sf.basename+'.IN7'
    with open(os.path.join(sf.path, ifile), 'w') as f:
        f.write(F)
    
    # Needed so that the fields aren't normalized to 1 MV/m average
    inifile = os.path.join(sf.path, 'SF.INI')
    with open(inifile, 'w') as f:
        f.write("""[global]
Force1MVperMeter=No""")

    # Needed on WSL, otherwise optional
    t35file = sf.basename+'.T35'
    
        
    # Run
    sf.run_cmd('sf7', ifile, t35file)
    
    # Get the filename
    t7file = get_t7(sf.path)
    assert len(t7file) == 1, f'T7 file is missing.'
    t7file = t7file[0]
    
    
    # Optional fieldmesh parsing
    if return_fieldmesh:
        # Parsing is different for each:
        if problem == 'fish':
            type = 'electric'
            return FieldMesh.from_superfish(t7file)
        else:
            if sf.output['sfo']['header']['variable']['XJFACT'] == 0:
                type = 'electric'
            else:
                type = 'magnetic'
            return FieldMesh.from_superfish(t7file, type=type)
        
        
    
    # Parsing is different for each:
    if problem == 'fish':
        type = 'electric'
        d =  parse_fish_t7(t7file)
    else:
        if sf.output['sfo']['header']['variable']['XJFACT'] == 0:
            type = 'electric'
        else:
            type = 'magnetic'
        d =  parse_poisson_t7(t7file, 'cylindrical', type=type)
    
    return d

    

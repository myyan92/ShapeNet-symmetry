import os, sys
import numpy as np
import glob
import pdb

def rotateAxis(axis, rotate_axis, degree):
    """Rotate vector axis around rotate_axis for degree

    Args:
        axis: A (3,) numpy array. The vector to be rotated.
        rotate_axis: A (3,) numpy array. The vector to rotate around.
        degree: Float. In radians, the angle to rotate.
    Return:
        axis: A (3,) numpy array, the rotated vector.
    """
    third_axis = np.cross(rotate_axis, axis)
    third_axis /= np.linalg.norm(third_axis)
    rotated_axis = np.cos(degree)*axis + np.sin(degree)*third_axis
    rotated_axis /= np.linalg.norm(rotated_axis)
    return rotated_axis

def printf(file, type, center, axis, degree=None):
    if degree == 20:
        degree = 0
    if type=='C':
        file.write('C, %f %f %f, %f %f %f, %d\n' % (
            center[0], center[1], center[2],
            axis[0], axis[1], axis[2], degree))
    else:
        file.write('R, %f %f %f, %f %f %f\n' % (
            center[0], center[1], center[2],
            axis[0], axis[1], axis[2]))
    return

def convert(input, output):
    """Convert symmetry representation from sym3 format to a list of 
        reflection / rotation axes.

    Args:
        input: A string, the filename ending with ".sym3".
        output: A string, the filename to write output.
    Raises:
        NotImplementedError: If symmetry type is 'D*d', 'T' or 'I'.
    """
    with open(input) as f:
        lines = f.readlines()

    sym_type = lines[0].strip()
    center = lines[1].strip().split()
    center = [float(c) for c in center]
    center = -np.array(center)

    matrix = [l.strip().split() for l in lines[2:5]]
    matrix = [[float(e) for e in row] for row in matrix]
    matrix = np.array(matrix)

    if sym_type == 'E':
        return
    if sym_type == 'Cs':
        with open(output, 'w') as f:
            printf(f, 'R', center, matrix[:,2])
        return
    if sym_type == 'O':
        convertOGroup(center, matrix, output)
        return
    if sym_type == 'O3':
        with open(output, 'w') as f:
            printf(f, 'C', center, np.zeros((3,)), 0)
        return
    if (sym_type in ['T', 'I'] or 
        (sym_type[0]=='D' and sym_type[-1]=='d')):
        raise NotImplementedError(
            'Symmetry type T, O, and D*d are not implemented.')
    fout = open(output, 'w')
    if sym_type[-1] in ['v', 'h', 'd']:
        degree = int(sym_type[1:-1])
    else:
        degree = int(sym_type[1:])
    printf(fout, 'C', center, matrix[:,1], degree)
    if sym_type[-1] == 'h':
        printf(fout, 'R', center, matrix[:,1])
    if ((sym_type[0]=='C' and sym_type[-1]=='v') or
        (sym_type[0]=='D' and sym_type[-1]=='h')):
        # A series of reflections.
        reflect_normal = matrix[:,2]
        rotation_axis = matrix[:,1]
        for _ in xrange(degree):
            printf(fout, 'R', center, reflect_normal)
            reflect_normal = rotateAxis(reflect_normal, rotation_axis, np.pi/degree)
    if sym_type[0]=='D':
        # A series of rotations.
        minor_axis = matrix[:,2]
        main_axis = matrix[:,1]
        for _ in xrange(degree):
            printf(fout, 'C', center, minor_axis, 2)
            minor_axis = rotateAxis(minor_axis, main_axis, np.pi/degree)
    fout.close()
    return

def convertOGroup(center, matrix, output):
    """Convert sym3 annotation to a list of reflections and rotations
        for the special case of O type (cubes).

    Args:
        center: A (3,) numpy array, position of center.
        matrix: A (3,3) numpy array, the rotation matrix from sym3 file.
        output: A string. Output file name.
    """
    fout = open(output, 'w')
    for i in xrange(3):
        # Face center axes.
        printf(fout, 'R', center, matrix[:,i])
        printf(fout, 'C', center, matrix[:,i], 4)
    for i in xrange(3):
        for w in [-1.0, 1.0]:
            # Edge center axes.
            axis = matrix[:,i] + w*matrix[:,(i+1) % 3]
            axis = axis / np.linalg.norm(axis)
            printf(fout, 'R', center, axis)
            printf(fout, 'C', center, axis, 2)
    for w1 in [-1.0, 1.0]:
        for w2 in [-1.0, 1.0]:
            # Vertices axis.
            axis = matrix[:,0] + w1*matrix[:,1] + w2*matrix[:,2]
            axis = axis / np.linalg.norm(axis)
            printf(fout, 'R', center, axis)
            printf(fout, 'C', center, axis, 3)
    return

def main(synset):
    # files = glob.glob(os.path.join('./Results/', synset, '*.sym3'))
    # files = glob.glob(os.path.join('/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/', synset, '38bdba5f6455c5ff91663a74ccd2338*.sym'))
    # files = glob.glob(os.path.join('/afs/cs.stanford.edu/u/anastasiad/Lingxiao/samplefilesforsymmetrydetection/Results/*.sym'))
    files = glob.glob(os.path.join('/afs/cs.stanford.edu/u/anastasiad/Lingxiao/repptforsiemensfridaymeeting/Results/*.sym'))
    for f in files:
        print(f)
        output = f.replace('.sym', '.generators')
        convert(f, output)
    return

if __name__ == '__main__':
    synset = sys.argv[1]
    main(synset)

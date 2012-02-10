"""
cube.py

Read a cube file to an array for manipulation. Use the translational symmtry of
a supercell to average a cube file into a smaller subsection, eg the unit cell.
The number of grid points must be exactly divisible by the folding factors, FX
FY FZ.

Modified version for faps, provides Cube object.

"""

import bz2
import gzip
from glob import glob
from logging import info, debug, error

import numpy as np
from numpy import array, zeros, matrix


class Cube(object):
    """Container for a .cube file. Very specific for folding symmetry."""
    def __init__(self, filename=None, fold=None, debug=False):
        # User defined params
        self.filename = filename
        self.fold = fold
        self.debug = debug
        # Attributes for internal use
        self.error = 0.0
        self.header_block = []
        self.natoms = 0
        self.grid = [0, 0, 0]
        self.rgrid = [0, 0, 0]
        self.cell = np.eye(3)
        self.datapoints = np.array([])
        # Do we initialise here?
        if self.filename is not None:
            self.read_file()

    def read_file(self, filename=None, fold=None, crop_atoms=True):
        """Read the gridded cubedata and fold"""

        if filename is not None:
            self.filename = filename
        if fold is not None:
            self.fold = fold
        elif self.fold is None:
            self.fold = (1, 1, 1)
        fold = self.fold
        cube_temp = compressed_open(self.filename)
        top_bit = cube_temp.readlines(1024)
        cube_temp.seek(0)

        self.natoms = int(top_bit[2].split()[0])
        self.grid[0] = abs(int(top_bit[3].split()[0]))
        self.grid[1] = abs(int(top_bit[4].split()[0]))
        self.grid[2] = abs(int(top_bit[5].split()[0]))

        # Abort if we can't divide evenly!
        if self.grid[0] % fold[0]:
            raise ValueError("Grid points in x, %i, "
                             "not divisible by fold factor, %i" %
                             (self.grid[0], fold[0]))
        elif self.grid[1] % fold[1]:
            raise ValueError("Grid points in y, %i, "
                             "not divisible by fold factor, %i" %
                             (self.grid[1], fold[1]))
        elif self.grid[2] % fold[2]:
            raise ValueError("Grid points in z, %i, "
                             "not divisible by fold factor, %i" %
                             (self.grid[2], fold[2]))

        self.rgrid[0] = self.grid[0]/fold[0]
        self.rgrid[1] = self.grid[1]/fold[1]
        self.rgrid[2] = self.grid[2]/fold[2]

        # read in the top bits -- don't change for now
        for _lineno in range(self.natoms+6):
            self.header_block.append(cube_temp.readline())

        self.header_block[3:6] = [
            "%6i" % (self.rgrid[0]) + self.header_block[3][6:],
            "%6i" % (self.rgrid[1]) + self.header_block[4][6:],
            "%6i" % (self.rgrid[2]) + self.header_block[5][6:]]

        self.cell = np.array(
            [[float(x) for x in self.header_block[3].split()[1:]],
             [float(x) for x in self.header_block[4].split()[1:]],
             [float(x) for x in self.header_block[5].split()[1:]]])*0.529177249

        if crop_atoms:
            self.header_block = in_cell(self.header_block, self.rgrid, self.cell)
            self.natoms = len(self.header_block) - 6


        # Fromfile might be slower and can't deal with zipped files,
        # but uses less memory
        data_start_pos = cube_temp.tell()
        localdata = zeros([np.prod(fold)] + self.rgrid)
        stacked_data = True
        try:
            cube_data = np.fromstring(cube_temp.read(), sep=' ').reshape(self.grid)
        except MemoryError:
            try:
                cube_temp.seek(data_start_pos)
                cube_data = np.fromfile(cube_temp, sep=' ').reshape(self.grid)
            except MemoryError:
                # Read line by line directly into folded grid,
                # slowest but with least memory
                stacked_data = False
                cube_temp.seek(data_start_pos)
                localdata = zeros(self.rgrid)
                xidx, yidx, zidx = 0, 0, 0
                for line in cube_temp:
                    for point in line.split():
                        point = float(point)
                        localdata[xidx % self.rgrid[0],
                                  yidx % self.rgrid[1],
                                  zidx % self.rgrid[2]] += point
                        zidx += 1
                        if zidx == self.grid[2]:
                            zidx = 0
                            yidx += 1
                            if yidx == self.grid[1]:
                                yidx = 0
                            xidx += 1

        cube_temp.close()

        if stacked_data:
            for xidx in range(fold[0]):
                for yidx in range(fold[1]):
                    for zidx in range(fold[2]):
                        grid_idx = zidx+yidx*fold[2]+xidx*fold[2]*fold[1]
                        localdata[grid_idx] = cube_data[
                            (xidx*self.grid[0])/fold[0]:((xidx+1)*self.grid[0])/fold[0],
                            (yidx*self.grid[1])/fold[1]:((yidx+1)*self.grid[1])/fold[1],
                            (zidx*self.grid[2])/fold[2]:((zidx+1)*self.grid[2])/fold[2]]

            del cube_data
            self.datapoints = np.mean(localdata, axis=0)
            stdev = np.std(localdata, axis=0)
            self.error = ((np.sum(stdev)/np.sum(self.datapoints))/
                          np.flatnonzero(self.datapoints).size)
            info("Estimated error in cube file: %g" % self.error)
            if self.debug:
                self.write_generic(rstdev, self.error_name)
        else:
            self.datapoints = localdata/float(fold[0]*fold[1]*fold[2])

    def write_cube(self, outname=None):
        """Write out the data, already folded"""

        if outname is None:
            outname = self.folded_name

        outfile = open(outname, "w")
        outfile.writelines(self.header_block)
        for xidx in range(self.rgrid[0]):
            for yidx in range(self.rgrid[1]):
                for zidx in range(self.rgrid[2]):
                    outfile.write("%13.5E" % self.datapoints[xidx, yidx, zidx])
                    if zidx % 6 == 5:
                        outfile.write("\n")
                outfile.write("\n")

        outfile.flush()
        outfile.close()

    def write_generic(self, data, outname):
        """Write data to a Gaussian '.cube' file."""

        outfile = open(outname, "w")
        outfile.writelines(self.header_block)
        for xidx in range(self.rgrid[0]):
            for yidx in range(self.rgrid[1]):
                for zidx in range(self.rgrid[2]):
                    outfile.write("%13.5E" % data[xidx, yidx, zidx])
                    if zidx % 6 == 5:
                        outfile.write("\n")
                outfile.write("\n")

        outfile.flush()
        outfile.close()

    def maxima(self):
        """
        Use a median filter to reduce noise, smooth with gaussian blur
        then use the 124 nearest neighbours to estimate positions of maxima.
        Return the cartesian positions of maxima in a tuple with their
        magnitudes from the smoothed data.
        """
        try:
            from scipy.ndimage.filters import gaussian_filter, maximum_filter
            from scipy.ndimage.filters import median_filter
            from scipy.ndimage.morphology import generate_binary_structure
            from scipy.ndimage.morphology import binary_erosion
            from scipy.ndimage.morphology import iterate_structure
        except ImportError:
            error("Scipy not found, skipping maxima finding")
            return []

        temp_data = self.datapoints
        normalising_sum = sum(temp_data)

        # First run a median filter to remove any exceptionally low or high
        # Size seems to work best at about 0.2 A
        temp_data = median_filter(temp_data, 4, mode='wrap')
        # Gaussian filter smoothes out the data
        temp_data = gaussian_filter(temp_data, 4, mode="wrap")

        # Renormalise to pre-filtered values
        temp_data *= normalising_sum/sum(temp_data)

        # define a connectivity neighborhood
        neighborhood = generate_binary_structure(np.ndim(temp_data), 3)
        # expand it to a 5x5x5 grid
        neighborhood = iterate_structure(neighborhood, 2)

        #apply the local maximum filter; all pixel of maximal value
        #in their neighborhood are set to 1
        local_max = maximum_filter(temp_data, footprint=neighborhood, mode='wrap') == temp_data
        #local_max is a mask that contains the peaks we are
        #looking for, but also the background.
        #In order to isolate the peaks we must remove the background from the mask.

        #we create the mask of the background
        background = (temp_data == 0)

        #a little technicality: we must erode the background in order to
        #successfully subtract it form local_max, otherwise a line will
        #appear along the background border (artifact of the local maximum filter)
        eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)

        #we obtain the final mask, containing only peaks,
        #by removing the background from the local_max mask
        detected_peaks = local_max - eroded_background

        cell = np.array(
            [[float(x) for x in self.header_block[3].split()[1:]],
             [float(x) for x in self.header_block[4].split()[1:]],
             [float(x) for x in self.header_block[5].split()[1:]]])*0.529177249

        peaks = np.where(detected_peaks)
        cartesian_peaks = []
        for point in zip(peaks[0], peaks[1], peaks[2]):
            # Some random peaks with no magnitude were showing up at the PBCs
            if temp_data[point] > 0.0:
                cartesian_peaks.append((np.dot(point, cell).tolist(),
                                        temp_data[point]))

        return cartesian_peaks

    @property
    def folded_name(self):
        """File name with _folded inserted for output."""
        if ".cube" in self.filename:
            return self.filename.replace(".cube", "_folded.cube")
        else:
            return self.filename + "_folded.cube"

    @property
    def error_name(self):
        """File name with _folded inserted for output."""
        if ".cube" in self.filename:
            return self.filename.replace(".cube", "_error.cube")
        else:
            return self.filename + "_error.cube"

def in_cell(header_block, grid, cell):
    """Cut any atoms that are not in the box from the header block"""
    rcell = matrix(cell).I
    newlines = []
    for line in header_block[6:]:
        cpos = [float(x) for x in line.split()[2:]]
        fpos = array(cpos*rcell)[0]
        if (fpos[0] < grid[0]) and (fpos[1] < grid[1]) and (fpos[2] < grid[2]):
            newlines.append(line)

    header_block[2] = "%6i" % (len(newlines)) + header_block[2][6:]
    return header_block[:6] + newlines


def compressed_open(filename):
    """Return file objects for either compressed and uncompressed files"""
    filenames = glob(filename) + glob(filename+".gz") + glob(filename+".bz2")
    try:
        filename = filenames[0]
    except IndexError:
        raise IOError("File not found: %s" % filename)
    if filename[-4:] == ".bz2":
        return bz2.BZ2File(filename)
    elif filename[-3:] == ".gz":
        return gzip.open(filename)
    else:
        return open(filename, "r")

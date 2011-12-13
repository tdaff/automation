#!/usr/bin/env python

"""
Use the translational symmtry of a supercell to average a cube file into a
smaller subsection, eg the unit cell. The number of grid points must be
exactly divisible by the folding factors, FX FY FZ. Optionally, produce cubes
with data averaged along an axis. Output is written to OUTPREFIX.cube etc.

"""

# TD 2011-05
# Bugs aplenty

import numpy as np
from numpy import array, zeros, matrix
from optparse import OptionParser
import bz2
import gzip


GO = "\033[1;32m>> \033[0m\033[1m"
OK = "\033[1;32m * \033[0m"
WARN = "\033[1;33m * \033[0m"
ERR = "\033[1;31m * \033[0m\033[1m"
RESET = "\033[0m"

class Cube(object):
    """Container for a .cube file. Very specific for folding symmetry."""
    def __init__(self):
        self.filename = None
        self.header_block = []
        self.natoms = 0
        self.grid = [0, 0, 0]
        self.rgrid = [0, 0, 0]
        self.datapoints = np.array([])
        self.fold = (1, 1, 1)

    def read_file(self, filename="in.cube", fold=(1, 1, 1), crop_atoms=True,
                  verbose=False):
        """Read the gridded cubedata and fold"""

        if verbose:
            print(GO + "Opening file %s" % filename + RESET)
        self.filename = filename
        self.fold = fold
        filetemp = compressed_open(self.filename)
        top_bit = filetemp.readlines(1024)
        filetemp.seek(0)

        self.natoms = int(top_bit[2].split()[0])
        self.grid[0] = abs(int(top_bit[3].split()[0]))
        self.grid[1] = abs(int(top_bit[4].split()[0]))
        self.grid[2] = abs(int(top_bit[5].split()[0]))

        # Abort if we can't divide evenly!
        if self.grid[0] % fold[0]:
            abort("Grid points in x, %i, not divisible by fold factor, %i" %
                  (self.grid[0], fold[0]))
        elif self.grid[1] % fold[1]:
            abort("Grid points in y, %i, not divisible by fold factor, %i" %
                  (self.grid[1], fold[1]))
        elif self.grid[2] % fold[2]:
            abort("Grid points in z, %i, not divisible by fold factor, %i" %
                  (self.grid[2], fold[2]))

        self.rgrid[0] = self.grid[0]/fold[0]
        self.rgrid[1] = self.grid[1]/fold[1]
        self.rgrid[2] = self.grid[2]/fold[2]

        if verbose:
            print(OK +
                  "Folding gridpoints %s -> %s" % (self.grid, self.rgrid) +
                  RESET)

        localdata = zeros(self.rgrid)
        # read in the top bits -- don't change for now
        for _lineno in range(self.natoms+6):
            self.header_block.append(filetemp.readline())

        self.header_block[3:6] = [
            "%6i" % (self.rgrid[0]) + self.header_block[3][6:],
            "%6i" % (self.rgrid[1]) + self.header_block[4][6:],
            "%6i" % (self.rgrid[2]) + self.header_block[5][6:]]

        if crop_atoms:
            self.header_block = in_cell(self.header_block, self.rgrid)
            self.natoms = len(self.header_block) - 6
            if verbose:
                print(OK + "Cropped to %i in-cell atoms" % self.natoms + RESET)

        # read the rest of the file
        xidx, yidx, zidx = 0, 0, 0

        for line in filetemp:
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
                        if verbose:
                            print(OK + "Read: %i/%i\r" % (xidx, self.grid[0])),

        if verbose:
            print("\n" + OK + "Done reading!" + RESET)
        filetemp.close()
        self.datapoints = localdata/float(fold[0]*fold[1]*fold[2])

    def smooth_data(self, factor=1, verbose=False):
        """Apply a gaussian filter to the data."""
        if verbose:
            print(GO + "Smoothing data" + RESET)
            print(WARN + "Smooting will reduce accuracy!" + RESET)
        try:
            import scipy.ndimage
        except ImportError:
            print(WARN + "scipy not installed, skipping smoothing" + RESET)
            return 1
        self.datapoints = scipy.ndimage.gaussian_filter(
            self.datapoints, factor, mode="wrap")

    def zoom_data(self, zoom_factor, verbose=False):
        """Spline interpolate the data on a new grid."""
        if verbose:
            print(GO + "Zooming data by a factor of %f" % zoom_factor + RESET)
            print(WARN + "Zooming routines will make your data noisy\n" +
                  WARN + "and reduce accuracy!" + RESET)
        try:
            import scipy.ndimage
        except ImportError:
            print(WARN + "scipy not installed, skipping zooming" + RESET)
            return 1
        self.datapoints = scipy.ndimage.interpolation.zoom(
            self.datapoints, zoom_factor, mode='wrap')
        new_rgrid = self.datapoints.shape
        gridvector = [float(x)*(float(self.rgrid[0])/new_rgrid[0])
                      for x in self.header_block[3].split()[1:]]
        self.header_block[3] = "%6i%12.6f%12.6f%12.6f\n" % tuple([new_rgrid[0]]
                                                                  + gridvector)
        gridvector = [float(x)*(float(self.rgrid[1])/new_rgrid[1])
                      for x in self.header_block[4].split()[1:]]
        self.header_block[4] = "%6i%12.6f%12.6f%12.6f\n" % tuple([new_rgrid[1]]
                                                                  + gridvector)
        gridvector = [float(x)*(float(self.rgrid[2])/new_rgrid[2])
                      for x in self.header_block[5].split()[1:]]
        self.header_block[5] = "%6i%12.6f%12.6f%12.6f\n" % tuple([new_rgrid[2]]
                                                                  + gridvector)
        if verbose:
            print(OK + "Zoom grid %s -> %s" % (self.rgrid, new_rgrid) + RESET)
        self.rgrid = new_rgrid
        return 0


    def write_cube(self, basename="mini", verbose=False):
        """Write out the data, already folded"""

        outname = basename + ".cube"
        if verbose:
            print(GO + "Writing file %s" % outname + RESET)

#        rescale = 1.0/(self.fold[0]*self.fold[1]*self.fold[2])

        outfile = open(outname, "w")
        outfile.writelines(self.header_block)
        for xidx in range(self.rgrid[0]):
            if verbose:
                print(OK + "%i/%i\r" % (xidx, self.rgrid[0])),
            for yidx in range(self.rgrid[1]):
                for zidx in range(self.rgrid[2]):
                    outfile.write("%13.5E" % self.datapoints[xidx, yidx, zidx])
                    if zidx % 6 == 5:
                        outfile.write("\n")
                outfile.write("\n")
            #outfile.write("\n")

        outfile.flush()
        outfile.close()
        if verbose:
            print("\n" + OK + "Done writing!" + RESET)

    def flat_cube_a(self, basename="mini", verbose=False):
        """Write out the data, already folded"""

        suffix = "-a.cube"
        flatdata = np.average(self.datapoints, axis=0)

        outname = basename + suffix
        if verbose:
            print(GO + "Writing file %s" % outname + RESET)

        outfile = open(outname, "w")
        outfile.writelines(self.header_block)
        for xidx in range(self.rgrid[0]):
            if verbose:
                print(OK + "%i/%i\r" % (xidx, self.rgrid[0])),
            for yidx in range(self.rgrid[1]):
                for zidx in range(self.rgrid[2]):
                    outfile.write("%13.5E" % flatdata[yidx, zidx])
                    if zidx % 6 == 5:
                        outfile.write("\n")
                outfile.write("\n")
#            outfile.write("\n")

        outfile.flush()
        outfile.close()
        if verbose:
            print("\n" + OK + "Done writing!" + RESET)

    def flat_cube_b(self, basename="mini", verbose=False):
        """Write out the data, already folded"""

        suffix = "-b.cube"
        flatdata = np.average(self.datapoints, axis=1)

        outname = basename + suffix
        if verbose:
            print(GO + "Writing file %s" % outname + RESET)

        outfile = open(outname, "w")
        outfile.writelines(self.header_block)
        for xidx in range(self.rgrid[0]):
            if verbose:
                print(OK + "%i/%i\r" % (xidx, self.rgrid[0])),
            for _yidx in range(self.rgrid[1]):
                for zidx in range(self.rgrid[2]):
                    outfile.write("%13.5E" % flatdata[xidx, zidx])
                    if zidx % 6 == 5:
                        outfile.write("\n")
                outfile.write("\n")
#            outfile.write("\n")

        outfile.flush()
        outfile.close()
        if verbose:
            print("\n" + OK + "Done writing!" + RESET)

    def flat_cube_c(self, basename="mini", verbose=False):
        """Write out the data, already folded"""

        suffix = "-c.cube"
        flatdata = np.average(self.datapoints, axis=2)

        outname = basename + suffix
        if verbose:
            print(GO + "Writing file %s" % outname + RESET)

        outfile = open(outname, "w")
        outfile.writelines(self.header_block)
        for xidx in range(self.rgrid[0]):
            if verbose:
                print(OK + "%i/%i\r" % (xidx, self.rgrid[0])),
            for yidx in range(self.rgrid[1]):
                for zidx in range(self.rgrid[2]):
                    outfile.write("%13.5E" % flatdata[xidx, yidx])
                    if zidx % 6 == 5:
                        outfile.write("\n")
                outfile.write("\n")
#            outfile.write("\n")

        outfile.flush()
        outfile.close()
        if verbose:
            print("\n" + OK + "Done writing!" + RESET)

    def maxima(self):
        from scipy.ndimage.filters import gaussian_filter, maximum_filter, median_filter
        from scipy.ndimage.morphology import generate_binary_structure, binary_erosion

        temp_data = self.datapoints

        #temp_data = np.where(temp_data > 0.1*temp_data.max(), temp_data, 0)

        temp_data = median_filter(temp_data, 4, mode='wrap')
        temp_data = gaussian_filter(temp_data, 4, mode="wrap")
        temp_data = median_filter(temp_data, 4, mode='wrap')
        temp_data = gaussian_filter(temp_data, 4, mode="wrap")

        # define a connectivity neighborhood
        neighborhood = generate_binary_structure(np.ndim(temp_data), 3)

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
        max_file = open("maxima.xyz", "wb")
        max_file.write(" %i\nMaxima\n" % len(peaks[0]))
        print(" %i\nMaxima" % len(peaks[0]))
        for point in zip(peaks[0], peaks[1], peaks[2]):
            print("Xx    %12.8f %12.8f %12.8f" % tuple(np.dot(point, cell)))
            max_file.write("Xx    %12.8f %12.8f %12.8f\n" % tuple(np.dot(point, cell)))


def in_cell(header_block, grid):
    """Cut any atoms that are not in the box from the header block"""
    rcell = matrix([[float(x) for x in header_block[3].split()[1:]],
                    [float(x) for x in header_block[4].split()[1:]],
                    [float(x) for x in header_block[5].split()[1:]]]).I
    newlines = []
    for line in header_block[6:]:
        cpos = [float(x) for x in line.split()[2:]]
        fpos = array(cpos*rcell)[0]
        if (fpos[0] < grid[0]) and (fpos[1] < grid[1]) and (fpos[2] < grid[2]):
            newlines.append(line)

    header_block[2] = "%6i" % (len(newlines)) + header_block[2][6:]
    return header_block[:6] + newlines


def get_opts():
    """Define options, some error checking"""
    usage = "usage: %prog [options] INFILE FX FY FZ OUTPREFIX"
    parser = OptionParser(usage=usage, version="%prog 0.1",
                          description=__doc__)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="write extra information to stdout [default]",
                      default=True)
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                      help="silence all output")
    parser.add_option("-c", "--crop", action="store_true", dest="crop",
                      help="crop atoms to those in the cell [default]",
                      default=True)
    parser.add_option("-n", "--no-crop", action="store_false", dest="crop",
                      help="do not crop atoms; use original configuration")
    parser.add_option("-z", "--zoom", dest="zoom", type='float',
                      help="resize the grid with a spline interpolation")
    parser.add_option("-s", "--smooth", action="store_true", dest="smooth",
                      help="Smooth data with a gaussian filter")
    parser.add_option("-A", action="store_true", dest="flatten_a",
                      help="flatten along the a axis")
    parser.add_option("-B", action="store_true", dest="flatten_b",
                      help="flatten along the b axis")
    parser.add_option("-C", action="store_true", dest="flatten_c",
                      help="flatten along the c axis")
    parser.add_option("-m", "--maxima", action="store_true", dest="maxima",
                      help="search for local maxima in the density")
    (local_options, local_args) = parser.parse_args()

    if len(local_args) != 5:
        parser.error("Five arguments are required (try %prog --help)")

    return (local_options, local_args)

def abort(message):
    """Deliver the message and exit"""
    print(ERR + message + RESET)
    raise SystemExit(1)


def compressed_open(filename):
    """Return file objects for either compressed and uncompressed files"""
    if filename[-4:] == ".bz2":
        return bz2.BZ2File(filename)
    elif filename[-3:] == ".gz":
        return gzip.open(filename)
    else:
        return open(filename, "r")


def main():
    """Main program logic"""
    options, args = get_opts()
    mydata = Cube()
    mydata.read_file(filename=args[0], fold=tuple(int(x) for x in args[1:4]),
                     verbose=options.verbose, crop_atoms=options.crop)
    # Smooth before anything else if we must
    if options.smooth:
        mydata.smooth_data(verbose=options.verbose)
    # Rescale the data if we need to
    if options.zoom:
        mydata.zoom_data(options.zoom, verbose=options.verbose)
    if options.maxima:
        mydata.maxima()
        raise SystemExit
    # Write the folded cube
    mydata.write_cube(basename=args[4], verbose=options.verbose)
    # Write any flattened data requested
    if options.flatten_a:
        mydata.flat_cube_a(basename=args[4], verbose=options.verbose)
    if options.flatten_b:
        mydata.flat_cube_b(basename=args[4], verbose=options.verbose)
    if options.flatten_c:
        mydata.flat_cube_c(basename=args[4], verbose=options.verbose)


if __name__ == "__main__":
    main()

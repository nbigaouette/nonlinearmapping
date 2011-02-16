#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# ***************************************************************************
# Parse command line arguments
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-v", "--verbose",            action="store_true",    dest="verbose",     default=False,      help="Show all output. [default: %default]")
parser.add_option("-m", "--mapping",type=str,                           dest="mapping_type",default="potential",help="Mapping type. [default: %default]")
parser.add_option("-e", "--even",               action="store_true",    dest="equal",       default=False,      help="Spread the atoms evenly in domain (overwrites -x). [default: %default]")
parser.add_option("-x",             type=float, action="append",        dest="x0s",         default=None,       help="Ion positions [default: [3,5,8] Bohr]")
parser.add_option("-d", "--dxmin",  type=float,                         dest="dxmin",       default=0.05,       help="Minimum cell size. [default: %default Bohr]")
parser.add_option("--xmin",         type=float,                         dest="xmin",        default=0.0,        help="Domain lower bound. [default: %default Bohr]")
parser.add_option("--xmax",         type=float,                         dest="xmax",        default=10.0,       help="Domain upper bound. [default: %default Bohr]")
parser.add_option("-n",             type=int,                           dest="ni",          default=50,         help="Total number of grid points. [default: %default Bohr]")

(options, args) = parser.parse_args()
# ***************************************************************************

# Default values
if (options.x0s == None):
    options.x0s = [3.0, 5.0, 8.0]

# Imported necessary python modules
import numpy, sys, math
import matplotlib.pyplot as plt

# Import specific figure parameters
import matplotlib_params

# ***************************************************************************
class NonLinearMapping:
    """
        Generic class to build the mapping.
    """

    def Initialize(self, imin, imax, xmin, xmax, x0, dxmin):
        """
            Initialize generic mapping.
            Arguments:
                imin:   Minimum value for "i"
                imax:   Maximum value for "i"
                xmin:   Minimum value for "x"
                xmax:   Maximum value for "x"
                x0:     Location of center of interest
                dxmin:  Minimum cell size
        """

        self.imin  = imin
        self.imax  = imax
        self.xmin  = xmin
        self.xmax  = xmax
        self.x0    = x0
        self.dxmin = dxmin

        self.n     = imax - imin
        self.d     = 0.0
        self.A     = 0.0
        self.i_x0md= 0.0
        self.i_x0pd= 0.0
    #

    def Allocate_and_Set_is(self):
        """
            Allocate memory for arrays (continuous and discrete) in the three different regions (1, 2 and 3).
        """

        # Discrete
        if (self.i_x0md > self.imin):
            self.i1 = numpy.arange(int(math.ceil(self.imin)),       int(math.floor(self.i_x0md))+1)
        else:
            self.i1 = numpy.asarray([])
        self.i2     = numpy.arange(int(math.ceil(self.i_x0md)),     int(math.floor(self.i_x0pd))+1)
        if (self.i_x0pd < self.imax):
            self.i3     = numpy.arange(int(math.ceil(self.i_x0pd)), int(math.floor(self.imax))+1)  # TODO: Fix im-1
        else:
            self.i3 = numpy.asarray([])

        # Continuous
        self.nii    = 1e5
        if (self.i_x0md > 0.0):
            self.ii1= numpy.arange(int(math.ceil(self.imin)),       self.i_x0md,    (self.i_x0md - self.imin)   / self.nii)
        else:
            self.ii1= numpy.asarray([])
        self.ii2    = numpy.arange(self.i_x0md,                     self.i_x0pd,    (self.i_x0pd - self.i_x0md) / self.nii)
        if (self.i_x0pd < self.imax):
            self.ii3= numpy.linspace(self.i_x0pd,                   self.imax,      self.nii)
        else:
            self.ii3= numpy.asarray([])
    #

    def Asserts(self):
        """
            Verify all calculated values to prevent breakage.
        """

        if (self.x0 + self.d > self.xmax):
            raise ValueError( \
                "\n  Linear region of ion escape the domain!" + \
                "\n  x0 + d = "+str(self.x0)+" + "+str(self.d)+" = "+str(self.x0+self.d)+" >= xmax = "+str(self.xmax) + \
                "\n  Try a smaller cell size.")
        if (self.x0 - self.d < self.xmin):
            raise ValueError( \
                "\n  Linear region of ion escape the domain!" + \
                "\n  x0 - d = "+str(self.x0)+" - "+str(self.d)+" = "+str(self.x0-self.d)+" <= xmin = "+str(self.xmin) + \
                "\n  Try a smaller cell size.")
        if (self.i_x0md < self.imin):
            raise ValueError( \
                "\n  Linear region of ion escape the domain!" + \
                "\n  i_x0md = "+str(self.i_x0md)+" <= imin = "+str(self.imin) + \
                "\n  Try a smaller cell size.")
        if (self.i_x0pd > self.imax):
            raise ValueError( \
                "\n  Linear region of ion escape the domain!" + \
                "\n  i_x0pd = "+str(self.i_x0pd)+" <= imax = "+str(self.imax) + \
                "\n  Try a smaller cell size.")

        assert(self.A > 0.0)

        # Verify that x(i(x)) == x and i(x(i)) == i
        tmp_x = self.xmin*0.99999
        assert(abs(tmp_x - self.Calculate_x(self.Calculate_i(tmp_x))) < 1.0e-10)

        tmp_x = self.xmax*0.99999
        assert(abs(tmp_x - self.Calculate_x(self.Calculate_i(tmp_x))) < 1.0e-10)

        tmp_i = self.imin*0.99999
        assert(abs(tmp_i - self.Calculate_i(self.Calculate_x(tmp_i))) < 1.0e-10)

        tmp_i = self.imax*0.99999
        assert(abs(tmp_i - self.Calculate_i(self.Calculate_x(tmp_i))) < 1.0e-10)

        # Verify i_x0md and i_x0pd
        assert(self.i_x0md >= self.imin)
        assert(self.i_x0pd >= self.imin)
        assert(self.i_x0md <= self.imax)
        assert(self.i_x0pd <= self.imax)
    #

    def Print(self):
        """
            Print generic mapping parameters.
        """

        print "i     = [" + str(self.imin) + ", " + str(self.imax) + "["
        print "x     = [" + str(self.xmin) + ", " + str(self.xmax) + "["
        print "dxmin =", self.dxmin
        print "d     =", self.d
        print "i_x0md=", self.i_x0md
        print "i_x0pd=", self.i_x0pd
        print "A     =", self.A
        print "x0 - d=", self.x0-self.d
        print "x0    =", self.x0
        print "x0 + d=", self.x0+self.d
        print "Calculate_x(imin) =", self.Calculate_x(self.imin)
        try:
            print "Calculate_x(i_x0md) =", self.Calculate_x(self.i_x0md)
        except:
            pass
        try:
            print "Calculate_x(i_x0pd) =", self.Calculate_x(self.i_x0pd)
        except:
            pass
        print "Calculate_x(imax) =", self.Calculate_x(self.imax)
        try:
            print "i1    =", self.i1
            print "i2    =", self.i2
            print "i3    =", self.i3
            print "x1    =", self.x1
            print "x2    =", self.x2
            print "x3    =", self.x3
        except:
            pass
    #

    def Calculate_i_x0md_i_x0pd(self):
        """ Calculate i(x0 - d) and i(x0 + d) """

        assert(self.d >= 1.0e-8)
        assert(self.dxmin >= 1.0e-8)
        self.i_x0md = self.Calculate_i1(self.x0-self.d)
        self.i_x0pd = self.i_x0md + 2.0*self.d/self.dxmin
    #

    def Calculate_i2(self, x):
        """ Calculate i2(x) (linear region 2). See equation (21b) """
        return (x - (self.x0 - self.d)) / self.dxmin + self.imin
    def Calculate_x2(self, i):
        """ Calculate x2(x) (linear region 2) See equation (23b) """
        return self.dxmin * (i - self.i_x0md) + self.xmin
    def Calculate_dx2di(self, i):
        """ Calculate del x2/del i (linear region 2) """
        return self.dxmin * numpy.ones(len(i))
    def Calculate_d2x2di2(self, i):
        """ Calculate del^2 x2/del i^2 (linear region 2) """
        return numpy.zeros((len(i)))
    #

    def Calculate_Mapping(self):
        """ Populate the (previously) allocated arrays with the mapped values. See equations (19) """

        # Make sure the linear region is inside the domain
        if (self.x0 + self.d > self.xmax):
            raise ValueError("ERROR: Linear region of ion escape the domain!")
        if (self.x0 - self.d < self.xmin):
            raise ValueError("ERROR: Linear region of ion escape the domain!")

        # Discrete
        self.x1  = self.Calculate_x1(self.i1)
        self.x2  = self.Calculate_x1(self.i_x0md) + self.Calculate_x2(self.i2) - self.xmin
        self.x3  = self.Calculate_x1(self.i_x0md) + self.Calculate_x2(self.i_x0pd) + self.Calculate_x3(self.i3) - 2.0*self.xmin

        self.dx1 = self.Calculate_dx1di(self.i1)
        self.dx2 = self.Calculate_dx2di(self.i2)
        self.dx3 = self.Calculate_dx3di(self.i3)

        self.ddx1 = self.Calculate_d2x1di2(self.i1)
        self.ddx2 = self.Calculate_d2x2di2(self.i2)
        self.ddx3 = self.Calculate_d2x3di2(self.i3)

        # Continuous
        self.xx1 = self.Calculate_x1(self.ii1)
        self.xx2 = self.Calculate_x1(self.i_x0md) + self.Calculate_x2(self.ii2) - self.xmin
        self.xx3 = self.Calculate_x1(self.i_x0md) + self.Calculate_x2(self.i_x0pd) + self.Calculate_x3(self.ii3) - 2.0*self.xmin

        self.dxx1 = self.Calculate_dx1di(self.ii1)
        self.dxx2 = self.Calculate_dx2di(self.ii2)
        self.dxx3 = self.Calculate_dx3di(self.ii3)

        self.ddxx1 = self.Calculate_d2x1di2(self.ii1)
        self.ddxx2 = self.Calculate_d2x2di2(self.ii2)
        self.ddxx3 = self.Calculate_d2x3di2(self.ii3)
    #

    def Calculate_i(self, x):
        """ Calculate i(x) by testing in which region x is. See equation (16). """

        if (x < (self.xmin - 1.0e-2)):
            raise ValueError("Calculate_i(x) called with a value x = "+ str(x) +" lower then xmin = " + str(self.xmin))
        elif (x < self.x0-self.d):
            return self.Calculate_i1(x)
        elif (x < self.x0+self.d):
            return self.Calculate_i1(self.x0-self.d) + self.Calculate_i2(x) - self.imin
        elif (x <= self.xmax*1.0001):
            return self.Calculate_i1(self.x0-self.d) + self.Calculate_i2(self.x0+self.d) + self.Calculate_i3(x) - 2.0*self.imin
        else:
            raise ValueError("Calculate_i(x) called with a value x = "+ str(x) +" greater then xmax = " + str(self.xmax))
        #
    #

    def Calculate_x(self, i):
        """ Calculate x(i) by testing in which region i is. """

        if (i < (self.imin - 1.0e-2)):
            raise ValueError("Calculate_x(i) called with a value i = "+ str(i) +" lower then imin = " + str(self.imin))
        elif (i < self.i_x0md):
            return self.Calculate_x1(i)
        elif (i < self.i_x0pd):
            return self.Calculate_x1(self.i_x0md) + self.Calculate_x2(i) - self.xmin
        elif (i < self.imax*1.0001):
            return self.Calculate_x1(self.i_x0md) + self.Calculate_x2(self.i_x0pd) + self.Calculate_x3(i) - 2.0*self.xmin
        else:
            raise ValueError("Calculate_x(i) called with a value i = "+ str(i) +" greater then imax = " + str(self.imax))
        #
    #

    def Get_i(self):
        """ Return i(x) for the whole subdomain (discrete values). """
        return numpy.concatenate((self.i1, self.i2, self.i3))
    def Get_x(self):
        """ Return x(i) for the whole subdomain (discrete values). """
        return numpy.concatenate((self.x1, self.x2, self.x3))
    def Get_dx(self):
        """ Return del x / del i for the whole subdomain (discrete values). """
        return numpy.concatenate((self.dx1, self.dx2, self.dx3))
    def Get_ddx(self):
        """ Return del^2 x / del i^2 for the whole subdomain (discrete values). """
        return numpy.concatenate((self.ddx1, self.ddx2, self.ddx3))
    def Get_ii(self):
        """ Return i(x) for the whole subdomain (continuous values). """
        return numpy.concatenate((self.ii1, self.ii2, self.ii3))
    def Get_xx(self):
        """ Return x(i) for the whole subdomain (continuous values). """
        return numpy.concatenate((self.xx1, self.xx2, self.xx3))
    def Get_dxx(self):
        """ Return del x / del i for the whole subdomain (continuous values). """
        return numpy.concatenate((self.dxx1, self.dxx2, self.dxx3))
    def Get_ddxx(self):
        """ Return del^2 x / del i^2 for the whole subdomain (continuous values). """
        return numpy.concatenate((self.ddxx1, self.ddxx2, self.ddxx3))
    #
#

# ***************************************************************************
class SqrtMapping(NonLinearMapping):
    """
        Mapping with square root as source function S(x'). See equation (15).
    """

    def Initialize(self, imin, imax, xmin, xmax, x0, dxmin):
        """
            Initialize mapping.
            Arguments:
                imin:   Minimum value for "i"
                imax:   Maximum value for "i"
                xmin:   Minimum value for "x"
                xmax:   Maximum value for "x"
                x0:     Location of center of interest
                dxmin:  Minimum cell size
        """

        # Initialize generic mapping
        NonLinearMapping.Initialize(self, imin, imax, xmin, xmax, x0, dxmin)

        # Needed parameters for square root mapping.
        if (2.0 * self.x0 > self.xmax):
            discriminant = (numpy.sqrt(self.x0) + numpy.sqrt(self.xmax - self.x0))**2 - 2.0*self.imax*self.dxmin
            self.u1 = ((numpy.sqrt(self.x0) + numpy.sqrt(self.xmax - self.x0)) + numpy.sqrt(discriminant)) / 2.0
            self.u2 = ((numpy.sqrt(self.x0) + numpy.sqrt(self.xmax - self.x0)) - numpy.sqrt(discriminant)) / 2.0
            self.d1 = self.u1**2
            self.d2 = self.u2**2
            self.d = self.d2
            self.A = (self.imax - 2.0*self.d/self.dxmin) / (numpy.sqrt(self.x0) - 2.0*numpy.sqrt(self.d) + numpy.sqrt(self.xmax - self.x0))
        elif (2.0 * self.x0 < self.xmax):
            discriminant = (numpy.sqrt(self.xmax - self.x0) + numpy.sqrt(self.x0))**2 - 2.0*self.imax*self.dxmin
            self.u1 = ( (numpy.sqrt(self.xmax - self.x0) + numpy.sqrt(self.x0)) + numpy.sqrt(discriminant) ) / 2.0
            self.u2 = ( (numpy.sqrt(self.xmax - self.x0) + numpy.sqrt(self.x0)) - numpy.sqrt(discriminant) ) / 2.0
            self.d1 = self.u1**2
            self.d2 = self.u2**2
            self.d = self.d2
            self.A = (self.imax - 2.0*self.d/self.dxmin) / ( numpy.sqrt(self.xmax - self.x0) + numpy.sqrt(self.x0) - 2.0*numpy.sqrt(self.d) )
        else:
            discriminant = self.x0 - self.imax*self.dxmin/2.0
            self.u1 = numpy.sqrt(self.x0) + numpy.sqrt(discriminant)
            self.u2 = numpy.sqrt(self.x0) - numpy.sqrt(discriminant)
            self.d1 = self.u1**2
            self.d2 = self.u2**2
            self.d = self.d2
            self.A = (self.imax - 2.0*self.d/self.dxmin) / (2.0 * (numpy.sqrt(self.x0) - numpy.sqrt(self.d)))
        #
        self.B = self.A * numpy.sqrt(self.x0)
        self.D = -numpy.sqrt(self.d)

        self.Calculate_i_x0md_i_x0pd()

        NonLinearMapping.Allocate_and_Set_is(self)
        self.Asserts()

        NonLinearMapping.Calculate_Mapping(self)
        print "ddx1 =", self.ddx1
        print "ddx2 =", self.ddx2
        print "ddx3 =", self.ddx3
        print "ddx1.shape =", self.ddx1.shape," ddx2 =", self.ddx2.shape," ddx3 =", self.ddx3.shape
    #

    def Print(self):
        """
            Print specific mapping parameters.
        """

        print "u1    =", self.u1
        print "u2    =", self.u2
        print "d1    =", self.d1
        print "d2    =", self.d2
        print "B     =", self.B
        print "D     =", self.D

        # Print generic parameters
        NonLinearMapping.Print(self)

    def Asserts(self):
        """
            Verify specific parameters.
        """

        assert(not numpy.isnan(self.u1))
        assert(not numpy.isnan(self.u2))
        assert(not numpy.isnan(self.d1))
        assert(not numpy.isnan(self.d2))
        assert(not numpy.isnan(self.d))
        assert(not numpy.isnan(self.A))

        # Verify generic parameters
        NonLinearMapping.Asserts(self)

    def Calculate_i1(self, x):
        """ Calculate i1(x) (region 1). """
        return -self.A * numpy.sqrt(self.x0 - x) + self.B
    def Calculate_i3(self, x):
        """ Calculate i3(x) (region 3). """
        return self.A * (numpy.sqrt(x - self.x0) + self.D)
    #

    def Calculate_x1(self, i):
        """ Calculate x1(i) (region 1). """
        return self.x0 - (numpy.sqrt(self.x0) - i/self.A)**2
    def Calculate_x3(self, i):
        """ Calculate x3(i) (region 3). """
        return ((i-self.i_x0pd)/self.A - self.D)**2 - self.d
    #

    def Calculate_dx1di(self, i):
        """ Calculate del x1 / del i (region 1). """
        return 2.0/self.A * (numpy.sqrt(self.x0) - i/self.A)
    def Calculate_dx3di(self, i):
        """ Calculate del x3 / del i (region 3). """
        return 2.0/self.A * ( (i-self.i_x0pd)/self.A - self.D )
    #

    def Calculate_d2x1di2(self, i):
        """ Calculate del^2 x1 / del i^2 (region 1). """
        return -2.0/self.A**2 * numpy.ones(len(i))
    def Calculate_d2x3di2(self, i):
        """ Calculate del^2 x3 / del i^2 (region 3). """
        return 2.0/self.A**2 * numpy.ones(len(i))
    #
# class SqrtMapping(NonLinearMapping):

# ***************************************************************************
class PotentialMapping(NonLinearMapping):
    """
        Mapping with electrostatic potential as source function S(x'). See equation (15).
    """

    def fct_of_d(self, d):
        """ Calculate equation (22) as a function of "d" """
        # a is box's lower bound
        # b is box's upper bound
        if (d == 0.0):
            # limit(f, d-> 0) == imax * dxmin
            return (self.imax - self.imin)*self.dxmin
        else:
            return (self.imax - self.imin)*self.dxmin - d*( 2.0 + numpy.log( (self.x0-self.xmin)*(self.xmax-self.x0) / d**2 ) )
        #
    #

    def diff_fct_of_d(self, d):
        """ Calculate derivative of equation (22) with respect to "d" """
        # a is box's lower bound
        # b is box's upper bound
        if (d == 0.0):
            # Infinity, but don't overflow
            # limit(diff(f,d), d-> 0) = -infity
            return -1.0e100
        else:
            return -numpy.log( (self.x0-self.xmin)*(self.xmax-self.x0) / d**2 )
    #

    def d_bisection(self):
        """ Solve equation (22) for "d" using the bisection method. The initial guess is d = 0. """

        # We know the shape of fct_of_d(): For positive "d", it looks like
        # a "U". If the lowest part of the function is negative, there is a root,
        # else, there isn't. To start the bisection, we find the lowest point
        # of the function by finding where the derivative is zero. The
        # derivative diff_fct_of_d() is zero at d = sqrt( (x0-xmin)*(xmax-x0) )
        d_left = 0.0
        inside_sqrt = (self.x0-self.xmin)*(self.xmax-self.x0)
        d_right = math.sqrt(inside_sqrt)
        assert(-1.0e-6 < self.diff_fct_of_d(d_right))
        assert(self.diff_fct_of_d(d_right) < 1.0e-6)

        # But then, to have a root, fct_of_d(d_right) MUST be negative!
        if (self.fct_of_d(d_right) > 0.0):
            print "*********************************************************"
            print "ERROR: There is no root to be found with these parameters!"
            print "x0 =", self.x0
            print "dxmin =", self.dxmin
            print "imin =", self.imin
            print "imax =", self.imax
            print "xmin =", self.xmin
            print "xmax =", self.xmax
            print "Location of the minimum of the function:"
            print "Function's minimum's location:", d_right
            print "Value of the function there:", self.fct_of_d(d_right), "but we are expecting 0.0"
            sys.exit(0)

        d = d_right + (d_left - d_right) / 2.0
        while (d != d_left and d != d_right):
            #print "d =", d, " d_left =", d_left, " d_right =", d_right, " f(d) =", self.fct_of_d(d), " diff(f,d) =", self.diff_fct_of_d(d)
            if (self.fct_of_d(d) <= 0.0):
                d_right = d
            else:
                d_left = d
            d = d_right + (d_left - d_right) / 2.0
        #

        #print "fct_of_d(d="+str(d)+") =", self.fct_of_d(d)

        assert(abs(self.fct_of_d(d)) < 1.0e-4)
        assert(d > 0.0)

        return d
    #

    def Initialize(self, imin, imax, xmin, xmax, x0, dxmin):
        """
            Initialize mapping.
            Arguments:
                imin:   Minimum value for "i"
                imax:   Maximum value for "i"
                xmin:   Minimum value for "x"
                xmax:   Maximum value for "x"
                x0:     Location of center of interest
                dxmin:  Minimum cell size
        """

        # Initialize generic mapping
        NonLinearMapping.Initialize(self, imin, imax, xmin, xmax, x0, dxmin)

        # Needed parameters for electrostatic potential mapping.
        self.d = self.d_bisection()
        self.A = self.d / self.dxmin

        self.Calculate_i_x0md_i_x0pd()

        NonLinearMapping.Allocate_and_Set_is(self)
        self.Asserts()
        NonLinearMapping.Calculate_Mapping(self)
    #

    def Asserts(self):
        """
            Verify specific parameters.
        """

        j = self.A * numpy.log((self.x0-self.xmin)/self.d) + self.imin
        assert(abs(self.i_x0md - j) < 1.0e-10)

        # Verify generic parameters
        NonLinearMapping.Asserts(self)
    #

    def Print(self):
        """
            Print specific mapping parameters.
        """
        # Print generic parameters
        NonLinearMapping.Print(self)
    #

    def Calculate_i1(self, x):
        """ Calculate i1(x) (region 1). See equation (21a). """
        return self.A * numpy.log((self.x0 - self.xmin) / (self.x0 - x)) + self.imin
    def Calculate_i3(self, x):
        """ Calculate i3(x) (region 3). See equation (21c). """
        return self.A * numpy.log((x - self.x0) / self.d) + self.imin
    #

    def Calculate_x1(self, i):
        """ Calculate x1(i) (region 1). See equation (23a). """
        return self.x0 + (self.xmin - self.x0)*numpy.exp((self.imin - i)/self.A)
    def Calculate_x3(self, i):
        """ Calculate x3(i) (region 3). See equation (23c). """
        return self.d * (numpy.exp((i-self.i_x0pd) / self.A) - 1.0) + self.xmin
    #

    def Calculate_dx1di(self, i):
        """ Calculate del x1 / del i (region 1). """
        return -(self.xmin - self.x0)*numpy.exp((self.imin - i) / self.A) / self.A
    def Calculate_dx3di(self, i):
        """ Calculate del x3 / del i (region 3). """
        return self.d/self.A * numpy.exp((i-self.i_x0pd)/self.A)
    #

    def Calculate_d2x1di2(self, i):
        """ Calculate del^2 x1 / del i^2 (region 1). """
        return (self.xmin - self.x0)*numpy.exp(-i / self.A) / self.A**2
    def Calculate_d2x3di2(self, i):
        """ Calculate del^2 x3 / del i^2 (region 3). """
        return self.d/self.A**2 * numpy.exp((i-self.i_x0pd)/self.A)
    #
# class PotentialMapping(NonLinearMapping):

# ***************************************************************************
class FieldMapping(NonLinearMapping):
    """
        Mapping with electrostatic field as source function S(x'). See equation (15).
    """

    def Initialize(self, imin, imax, xmin, xmax, x0, dxmin):
        """
            Initialize mapping.
            Arguments:
                imin:   Minimum value for "i"
                imax:   Maximum value for "i"
                xmin:   Minimum value for "x"
                xmax:   Maximum value for "x"
                x0:     Location of center of interest
                dxmin:  Minimum cell size
        """

        # Initialize generic mapping
        NonLinearMapping.Initialize(self, imin, imax, xmin, xmax, x0, dxmin)

        # Needed parameters for electrostatic field mapping.
        tmp = 1.0/(self.xmax - self.x0) + 1.0/(self.x0 - self.xmin)
        discriminant = 4.0 - tmp*self.imax*self.dxmin
        assert(discriminant > 0.0)
        self.d1 = (2.0 + numpy.sqrt(discriminant)) / tmp
        self.d2 = (2.0 - numpy.sqrt(discriminant)) / tmp
        self.d = self.d2
        self.A = self.d**2 / self.dxmin

        NonLinearMapping.Calculate_i_x0md_i_x0pd(self)
        NonLinearMapping.Allocate_and_Set_is(self)
        self.Asserts()
        NonLinearMapping.Calculate_Mapping(self)
    #

    def Asserts(self):
        """
            Verify specific parameters.
        """

        # Verify generic parameters
        NonLinearMapping.Asserts(self)
    #

    def Print(self):
        """
            Print specific mapping parameters.
        """

        print "d1    =", self.d1
        print "d2    =", self.d2

        # Print generic parameters
        NonLinearMapping.Print(self)

    def Calculate_i1(self, x):
        """ Calculate i1(x) (region 1). """
        return self.A * ( 1.0/(self.x0 - x) - 1.0/(self.x0 - self.xmin) )
    def Calculate_i3(self, x):
        """ Calculate i3(x) (region 3). """
        return self.A * ( 1.0/(self.x0 - x) + 1.0/self.d )
    #

    def Calculate_x1(self, i):
        """ Calculate x1(i) (region 1). """
        return self.x0 - self.A / (i + self.A/(self.x0-self.xmin))
    def Calculate_x3(self, i):
        """ Calculate x3(i) (region 3). """
        return - 1.0 / ( (i-self.i_x0pd)/self.A - 1.0/self.d ) - self.d
    #

    def Calculate_dx1di(self, i):
        """ Calculate del x1 / del i (region 1). """
        return self.A / (self.A/(self.x0 - self.xmin) + i)**2
    def Calculate_dx3di(self, i):
        """ Calculate del x3 / del i (region 3). """
        return 1.0 / (self.A * ( (i-self.i_x0pd)/self.A - 1.0/self.d )**2)
    #

    def Calculate_d2x1di2(self, i):
        """ Calculate del^2 x1 / del i^2 (region 1). """
        return -2.0*self.A / (self.A/(self.x0 - self.xmin) + i)**3
    def Calculate_d2x3di2(self, i):
        """ Calculate del^2 x3 / del i^2 (region 3). """
        return -2.0 / (self.A**2 * ( (i-self.i_x0pd)/self.A - 1.0/self.d )**3)
    #
#class FieldMapping(NonLinearMapping):

# ***************************************************************************
class LinearMapping(NonLinearMapping):
    """
        Linear mapping class. Linear mapping is achieved by
        making the linear region 2 encompass the whole (sub)domain.
    """

    def Initialize(self, imin, imax, xmin, xmax, x0, dxmin):
        """
            Initialize mapping.
            Arguments:
                imin:   Minimum value for "i"
                imax:   Maximum value for "i"
                xmin:   Minimum value for "x"
                xmax:   Maximum value for "x"
                x0:     Location of center of interest
                dxmin:  Minimum cell size
        """

        # Initialize generic mapping
        NonLinearMapping.Initialize(self, imin, imax, xmin, xmax, x0, dxmin)

        # Needed parameters for linear mapping.
        self.A  = 1.0e100 # Just to shut up Asserts(), not used anywhere
        self.x0 = (self.xmax + self.xmin) / 2.0
        self.d  = self.xmax - self.x0
        self.dxmin = (self.xmax - self.xmin) / (self.imax - self.imin)

        self.i_x0md = self.imin
        self.i_x0pd = self.imax

        NonLinearMapping.Allocate_and_Set_is(self)
        self.Asserts()
        NonLinearMapping.Calculate_Mapping(self)
    #

    def Asserts(self):
        """
            Verify specific parameters.
        """

        # Verify generic parameters
        NonLinearMapping.Asserts(self)
    #

    def Print(self):
        """
            Print specific mapping parameters.
        """

        # Print generic parameters
        NonLinearMapping.Print(self)

    def Calculate_i1(self, x):
        """ Calculate i1(x) (region 1). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(x)))
        except:
            return 0.0
    def Calculate_i3(self, x):
        """ Calculate i3(x) (region 3). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(x)))
        except:
            return 0.0
    #

    def Calculate_x1(self, i):
        """ Calculate x1(i) (region 1). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(i)))
        except:
            return 0.0
    def Calculate_x3(self, i):
        """ Calculate x3(i) (region 3). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(i)))
        except:
            return 0.0
    #

    def Calculate_dx1di(self, i):
        """ Calculate del x1 / del i (region 1). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(i)))
        except:
            return 0.0
    def Calculate_dx3di(self, i):
        """ Calculate del x3 / del i (region 3). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(i)))
        except:
            return 0.0
    #

    def Calculate_d2x1di2(self, i):
        """ Calculate del^2 x1 / del i^2 (region 1). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(i)))
        except:
            return 0.0
    def Calculate_d2x3di2(self, i):
        """ Calculate del^2 x3 / del i^2 (region 3). This should not be called since in linear mapping region 1 is non-existent. """
        try:
            return numpy.zeros((len(i)))
        except:
            return 0.0
    #
#class LinearMapping(NonLinearMapping):

class Subdomain:
    """ Class to ease up subdomain boundaries """
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 0.0
        self.imin = 0.0
        self.imax = 0.0
    def Set(self, imin, imax, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax
        self.imin = imin
        self.imax = imax
    def Size_i(self):
        return (self.imax - self.imin)
    def Size_x(self):
        return (self.xmax - self.xmin)
    #
#


# ***************************************************************************
def mapping_nonlinear(xmin, xmax, ni, dxmin = 0.1, x0s = None,
                      i_x0md = None, i_x0pd = None, ds = None):
    """
        Calculate mapping of a domain using given parameters.
        Input arguments:
            xmin:   Minimum value of domain's x positions
            xmax:   Maximum value of domain's x positions
            ni:     Total number of cells
            dxmin:  Minimum cell size
            x0s:    List containing ions location
        Output arguments:
            i_x0md: List containing i(x0 - d) for each ion. Useful for plotting.
            i_x0pd: List containing i(x0 + d) for each ion. Useful for plotting.
            ds: List containing "d" for each ion. Useful for plotting.
    """

    if (x0s == None):
        nb_ions = 1
        x0s = [xmax / 2.0]
    else:
        nb_ions = len(x0s)
    #

    domain_size = xmax - xmin

    xstart = xmin
    xstop  = xmin
    istart = 0.0
    istop  = 0.0

    ii  = numpy.zeros((0))
    xx  = numpy.zeros((0))
    dxx = numpy.zeros((0))
    ddxx= numpy.zeros((0))
    i   = numpy.zeros((0))
    x   = numpy.zeros((0))
    dx  = numpy.zeros((0))
    ddx = numpy.zeros((0))

    imax_prev_ion = 0.0
    xmax_prev_ion = 0.0

    if (options.verbose):
        print "************************************"
        print "Domain: "
        print "domain_size =", domain_size
        print "x: [", xmin, ",", xmax, "["
        print "ni =", ni

    # If a subdomain is to receive less then 10% of cells, get some cells from
    # its neigbour.
    # Make sure each subdomain receive at least 10% of the cells
    subdomains = [Subdomain() for n in xrange(nb_ions)]
    prev_xstop = 0.0
    for n in xrange(nb_ions):
        # Set spacial limits
        subdomains[n].xmin = prev_xstop
        if (n+1 == nb_ions):
            # If last ion, take domain boundary
            subdomains[n].xmax = xmax
        else:
            # If not last ion, take middle point with the next ion
            subdomains[n].xmax = (x0s[n] + x0s[n+1]) / 2.0
        #

        prev_xstop = subdomains[n].xmax

    # Calculate ratio to total domain size for each subdomains
    subdomains_ratios = numpy.zeros((nb_ions), dtype=float)
    for n in xrange(nb_ions):
        subdomains_ratios[n] = subdomains[n].Size_x() / domain_size
        #print "subdomains_ratios["+str(n)+"] = ", subdomains_ratios[n]

    # Make sure each rations is at least 10%. If not, take 1% from its neigbours
    stop = False
    wanted_ratio = 0.10
    while (not stop):
        for n in xrange(nb_ions):
            if (subdomains_ratios[n] < wanted_ratio):
                ratio_missing = wanted_ratio - subdomains_ratios[n]
                nb_sd_to_remove = 0
                if (n != 0):
                    nb_sd_to_remove += 1
                if (n != nb_ions-1):
                    nb_sd_to_remove += 1
                if (n != 0):
                    subdomains_ratios[n-1] -= ratio_missing/nb_sd_to_remove
                if (n != nb_ions-1):
                    subdomains_ratios[n+1] -= ratio_missing/nb_sd_to_remove
                subdomains_ratios[n] += ratio_missing
                break
            if (n == nb_ions-1):
                stop = True
    assert(subdomains_ratios.sum() - 1.0 < 1.0e-5)

    for n in xrange(nb_ions):
        if (options.verbose):
            print "************************************"
            print "Ion #" + str(n), " position =", x0s[n]

        # Set x limits
        xstart = xstop
        if (n+1 == nb_ions):
            # If last ion, take domain boundary
            xstop = xmax
        else:
            # If not last ion, take middle point with the next ion
            xstop = (x0s[n] + x0s[n+1]) / 2.0
        #

        subdomain_size = xstop - xstart

        # Give each subdomain their fair share of the domain
        delta_i = ni * subdomains_ratios[n]

        # If there is no subdomain to the left, remove one half
        if (n == 0):
            delta_i -= 0.5
        # If there is no subdomain to the right, remove one half
        if (n+1 == nb_ions):
            delta_i -= 0.5

        # Start at previous' stop
        istart = istop
        istop += delta_i
        if (n+1 == nb_ions):
            # Make sure the last subdomain falls on the number of cells minus 1
            assert(abs(istop - (ni-1.0)) < 1.0e-5)
        #

        x0 = x0s[n]

        if (options.verbose):
            print "delta_i          =", delta_i
            print "subdomain_size   =", subdomain_size
        assert(subdomain_size > 0.0)


        # Choose which mapping to plot
        if (options.mapping_type == "linear"):
            mapping_obj = LinearMapping()
        elif (options.mapping_type == "field"):
            mapping_obj = FieldMapping()
        elif (options.mapping_type == "sqrt"):
            mapping_obj = SqrtMapping()
        elif (options.mapping_type == "potential"):
            mapping_obj = PotentialMapping()

        mapping_obj.Initialize(istart, istop, xstart, xstop, x0, dxmin)
        if (options.verbose):
            mapping_obj.Print()

        # Save values for later use.
        if (i_x0md != None):
            i_x0md.append(mapping_obj.i_x0md)
        if (i_x0pd != None):
            i_x0pd.append(mapping_obj.i_x0pd)
        if (ds != None):
            ds.append(mapping_obj.d)

        i   = numpy.concatenate((i,     mapping_obj.Get_i()+imax_prev_ion))
        x   = numpy.concatenate((x,     mapping_obj.Get_x()+xmax_prev_ion))
        dx  = numpy.concatenate((dx,    mapping_obj.Get_dx()))
        ddx = numpy.concatenate((ddx,   mapping_obj.Get_ddx()))

        ii  = numpy.concatenate((ii,    mapping_obj.Get_ii()+imax_prev_ion))
        xx  = numpy.concatenate((xx,    mapping_obj.Get_xx()+xmax_prev_ion))
        dxx = numpy.concatenate((dxx,   mapping_obj.Get_dxx()))
        ddxx= numpy.concatenate((ddxx,  mapping_obj.Get_ddxx()))
    #

    return i, x, dx, ddx, ii, xx, dxx, ddxx
# def mapping_nonlinear()


# ***************************************************************************
def main():
    show_figure_mapping_xi  = True
    show_figure_all_mapping = False

    nb_ions = len(options.x0s)

    distance = 2.0 # Distance between each ions [bohr]

    # Set ions' locations
    if (options.equal):
        x0s = numpy.zeros((nb_ions), dtype=float)
        Zs  = numpy.ones((nb_ions), dtype=float)
        xmiddle = (xmin + xmax) / 2.0
        xstart = xmiddle - (nb_ions * distance / 2.0)
        xcm = 0.0
        for n in xrange(nb_ions):
            x0s[n] = (float(n) * distance)
            xcm += x0s[n]
        xcm /= nb_ions
        for n in xrange(nb_ions):
            x0s[n] += -xcm + xmiddle
    else:
        x0s = options.x0s

    ds = []
    i_x0md = []
    i_x0pd = []

    i, x, J1, J2, ii, xx, dxx, ddxx = mapping_nonlinear(options.xmin, options.xmax, options.ni, options.dxmin, x0s, i_x0md, i_x0pd, ds)

    print "######################################################################################################"
    print "Domain range:   [" + str(options.xmin) + ", " + str(options.xmax) + "] Bohr"
    print "Ions positions:", x0s, "Bohr"
    print "Total number of cells:", int(options.ni)
    print "Subdomains' width 'd':", ds, "Bohr"

    # ***************************************************************************
    # Plot the mapping x(i)
    if (show_figure_mapping_xi):
        fig = plt.figure()
        axprops = dict()

        ax1 = plt.subplot(111)
        plt.plot(ii, xx, label="Continuous")
        plt.plot(i, x, "xr", label="Discrete (integers)")
        plt.xlabel(r"$i$")
        plt.ylabel(r"$x$ (Bohr)")

        assert(nb_ions == len(ds))
        assert(nb_ions == len(i_x0md))
        assert(nb_ions == len(i_x0pd))

        hl = 3.0 * ii[-1] / 32.0
        vl = xx[-1] / 8.0
        # Horizontal lines
        for j in xrange(nb_ions):
            plt.plot([i_x0md[j]-hl, i_x0pd[j]+hl], [x0s[j]-ds[j], x0s[j]-ds[j]], ':k')
            plt.plot([i_x0md[j]-hl, i_x0pd[j]+hl], [x0s[j]+ds[j], x0s[j]+ds[j]], ':k')
        # Vertical lines
        for j in xrange(nb_ions):
            plt.plot([i_x0md[j], i_x0md[j]], [x0s[j]-ds[j]-vl, x0s[j]+ds[j]+vl], ':k')
            plt.plot([i_x0pd[j], i_x0pd[j]], [x0s[j]-ds[j]-vl, x0s[j]+ds[j]+vl], ':k')

        # Add arrows surrounding linear regions
        arrow_length = (options.xmax-options.xmin)/100.0*7.0 # 7% of figure's vertical axis range
        head_length  = arrow_length/2.0
        gap = head_length / 2.0
        #alpha = 0.75
        alpha = 1.0
        for j in xrange(nb_ions):
            #arrow_x = (ni/16.0) * (1.0 + (ai%2))
            arrow_x = i_x0md[j] - (15.0*hl/16.0)
            plt.arrow(arrow_x, x0s[j]-ds[j]-arrow_length-head_length-gap, 0.0,  arrow_length, color = 'k', head_length=head_length, linewidth=3.0, head_width=1.0, alpha = alpha)
            plt.arrow(arrow_x, x0s[j]+ds[j]+arrow_length+head_length+gap, 0.0, -arrow_length, color = 'k', head_length=head_length, linewidth=3.0, head_width=1.0, alpha = alpha)
            plt.text( arrow_x, x0s[j], r'$2d_' + str(j) + '$.', horizontalalignment='right', verticalalignment='center')

        # Add ions' positions to ylabels
        yaxis, yaxis_label = plt.yticks()
        yaxis_label = [0]*len(yaxis)
        for ai in xrange(len(yaxis)):
            #print "ai = ", ai
            if (str('%.0f' % yaxis[ai]) == "8"):
                yaxis_label[ai] = r''
                continue
            yaxis_label[ai] = r'$'+str('%.0f' % yaxis[ai])+'$'
        for x0i in xrange(len(x0s)):
            #yaxis_label.append(r'$x_{\rm{ion}_' + str(x0i) + r'}$')
            yaxis_label.append(r'$x_{0,' + str(x0i) + r'}$')
            yaxis       = numpy.append(yaxis,       x0s[x0i])
        plt.yticks(yaxis, yaxis_label)

        ax1.set_ylim((options.xmin, options.xmax))
        ax1.set_xlim((0.0, options.ni-1.0))

        # By explicitly setting the yaxis labels, matplotlib will fail to detect
        # the mouse's vertical position (what's reported in the lower right corner of the window).
        # So clone the axis (and hide it), set the right limits so a mouse over
        # will correctly report the position.
        old_yaxis = plt.twinx()
        plt.setp(old_yaxis.get_yticklabels(), visible=False)
        old_yaxis.set_ylim((options.xmin, options.xmax))

        xaxis, xaxis_label = plt.xticks()
        xaxis[-1] = int(ii[-1])
        plt.xticks(xaxis)


    # ***************************************************************************
    if (show_figure_all_mapping):
        fig = plt.figure()
        axprops = dict()

        im_x0       = 0.125
        im_y0       = 0.125
        im_width    = 0.85
        #im_height   = 0.266667
        im_height   = 1.0/3.0 - 1.0/15.0
        im_gap      = 0.0

        ax1 = fig.add_axes([im_x0, im_y0+2*(im_height+im_gap), im_width, im_height], **axprops)
        axprops['sharex'] = ax1
        plt.plot(ii, xx, label=r'$x(i)$ (continuous)')
        plt.plot(i, x, 'xr', label=r'$x(i)$ (discrete)')
        for n in xrange(nb_ions):
            plt.plot([0.0, ii.max()], [x0s[n]-ds[n], x0s[n]-ds[n]], ':k')
            plt.plot([0.0, ii.max()], [x0s[n],       x0s[n]],       ':k')
            plt.plot([0.0, ii.max()], [x0s[n]+ds[n], x0s[n]+ds[n]], ':k')
        for j in xrange(len(i_x0md)):
            plt.plot([i_x0md[j], i_x0md[j]], [options.xmin, options.xmax], ':k')
            plt.plot([i_x0pd[j], i_x0pd[j]], [options.xmin, options.xmax], ':k')
        ax1.set_ylim((options.xmin, options.xmax))
        plt.grid()
        plt.ylabel(r"$x$")
        plt.legend()

        ax2 = fig.add_axes([im_x0, im_y0+im_height+im_gap, im_width, im_height], **axprops)
        plt.plot(ii, dxx, '-b', label=r'$J_1(i)$ (continuous)')
        plt.plot(i, J1, 'xr', label=r'$J_1(i)$ (discrete)')
        dii = numpy.concatenate(([ii[1]-ii[0]],ii[1:-1]-ii[0:-2],[ii[-2]-ii[-1]]))
        plt.plot(ii, numpy.gradient(xx, dii), ':m', label=r"$\frac{\Delta x(i)}{\Delta i}$")
        if (dxx.max() != 0.0):
            ax2.set_yscale('log')
        plt.grid()
        plt.ylabel(r"$J_1$")
        plt.legend()

        ax3 = fig.add_axes([im_x0, im_y0, im_width, im_height], **axprops)
        plt.plot(ii, ddxx)
        plt.plot(i, J2, 'xr')
        plt.plot(ii, numpy.gradient(dxx, dii), ':m')
        plt.xlabel(r"$i$")
        plt.ylabel(r"$J_2$")
        plt.grid()
        ax3.set_ylim((-5.0, 5.0))

        for ax in ax1, ax2:
            plt.setp(ax.get_xticklabels(), visible=False)
        for ax in ax1, ax2, ax3:
            ax.set_xlim((0.0, options.ni-1.0))

    plt.show()
# def main()


if __name__ == "__main__":
    main()
#


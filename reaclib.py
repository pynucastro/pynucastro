# parse the reaclib stuff

import numpy as np

class Tfactors(object):
    """ precompute temperature factors for speed """

    def __init__(self, T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)


class SingleSet(object):
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

        A single rate in Reaclib can be composed of multiple sets
    """

    def __init__(self, a, label=None):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.label = label


    def f(self):
        """
        return a function for this set -- note: Tf here is a Tfactors
        object
        """
        return lambda tf: np.exp(self.a[0] + 
                                 self.a[1]*tf.T9i +
                                 self.a[2]*tf.T913i +
                                 self.a[3]*tf.T913 +
                                 self.a[4]*tf.T9 +
                                 self.a[5]*tf.T953 +
                                 self.a[6]*tf.lnT9)


    def set_string(self, prefix="set", plus_equal=False):
        """ 
        return a string containing the python code for this set
        """
        if plus_equal:
            string =  "{} += np.exp( ".format(prefix)
        else:
            string =  "{} = np.exp( ".format(prefix)
        string += " {}".format(self.a[0])
        if not self.a[1] == 0.0: string += " + {}*tf.T9i".format(self.a[1])
        if not self.a[2] == 0.0: string += " + {}*tf.T913i".format(self.a[2])
        if not self.a[3] == 0.0: string += " + {}*tf.T913".format(self.a[3])
        string += "\n{}         ".format(len(prefix)*" ")
        if not self.a[4] == 0.0: string += " + {}*tf.T9".format(self.a[4])
        if not self.a[5] == 0.0: string += " + {}*tf.T953".format(self.a[5])
        if not self.a[6] == 0.0: string += " + {}*tf.lnT9".format(self.a[6])
        string += ")\n"
        return string


class Rate(object):
    """ a single Reaclib rate, which can be composed of multiple sets """

    def __init__(self, file):
        self.file = file
        self.chapter = None    # the Reaclib chapter for this reaction
        self.reactants = []
        self.products = []
        self.sets = []

        self.dens_exp = 1
        self.prefactor = 1.0    # this is 1/2 for rates like a + a (double counting)

        self.Q = 0.0

        # read in the file, parse the different sets and store them as
        # SingleSet objects in sets[]
        f = open(file, "r")
        lines = f.readlines()

        # first line is the chapter
        self.chapter = int(lines[0])

        # remove any black lines
        set_lines = [l for l in lines[1:] if not l.strip() == ""]

        # the rest is the sets
        first = 1
        while len(set_lines) > 0:
            # sets are 3 lines long
            s1 = set_lines.pop(0)
            s2 = set_lines.pop(0)
            s3 = set_lines.pop(0)

            # first line of a set has up to 6 nuclei, then the label,
            # and finally the Q value
            f = s1.split()
            Q = f.pop()
            label = f.pop()

            if first:
                self.Q = Q

                # what's left are the nuclei -- their interpretation
                # depends on the chapter
                if self.chapter == 1:
                    # e1 -> e2
                    self.reactants.append(f[0])
                    self.products.append(f[1])

                    self.string = "{} -> {}".format(*(self.reactants + self.products))

                elif self.chapter == 2:
                    # e1 -> e2 + e3
                    self.reactants.append(f[0])
                    self.products.append(f[1])
                    self.products.append(f[2])

                    self.string = "{} -> {} + {}".format(*(self.reactants + self.products))
                    
                elif self.chapter == 3:
                    # e1 -> e2 + e3 + e4
                    self.reactants.append(f[0])
                    self.products.append(f[1])
                    self.products.append(f[2])
                    self.products.append(f[3])

                    self.string = "{} -> {} + {} + {}".format(*(self.reactants + self.products))
                    
                elif self.chapter == 4:
                    # e1 + e2 -> e3
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.products.append(f[2])
                    
                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 2

                    self.string = "{} + {} -> {}".format(*(self.reactants + self.products))

                elif self.chapter == 5:
                    # e1 + e2 -> e3 + e4
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.products.append(f[2])
                    self.products.append(f[3])
                    
                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 2

                    self.string = "{} + {} -> {} + {}".format(*(self.reactants + self.products))

                elif self.chapter == 6:
                    # e1 + e2 -> e3 + e4 + e5
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.products.append(f[2])
                    self.products.append(f[3])
                    self.products.append(f[4])
                    
                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 2

                    self.string = "{} + {} -> {} + {} + {}".format(*(self.reactants + self.products))

                elif self.chapter == 7:
                    # e1 + e2 -> e3 + e4 + e5 + e6
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.products.append(f[2])
                    self.products.append(f[3])
                    self.products.append(f[4])
                    self.products.append(f[5])

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 2

                    self.string = "{} + {} -> {} + {} + {} + {}".format(*(self.reactants + self.products))

                elif self.chapter == 8:
                    # e1 + e2 + e3 -> e4
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.reactants.append(f[2])
                    self.products.append(f[3])

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./6.  # 1/3!
                    elif len(set(self.reactants)) == 2:
                        self.prefactor = 1./2.

                    self.dens_exp = 3

                    self.string = "{} + {} + {} -> {}".format(*(self.reactants + self.products))

                elif self.chapter == 9:
                    # e1 + e2 + e3 -> e4 + e5
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.reactants.append(f[2])
                    self.products.append(f[3])
                    self.products.append(f[4])
                    
                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./6.  # 1/3!
                    elif len(set(self.reactants)) == 2:
                        self.prefactor = 1./2.

                    self.dens_exp = 3

                    self.string = "{} + {} + {} -> {} + {}".format(*(self.reactants + self.products))

                elif self.chapter == 10:
                    # e1 + e2 + e3 + e4 -> e5 + e6
                    self.reactants.append(f[0])
                    self.reactants.append(f[1])
                    self.reactants.append(f[2])
                    self.reactants.append(f[3])
                    self.products.append(f[4])
                    self.products.append(f[5])
                    
                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./24.  # 1/4!
                    elif len(set(self.reactants)) == 2:
                        # there may be some instances where we have a + a + b + b, 
                        # so prefactor = 1/4?
                        self.prefactor = 1./6. # 1/3!   
                    elif len(set(self.reactants)) == 3:
                        self.prefactor = 1./2.

                    self.dens_exp = 4

                    self.string = "{} + {} + {} + {} -> {} + {}".format(*(self.reactants + self.products))

                elif self.chapter == 11:
                    # e1 -> e2 + e3 + e4 + e5 
                    self.reactants.append(f[0])
                    self.products.append(f[1])
                    self.products.append(f[2])
                    self.products.append(f[3])
                    self.products.append(f[4])

                    self.string = "{} -> {} + {} + {} + {}".format(*(self.reactants + self.products))

                first = 0
                
            # the second line contains the first 4 coefficients
            # the third lines contains the final 3
            # we can't just use split() here, since the fields run into one another
            n = 13  # length of the field
            a = [s2[i:i+n] for i in range(0, len(s2), n)]
            a += [s3[i:i+n] for i in range(0, len(s3), n)]
            
            a = [float(e) for e in a if not e.strip() == ""]
            self.sets.append(SingleSet(a, label=label))

            
    def eval(self, T):
        """ evauate the reaction rate for temperature T """
        tf = Tfactors(T)
        r = 0.0
        for s in self.sets:
            f = s.f()
            r += f(tf)

        return r


    def rate_string(self, indent=0, prefix="rate"):
        tstring = "# {}\n".format(self.string)
        tstring += "{} = 0.0\n\n".format(prefix)

        for s in self.sets:
            tstring += "# {}\n".format(s.label)
            tstring += "{}\n".format(s.set_string(prefix=prefix, plus_equal=True))


        string = ""
        for t in tstring.split("\n"):
            string += indent*" " + t + "\n"

        return string

    
    def function_string(self):
        idx = self.file.rfind("-")
        fname = self.file[:idx].replace("--","-").replace("-","_")

        string = ""
        string += "def {}(tf):\n".format(fname)
        string += "{}".format(self.rate_string(indent=4))
        string += "    return rate\n\n"

        return string


if __name__ == "__main__":
    r = Rate("c13-pg-n14-nacr")
    print r.rate_string(indent=3)
    print r.eval(1.0e9)


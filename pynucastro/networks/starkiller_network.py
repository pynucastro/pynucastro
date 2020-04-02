"""A Fortran reaction network for integration into the StarKiller
Microphysics set of reaction networks used by astrophysical hydrodynamics
codes"""

from __future__ import print_function

import glob
import os

from pynucastro.networks import BaseFortranNetwork

class StarKillerNetwork(BaseFortranNetwork):
    def __init__(self, *args, **kwargs):
        # Initialize BaseFortranNetwork parent class
        super(StarKillerNetwork, self).__init__(*args, **kwargs)

        # StarKiller-specific template processing functions
        self.ftags['<sparse_jac_nnz>'] = self._sparse_jac_nnz
        self.ftags['<csr_jac_metadata>'] = self._csr_jac_metadata
        self.ftags['<species_xin_test>'] = self._species_xin_test

    def _get_template_files(self):
        template_pattern = os.path.join(self.pynucastro_dir,
                                        'templates',
                                        'starkiller-microphysics',
                                        '*.template')
        return glob.glob(template_pattern)

    def get_sparse_jac_nnz(self):
        # Get the number of nonzero entries in the sparse Jacobian

        # Evaluate the species Jacobian if not already done
        if not self.solved_jacobian:
            self.compose_jacobian()

        nspec = len(self.unique_nuclei)
        number_nonzero = 0

        # Count the number of nonzero entries in the species Jacobian
        jac_idx = 0
        # Loop over rows
        for j, nj in enumerate(self.unique_nuclei):
            # Loop over columns
            for i, ni in enumerate(self.unique_nuclei):
                if (not self.jac_null_entries[jac_idx]) or i==j:
                    number_nonzero += 1
                jac_idx += 1

        # Add the df/dT contributions
        number_nonzero += nspec + 2

        # Add the dTdot/dX and dedot/dX contributions
        number_nonzero += 2*nspec

        # Add the dedot/de contribution
        # This is zero but we need it because we need the
        # entire diagonal of J for the integration
        number_nonzero += 1
        return number_nonzero

    def get_csr_jac_metadata(self):
        # Get row count and column index for the CSR Jacobian
        # Remember, we need the entire diagonal for the linear solver

        # Evaluate the species Jacobian if not already done
        if not self.solved_jacobian:
            self.compose_jacobian()

        row_count = []
        col_index = []

        nspec = len(self.unique_nuclei)
        
        # Start row_count at base 1
        row_count.append(1)

        jac_idx = 0
        # Loop over rows
        for j, nj in enumerate(self.unique_nuclei):
            num_in_row = row_count[-1]
            # Loop over columns
            for i, ni in enumerate(self.unique_nuclei):
                if (not self.jac_null_entries[jac_idx]) or i==j:
                    num_in_row += 1
                    col_index.append(i+1)
                jac_idx += 1
            # Done with species in row, do T
            num_in_row += 1
            col_index.append(nspec+1)
            row_count.append(num_in_row)

        # Done with species rows, do Tdot
        num_in_row = nspec + 1 + row_count[-1]
        row_count.append(num_in_row)
        for i in range(1, nspec+2):
            col_index.append(i)

        # Done with Tdot, do edot row
        num_in_row = nspec + 2 + row_count[-1]
        row_count.append(num_in_row)
        for i in range(1, nspec+3):
            col_index.append(i)

        # Return row count and column index lists
        return row_count, col_index

    def _initial_mass_fractions(self, n_indent, of):
        # Redefine initial mass fractions tag to set the
        # mass fractions in the burn_cell unit test inputs file.
        for i, n in enumerate(self.unique_nuclei):
            of.write("\n! {nuc: <5} initial mass fraction\n".format(nuc=str(n)))
            of.write("{}massfractions({}) = 0.0d0\n".format(
                self.indent*n_indent, i+1))

    def _sparse_jac_nnz(self, n_indent, of):
        of.write('{}integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = {}\n'.format(
            self.indent*n_indent, self.get_sparse_jac_nnz()))

    def _csr_jac_metadata(self, n_indent, of):
        row_count, col_index = self.get_csr_jac_metadata()

        of.write('{}csr_jac_col_index = [ &\n'.format(
            self.indent*n_indent))
        for ci in col_index[:-1]:
            of.write('{}{}, &\n'.format(
                self.indent*(n_indent+1), ci))
        of.write('{}{}  ]\n'.format(
            self.indent*(n_indent+1), col_index[-1]))

        of.write('\n')

        of.write('{}csr_jac_row_count = [ &\n'.format(
            self.indent*n_indent))
        for ri in row_count[:-1]:
            of.write('{}{}, &\n'.format(
                self.indent*(n_indent+1), ri))
        of.write('{}{}  ]\n'.format(
            self.indent*(n_indent+1), row_count[-1]))

    def _species_xin_test(self, size_test, of):
        xcomp = 1.0/float(len(self.unique_nuclei))
        for i, n in enumerate(self.unique_nuclei):
            if i!=0:
                of.write('#\n')
            of.write('# {}\n'.format(n))
            xin = [self.fmt_to_dp_f90(xcomp) for j in range(size_test)]
            of.write('{}\n'.format(' '.join(xin)))

    def _write_network(self, use_cse=False):
        """
        This writes the RHS, jacobian and ancillary files for the system of ODEs that
        this network describes, using the template files.
        """

        super()._write_network(use_cse=use_cse)

        # create a .net file with the nuclei properties
        with open("pynucastro.net", "w") as of:
            for nuc in self.unique_nuclei:
                of.write("{:25} {:6} {:6.1f} {:6.1f}\n".format(
                    nuc.spec_name, nuc.short_spec_name, nuc.A, nuc.Z))

        # write out some network properties
        with open("NETWORK_PROPERTIES", "w") as of:
            of.write("NSCREEN {}".format(self.num_screen_calls))



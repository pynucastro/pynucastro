"""A simple C++ reaction network for integrating into other C++ codes"""


from pathlib import Path

from pynucastro.networks.base_cxx_network import BaseCxxNetwork
from pynucastro.screening import get_screening_map


class SimpleCxxNetwork(BaseCxxNetwork):
    """A basic C++ network."""

    def __init__(self, *args, **kwargs):

        # Initialize BaseCxxNetwork parent class
        super().__init__(*args, **kwargs)

        self.function_specifier = "inline"
        self.dtype = "Real"
        self.array_namespace = ""

    def _get_template_files(self):

        path = self.pynucastro_dir/"templates/simple-cxx-network"

        return path.glob("*.template")

    def _compute_screening_factors(self, n_indent, of):
        if not self.do_screening:
            screening_map = []
        else:
            screening_map = get_screening_map(self.get_rates())
        for i, scr in enumerate(screening_map):

            nuc1_info = f'{float(scr.n1.Z)}_rt, {float(scr.n1.A)}_rt'
            nuc2_info = f'{float(scr.n2.Z)}_rt, {float(scr.n2.A)}_rt'

            if not (scr.n1.dummy or scr.n2.dummy):
                # Scope the screening calculation to avoid multiple definitions of scn_fac.
                of.write(f'\n{self.indent*n_indent}' + '{\n')

                of.write(f'{self.indent*(n_indent+1)}auto scn_fac = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n')

                of.write(f'{self.indent*(n_indent+1)}actual_screen(pstate, scn_fac, scor);\n')

                of.write(f'{self.indent*n_indent}' + '}\n')

            if scr.name == "He4_He4_He4":
                # we don't need to do anything here, but we want to avoid immediately applying the screening
                pass

            elif scr.name == "He4_He4_He4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "He4_He4_He4"
                # handle the second part of the screening for 3-alpha
                of.write(f'\n{self.indent*n_indent}' + '{\n')

                of.write(f'{self.indent*(n_indent+1)}auto scn_fac2 = scrn::calculate_screen_factor({nuc1_info}, {nuc2_info});\n')

                of.write(f'{self.indent*(n_indent+1)}actual_screen(pstate, scn_fac2, scor2);\n')

                of.write(f'{self.indent*n_indent}' + '}\n')

                # there should only be a single forward 3-alpha rate
                assert len(scr.rates) == 1

                of.write('\n')
                rr = scr.rates[0]
                of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.cname()}) *= scor * scor2;\n')

            else:
                # there might be several rates that have the same
                # reactants and therefore the same screening applies
                # -- handle them all now

                of.write('\n')
                for rr in scr.rates:
                    of.write(f'{self.indent*n_indent}rate_eval.screened_rates(k_{rr.cname()}) *= scor;\n')

            of.write('\n')

    def _write_network(self, odir=None):
        """Output the the RHS, jacobian and ancillary files for the
        system of ODEs that this network describes, using the template
        files.

        """

        # at the moment, we don't support TabularRates
        assert len(self.tabular_rates) == 0, "SimpleCxxNetwork does not support tabular rates"

        super()._write_network(odir=odir)

        if odir is None:
            odir = Path.cwd()
        # create a header file with the nuclei properties
        with open(Path(odir, "network_properties.H"), "w") as of:
            of.write("#ifndef NETWORK_PROPERTIES_H\n")
            of.write("#define NETWORK_PROPERTIES_H\n")
            of.write("#include <vector>\n")
            of.write("#include <string>\n")
            of.write("#include <amrex_bridge.H>\n\n")

            of.write(f"constexpr int NumSpec = {len(self.unique_nuclei)};\n\n")

            of.write("constexpr Real aion[NumSpec] = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    {nuc.A:6.1f}, // {n}\n")
            of.write(" };\n\n")

            of.write("constexpr Real aion_inv[NumSpec] = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    1.0/{nuc.A:6.1f}, // {n}\n")
            of.write(" };\n\n")

            of.write("constexpr Real zion[NumSpec] = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    {nuc.Z:6.1f}, // {n}\n")
            of.write(" };\n\n")

            of.write("static const std::vector<std::string> spec_names = {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                of.write(f"    \"{nuc.short_spec_name.capitalize()}\", // {n}\n")
            of.write(" };\n\n")

            of.write("namespace Species {\n")
            of.write("  enum NetworkSpecies {\n")
            for n, nuc in enumerate(self.unique_nuclei):
                if n == 0:
                    of.write(f"    {nuc.short_spec_name.capitalize()}=1,\n")
                else:
                    of.write(f"    {nuc.short_spec_name.capitalize()},\n")
            of.write("  };\n")
            of.write("}\n\n")

            of.write("#endif\n")

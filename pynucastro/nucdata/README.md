# Masses and Spins

We get the mass excess and spins from the Nubase 2020 evaluation:
https://www-nds.iaea.org/amdc/

Kondev, F. G., Wang, M., Huang, W. J., Naimi, S., & Audi, G.
Chinese Physics C, 45(3), 030001. (2021) doi:10.1088/1674-1137/abddae

See the scripts in `AtomicMassEvaluation/`

We compute the binding energy from the mass excesses.

# Partition Functions

The script `convert_rathpf_2000.py` converts the files `part_etfsiq.asc.txt`
and `part_frdm.asc.txt` to a convenient format to an standard pynucastro
format. The table `part_frdm.asc.txt` is depicted in:

T. Rauscher and F.-K.Thielemann, ADNDT 75 (2000) 1.

and computed by using the finite range droplet model (FRDM model). This model
is described in detail here:

P. Möller, J. R. Nix, W. D. Myers, and W. J. Swiateki, Atomic Data and Nuclear
Data Tables, 59, 185 (1995), doi:10.1006/adnd.1995.1002

The table `part_etfsiq.asc.txt` was not mentioned in the previous article, but
the authors included their computed partition function values by using the
Thomas-Fermi approach with Strutinski integral (ETFSI-Q), described here:

J. M. Pearson, R. C. Nayak, and S. Goriely, Phys. Lett. B 387, 455 (1996).

The script `convert_rathpf_2003.py` converts the `datafile2.txt` and `datafile3.txt`
partition function table files into a pynucastro standard. The files datafile2.txt
and datafile3.txt are the computed partition values from the FRDM and ETFSI-Q
models for higher temperatures, respectively. The details of these
calculations are exposed here:

T. Rauscher, Ap. J. Suppl. 147 (2003) 403

All the previous scripts and data files are located in the `PartitionFunction`
directory.

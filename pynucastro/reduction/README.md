# CS267 Final Project
Repository for collaborating on nuclear reaction network reduction project for CS267.

The *pynucastro* library has been added as a submodule, but needs to be *pip* installed
or something similar to be available for import. Would recommend running `pip install -e
pynucastro` from this base directory in case *pynucastro* itself needs to be updated.

Let's start with the smaller network (endpoint ni56), and work up to the larger one (te108).

## Checklist

We will almost certainly not get to all of this, but we should probably at least try to
implement the two graph theory algorithms in a data parallel manner.

- [ ] Randomly sample conditions from our last set of simulation runs to create a dataset
- [ ] Determine how we are going to map the compositions (mass fractions) to our larger
      networks, and implement the mapping
- Path flux analysis:
    - Serial version:
        - [ ] Localized search algorithm (BFS?) to find length 1 and 2 paths between nuclei
        - [ ] Interaction coefficient calculation
        - [ ] Removal of unnecessary nuclei and edges from the graph
    - Parallel version:
        - [ ] Parallelize over data so that each process handles a separate set of conditions
        - [ ] Add communication code to take the maximum interaction coefficient across all processes
              for each edge and send it back to all processes before doing the reduction
- Sensitivity analysis:
    - Serial version:
        - [ ] Implement modified Dijkstra's algorithm as part of DRGEP method
        - [ ] Calculate overall interaction coefficients to identify "limbo species"
        - [ ] Define loss/error function(s) (we want to minimize the error in nuclear energy generation rate,
              change in electron fraction, and change in average atomic weight) as functions of ydot
              values
        - [ ] Loop through limbo nuclei and evaluate loss/error function(s) upon the removal of each, and
              get rid of those for which the induced error is below some threshold. We may want to test
              their removal in batches so that the overall error remains below some threshold
    - Parallel version:
        - [ ] Parallelize over data and take the maximum of the induced error across all processes for each
              limbo nucleus or set of limbo nuclei
        - [ ] Try parallelizing over the limbo nuclei as well
        - [ ] *Lower Priority:* Implement a parallel Dijkstra's for the DRGEP
- Machine learning:
    - Serial version:
        - [ ] Use backpropagation to adjust the weights on the remaining edges to minimize the loss function,
              which should just be a sum of the errors in the quantities we want to conserve the most. PyTorch
              and TensorFlow are both options here.
        - [ ] Interleave these machine learning steps with sensitivity analysis to permit more steps of the latter
- Approximate (equilibrium) rates:
    - [ ] Allow *pynucastro* to replace entire paths through the network with single rates. We would want to use
          either assumption of equilibrium or the assumption that the slowest rate is the bottleneck to get a
          reasonable starting point for the rate coefficients. We can apply the machine learning step above
          to refine the estimate.

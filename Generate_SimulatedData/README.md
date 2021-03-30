# Simulated data generation

This simulated data generation tool has been implemented as a software prototype. It provides induced-error profile, groud-truth datasets and  error-injected datasets.

Aim to generate datasets with a close nature to wet-lab miRNA sequencing reads, we have two considerations in the process. One is to computationally replicate lab-verified miRNA sequences as templates to form the basic sequences of the simulated datasets, then we duplicate these basic sequences such that the copy counts of them follow a real distribution from a wet-lab dataset of miRNA sequencing reads. In fact, we replicated the mature miRNA sequences in miRBase as the templates, and made the copy count distribution of these template sequences to follow the distribution drawn from a typical miRNA dataset.

For users' convenience, users can use their own sequence templates and copy number distribution as well.


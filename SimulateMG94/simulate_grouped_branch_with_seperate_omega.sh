#!/bin/bash


hyphy SimulateMG94.bf \
--tree CD2.nwk \
--omega-Primates 1.5 \
--omega-Others 0.25 \
--sites 400 \
--replicates 1 \
--output data/example3/example3 \
--branch-variation partitioned \
--base-frequencies HIV


#!/bin/bash

hyphy SimulateMG94.bf \
--tree CD2.nwk \
--omega 0.25 \
--sites 25 \
--replicates 5 \
--output data/example2/example2 \
--frequency-estimator F3x4 \
--base-frequencies 0.3,0.15,0.23,0.2,0.31,0.19,0.4,0.05,0.11 \
--CT 2 \
--AC 0.1 

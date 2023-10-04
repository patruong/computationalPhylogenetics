#!/bin/bash

hyphy SimulateMG94.bf \
--tree CD2.dnds.nwk \
--sites 400 \
--replicates 10 \
--output data/example4/example4 \
--branch-variation dsdn \
--base-frequencies HIV

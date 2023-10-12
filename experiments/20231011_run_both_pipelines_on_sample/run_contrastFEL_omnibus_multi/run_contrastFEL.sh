#!/bin/bash

hyphy contrast-fel --alignment sims.0.settings.replicate.1 --tree sims.0.nwk --branch-set GROUP_0 --branch-set GROUP_1 --branch-set GROUP_2 --branch-set GROUP_3 --p-value 1.00 --q-value 1.00  > contrastFEL_log.txt


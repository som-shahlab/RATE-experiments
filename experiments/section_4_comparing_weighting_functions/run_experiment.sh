#!/bin/bash

# To reproduce Figure 2 in the main manuscript
python qini_vs_autoc.py --scoring_types Oracle --n_sims 10000

# To reproduce Figure 3 in the supplement
# python qini_vs_autoc.py --scoring_types IPW AIPW Oracle --n_sims 10000
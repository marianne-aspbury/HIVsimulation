import sys
sys.path.insert(1, '/Users/anthonywohns/Documents/mcvean_group/age_inference/tsdate')

import tsdate
import tskit
import tsinfer

ts = tskit.load('truncated_simulation_tree.trees')

tip_weights = tsdate.find_node_tip_weights_ts(ts)
prior = tsdate.make_prior(ts.num_samples, 10000)
mixture_prior = tsdate.get_mixture_prior_ts_new(tip_weights, prior)
tsdate.age_inference(ts)

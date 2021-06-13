#!/usr/bin/env python
import json
import munch
import os
import logging
import numpy as np

import equations as eqn
import PA_plot as pap

from absl import app
from absl import flags
from absl import logging as absl_logging

flags.DEFINE_string('config_path', '/Users/fahim/Desktop/Vijay_desktop/OOPA/configs/PA_configs.json',
                    """The path to load json file.""")
flags.DEFINE_string('config_path_', '/Users/fahim/Desktop/Vijay_desktop/OOPA/configs/PA_configs.json',
                    """The path to load json file.""")

FLAGS = flags.FLAGS

def main(argv):
    del argv
    with open(FLAGS.config_path) as json_data_file:
        config = json.load(json_data_file)
    with open(FLAGS.config_path_) as json_data_file:
        config_ = json.load(json_data_file)
    config = munch.munchify(config)
    config_ = munch.munchify(config_)
    PA = getattr(eqn, config.eqn_config.eqn_name)(config)
    print('[mu, gamma, r ,lmb ,sgm, rho]: ',PA.param)
    PA.bvp_S()
    print(PA.message)
    PA_ = getattr(eqn, config.eqn_config.eqn_name)(config_)
    print('[mu, gamma, r ,lmb ,sgm, rho]: ',PA_.param)
    PA_.bvp_S()
    print(PA_.message)
    x = pap.plotF([PA,PA_])
    x.plot_iT(' ')

if __name__ == '__main__':
    app.run(main)

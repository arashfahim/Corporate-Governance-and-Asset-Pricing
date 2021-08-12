#!/usr/bin/env python
import json
import munch
import os
import logging
import ast
import numpy as np

import equations as eqn
import PA_plot as pap
import gui

from absl import app
from absl import flags
from absl import logging as absl_logging

# flags.DEFINE_string('config_path', '/Users/fahim/Desktop/Vijay_desktop/OOPA/configs/PA_configs.json',
#                     """The path to load json file.""")
# flags.DEFINE_string('config_path_', '/Users/fahim/Desktop/Vijay_desktop/OOPA/configs/PA_configs.json',
#                     """The path to load json file.""")

flags.DEFINE_string('config_path', "/Users/fahim/Desktop/Vijay_desktop/OOPA/configs/",
                    """The path to write json file.""")

FLAGS = flags.FLAGS

def main(argv):
    del argv
    pa_gui = gui.ParamEntry()
    pa_gui.initiate()
    pa_gui.window.protocol("WM_DELETE_WINDOW", pa_gui.on_closing)
    pa_gui.window.mainloop()
    param_list = pa_gui.dict_list
    PA_list=[]
    for indx, config_dict in enumerate(param_list):
        with open(FLAGS.config_path + "PA_configs" + str(indx) + ".json", "w") as outfile:
            json.dump(config_dict, outfile, indent=1)
        print("Configuration is saved in " + FLAGS.config_path + "PA_configs" + str(indx) + ".json")
        config = munch.munchify(config_dict)
        PA_tmp = getattr(eqn, "eqn_FfS")(config)
        print('[mu, gamma, r ,lmb ,sgm, rho]: ',PA_tmp.param)
        PA_tmp.bvp_S()
        print(PA_tmp.message)
        PA_list.append(PA_tmp)
    x = pap.plotF(PA_list)
    x.plot_iT(' ')
    # with open(FLAGS.config_path) as json_data_file:
    #     config = json.load(json_data_file)
    # with open(FLAGS.config_path_) as json_data_file:
    #     config_ = json.load(json_data_file)
    # config = munch.munchify(config)
    # config_ = munch.munchify(config_)



    # PA = getattr(eqn, config.eqn_config.eqn_name)(config)
    # print('[mu, gamma, r ,lmb ,sgm, rho]: ',PA.param)
    # PA.bvp_S()

    # PA_ = getattr(eqn, config.eqn_config.eqn_name)(config_)
    # print('[mu, gamma, r ,lmb ,sgm, rho]: ',PA_.param)
    # PA_.bvp_S()
    # print(PA_.message)
    # x = pap.plotF([PA,PA_])
    # x.plot_iT(' ')

if __name__ == '__main__':
    app.run(main)

#!/usr/bin/env python3
import matplotlib
import argparse

from plot_modes_bp import plot_all_modes

matplotlib.use("Agg")
from pyredlum import pyRedLUM
from plot_bias_bp import plot_bias
import plot_modes_bp
import os

##
def parse_time(value):
    try:
        parts = value.split(',')
        if len(parts) == 1:
            # One number: 0 to X
            x = float(parts[0])
            return (0, x)
        elif len(parts) == 2:
            # Two numbers : X to Y
            x, y = map(float, parts)
            if x > y:
                raise ValueError("First number needs to be equal or higher than the first")
            return (x, y)
        else:
            raise ValueError("Incorrect format, use 'X' or 'X,Y'.")
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Wrong argument '{value}': {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="redlum-plot",
        description="redlum-plot is a script that is using the pyRedLUM framework allowing to quickly produce plots of the\
        RedLUM algorithm"
    )
    parser.add_argument("-v","--verbose",type=int,default=0,help="increase output verbosity")
    parser.add_argument("-n","--name",type=str,default=None,help="Specify name used when saving plots")
    parser.add_argument("-t","--time",type=parse_time,default=None,help="Specify the time range (if only one value, starts from 0)")
    working_dir = os.getcwd()
    plotting_dir = f"{working_dir}"

    # Getting arguments stored in parser
    args = parser.parse_args()
    verbose = args.verbose
    name = args.name
    time = args.time

    if name is None:
        name = os.getcwd().split("/")[-1]

    case = pyRedLUM(
        res_folder=working_dir,
        save_dir=plotting_dir,
        name = name,
        verbose=verbose,
        time=time
    )

    plot_bias(case)
    plot_all_modes(case)
##

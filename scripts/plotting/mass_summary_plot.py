import argparse
from utilities.io_tools import combinetf_input
from utilities import common
from wremnants import plot_tools
import pandas as pd

parser = argparse.ArgumentParser()
parser = common.plot_parser()
parser.add_argument("-n", "--nominal", required=True, type=str, help="Combine fitresult file for nominal result")
parser.add_argument("-m", "--mll", required=True, type=str, help="Combine fitresult file for mll result")
parser.add_argument("-l", "--flip", required=True, type=str, help="Combine fitresult file for even-odd flip result")

args = parser.parse_args()

pdg_massz = 91188
ref_massz = 91187.6

ref_massw = 80385
pdg_massw = 80385

def read_res_df(filename, name):
    fitresult = combinetf_input.get_fitresult(filename)
    df = combinetf_input.read_impacts_pois(fitresult, poi_type="nois", group=True, uncertainties=["stat",])

    df.iloc[0,1:] = df.iloc[0,1:]*100
    df.iloc[0,1] += ref_massz if df.loc[0, "Name"] == "massShiftZ100MeV_noi" else pdg_massw
    df.loc[0, "Name"] = name

    df.rename(columns={"err_stat" : "Stat unc."}, inplace=True)
    return df

nominal_df = read_res_df(args.nominal, "Nominal W-like $m_{Z}$ fit")
mll_df = read_res_df(args.mll, "$m_{\\ell\\ell}$ fit")
flip_df = read_res_df(args.flip, "W-like $m_{Z}$ fit (even $\\leftrightarrow$ odd)")

dfs = pd.concat([mll_df, nominal_df, flip_df])
print(dfs)

plot_tools.make_summary_plot(pdg_massz, 2.0, "PDG average",
    dfs,
    colors="auto",
    xlim=[91160, 91220],
    xlabel="$m_{Z}$ (MeV)",
    out=args.outpath, outfolder=args.outfolder, name="Wlike_unblinding_summary",
    offset=2,
    fontsize=22,
    #center_colors=["black"]*len(dfs),
)

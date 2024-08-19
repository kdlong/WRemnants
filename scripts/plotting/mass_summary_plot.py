import argparse
from utilities.io_tools import combinetf_input
from utilities import common
from wremnants import plot_tools

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

def read_res_df(filename):
    fitresult = combinetf_input.get_fitresult(filename)
    df = combinetf_input.read_impacts_pois(fitresult, poi_type="nois", group=True, uncertainties=["stat",])

    df.iloc[0,1:] = df.iloc[0,1:]*100
    df.iloc[0,1] += ref_massz if df.loc[0, "Name"] == "massShiftZ100MeV_noi" else pdg_massw

    df.rename(columns={"err_stat" : "Stat unc."}, inplace=True)
    return df

nominal_df = read_res_df(args.nominal)
mll_df = read_res_df(args.mll)
flip_df = read_res_df(args.flip)

label = "massShiftZ100MeV_noi" 

plot_tools.make_summary_plot(pdg_massz, 2.0, "PDG average",
    [nominal_df, flip_df, mll_df, None, None], 
    colors="auto",
    labels=["W-like even $\\leftrightarrow$ odd fit", "Nominal W-like fit", "$m_{\\ell\\ell}$ fit", None, None],
    xlim=[91160, 91220] if label == "massShiftZ100MeV_noi" else [80360, 80410],
    xlabel="$m_{Z}$ (MeV)",
    out=args.outpath, outfolder=args.outfolder, name="Wlike_unblinding_summary",
)

import pandas as pd
import ROOT
from scipy.stats import chi2

from utilities import logging, parsing
from utilities.io_tools import output_tools, tex_tools

translate = {}


def read_fitresult(filename):
    try:
        rfile = ROOT.TFile.Open(filename)
        ttree = rfile.Get("fitresults")
        ttree.GetEntry(0)

        if hasattr(ttree, "massShiftW100MeV"):
            m = ttree.massShiftW100MeV
            merr = ttree.massShiftW100MeV_err
        else:
            m = 0
            merr = 0

        val = 2 * (ttree.nllvalfull - ttree.satnllvalfull)
        ndf = rfile.Get("obs;1").GetNbinsX() - ttree.ndofpartial
        p = (1 - chi2.cdf(val, ndf)) * 100

        status = ttree.status
        errstatus = ttree.errstatus
        edmval = ttree.edmval

    except IOError as e:
        return 0, 1, 0, 0, 0, 0, 0, 0

    return val, ndf, p, status, errstatus, edmval, m, merr


parser = parsing.plot_parser()
parser.add_argument("inputs", nargs="+", type=str, help="Paths to fitresult files")
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

df = pd.DataFrame(args.inputs, columns=["path"])

df[
    ["chi2", "ndf", "pvalue", "status", "errstatus", "edmval", "mass_obs", "mass_err"]
] = (df["path"].apply(read_fitresult).apply(pd.Series))

df["name_parts"] = df["path"].apply(
    lambda x: [y for y in filter(lambda z: z, x.split("/"))]
)

df["channel"] = df["name_parts"].apply(lambda x: x[-2].split("_")[0])

df["column_name"] = df["name_parts"].apply(lambda x: x[-2].split("_")[-1])
df["column_name_ndf"] = df["column_name"] + df["ndf"].apply(lambda x: f" ({x})")

df["dataset"] = df["name_parts"].apply(
    lambda x: "_".join(x[-1].split("_")[2:]).replace(".root", "")
)

order = ["simple", "extended1D", "extended2D", "closure"]
cat_dtype = pd.CategoricalDtype(categories=order, ordered=True)
df["column_name"] = df["column_name"].astype(cat_dtype)
df["dataset"] = df["dataset"].astype(cat_dtype)

df = df.sort_values(by=["column_name", "dataset"])

for channel, df_c in df.groupby("channel"):
    tex_tools.make_latex_table(
        df_c,
        output_dir=outdir,
        output_name=f"table_{channel}",
        column_title="Pseudodata",
        caption=r"Resulting $\chi^2$ values (and p-values) using the saturated model test from fits on different data, and pseudodata sets.",
        label="Model",
        sublabel="",
        column_name="dataset",
        row_name="column_name",
        cell_columns=["chi2", "pvalue"],
        color_condition=lambda x, y: y < 5,
        cell_format=lambda x, y: rf"${round(x,1)} ({round(y,1)}\%)$",
    )

    tex_tools.make_latex_table(
        df_c,
        output_dir=outdir,
        output_name=f"table_mass_{channel}",
        column_title="Pseudodata",
        caption="Mass and uncertainty.",
        label="Model",
        sublabel="",
        column_name="dataset",
        row_name="column_name",
        cell_columns=["mass_obs", "mass_err"],
        color_condition=lambda x, y: abs(x) > y,
        cell_format=lambda x, y: rf"${round(x*100,2)}\, \pm {round(y*100,2)}$",
    )

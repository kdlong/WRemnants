import argparse
from utilities.io_tools import combinetf_input
from utilities import common
from wremnants import plot_tools,theory_tools

parser = argparse.ArgumentParser()
parser = common.plot_parser()
parser.add_argument("-r", "--reffile", required=True, type=str, help="Combine fitresult file for nominal result")
parser.add_argument("--pdfs", default=["ct18z", "ct18", "ct18z", "herapdf20",  "msht20",  "msht20an3lo",  "nnpdf31",  "nnpdf40",  "pdf4lhc21"],
    type=str, help="PDF to plot")
parser.add_argument("--print", action='store_true', help="Print results")

args = parser.parse_args()

pdg_massz = 91188
pdg_unc = 2.0

pdf_name = lambda p: theory_tools.pdfMap[p]["name"]

pdf_dfs = {pdf_name(pdf)[3:] : combinetf_input.read_groupunc_df(args.reffile.format(pdf=pdf), pdf_name(pdf)) for pdf in args.pdfs}
print(pdf_dfs)

label = "massShiftZ100MeV_noi"
outname = "Wlike_pdf_summary"
if args.postfix:
    outname += f"_{args.postfix}"

if args.print:
    for k,v in pdf_dfs.items():
        print(k, round(v.iloc[0,1], 1), round(v.iloc[0,3], 1) , round(v.iloc[0,2], 1))

plot_tools.make_summary_plot(pdg_massz, pdg_unc, "PDG average",
    pdf_dfs.values(), 
    colors="auto",
    labels=pdf_dfs.keys(),
    xlim=[91160, 91220] if label == "massShiftZ100MeV_noi" else [80360, 80410],
    xlabel="$m_{Z}$ (MeV)",
    out=args.outpath, outfolder=args.outfolder, name=outname,
)

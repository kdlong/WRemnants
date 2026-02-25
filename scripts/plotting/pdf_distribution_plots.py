import argparse

import numpy as np
import matplotlib.pyplot as plt
import lhapdf
import os

from wums import output_tools, plot_tools
from wremnants import theory_tools
from rabbit import io_tools

PARTON_FLAVOR_NAMES = {
    "uv": "u_{V}",
    "1": "d",
    "-1": r"\bar{d}",
    "2": "u",
    "-2": r"\bar{u}",
    "3": "s",
    "-3": r"\bar{s}",
    "dv": "d_{v}",
    "rs": "r_{s}",
}

def read_pdf_vals_and_errors(x_range, flavor, Q_scale, pdf_sets):
    values = []
    errors = []
    for i, name in enumerate(pdf_sets):
        print(name)
        vals = theory_tools.get_pdf_data(name, flavor, Q_scale, x_range)
        central = vals[0]
        err = np.sqrt(np.sum((vals[1:] - central) ** 2, axis=0))
        values.append(vals[0])
        errors.append(err)

    return values,errors

def read_vals_and_errors_from_fit(fitresult_file, fit_types):
    fitresult = io_tools.get_fitresult(fitresult_file)
    values = []
    errors = []
    for fit in fit_types:
        h = fitresult["mappings"]["BaseMapping"]["channels"]["ch0"][f"hist_{fit}_inclusive"].get()
        values.append(h.values())
        errors.append(np.sqrt(h.variances()))
    return values,errors

def make_pdf_plot(values, errors, x_range, flavor, Q_scale, labels, colors, outdir, args):
    # --- Setup Plot ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, 
                                   gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    
    reference_central = values[0]

    for central, err, label, color in zip(values, errors, labels, colors):
        # 1. Main Plot (Top Panel)
        ax1.plot(x_range, central, color=color, label=label)
        ax1.fill_between(x_range, central - err, central + err, color=color, alpha=0.2)
        
        # 2. Ratio Plot (Bottom Panel)
        ratio_central = central / reference_central
        ratio_err = err / reference_central
        
        ax2.plot(x_range, ratio_central, color=color)
        ax2.fill_between(x_range, ratio_central - ratio_err, ratio_central + ratio_err, 
                         color=color, alpha=0.2)
    
    # Formatting Top Panel
    flav_label = PARTON_FLAVOR_NAMES.get(str(flavor), flavor)
    ax1.set_ylabel(f'$x {flav_label}(x, Q^2)$', fontsize=16)
    ax1.set_title(f'PDF at $Q = {Q_scale}$ GeV', fontsize=14)
    ax1.legend(loc="upper left")
    ax1.grid(True, which="both", alpha=0.3)
    
    # Formatting Ratio Panel
    ax2.axhline(1.0, color='black', lw=1, ls='--') # Reference line
    ax2.set_ylabel('Ratio to central', fontsize=14)
    ax2.set_xlabel(r'$x$', fontsize=12)
    ax2.set_xscale('log')
    ax2.set_ylim(0.8, 1.2) # Typically +/- 20% for ratio plots
    ax2.grid(True, which="both", alpha=0.3)
    
    outfile = f"pdf_{flavor}_Q{int(Q_scale)}"
    if args.postfix:
        outfile += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outdir, outfile)
    output_tools.write_index_and_log(
        outdir,
        outfile,
        args=args,
    )

    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate PDF plots from Rabbit tensor")
    parser.add_argument("-p", "--postfix", help="Label to append to the plot name")
    parser.add_argument(
        "-s", "--pdf-sets", nargs="+", required=True, help="LHAPDF set names"
    )
    parser.add_argument("-l", "--labels", nargs="+", help="Labels for the legend")
    parser.add_argument(
        "-f", "--flavors", nargs="+", help="Flavors (uv, dv, rs, or PDG ID)"
    )
    parser.add_argument(
        "-q", "--q-scale", type=float, default=80.360, help="Q scale in GeV"
    )
    parser.add_argument("-o", "--outpath", required=True, help="Output filename")
    parser.add_argument(
        "--lhapdf-path",
        default="/scratch/submit/cms/wmass/PostfitPDF/",
        help="Additional path to LHAPDF data files (for custom sets)",
    )
    parser.add_argument("--colors", nargs="+", help="List of colors for the sets")
    parser.add_argument("--fromFit", action='store_true', help="Plot from Rabbit fit")
    parser.add_argument("--fitResultFile", type=str, help="Result of Rabbit fit")
    parser.add_argument("--fitTypes", nargs="+", help="Prefit and/or postfit results from Rabbit fit")

    args = parser.parse_args()

    if args.fromFit and not args.fitResultFile:
        parser.error("--fromFit requires --fitResultFile")
    if args.fromFit and not args.fitTypes:
        parser.error("--fromFit requires --fitTypes")

    # Set LHAPDF path 
    lhapdf.pathsAppend(args.lhapdf_path)
    
    # If labels aren't provided, use PDF set names
    labels = args.labels if args.labels else args.pdf_sets

    outdir = output_tools.make_plot_dir(args.outpath, "", eoscp=True)

    x_range = np.logspace(-4, -0.01, 201)[:-1]

    for flavor in args.flavors:
        if args.fromFit :
            vals,errors = read_pdf_vals_and_errors(
                x_range,
                flavor=flavor,
                Q_scale=args.q_scale,
                pdf_sets=args.pdf_sets,
            )             
            
            vals_temp,errs_temp = read_vals_and_errors_from_fit( 
                fitresult_file=args.fitResultFile,
                fit_types=args.fitTypes,
            )
            vals.extend(vals_temp)
            errors.extend(errs_temp)  
        else :

            vals,errors = read_pdf_vals_and_errors(
                x_range,
                flavor=flavor,
                Q_scale=args.q_scale,
                pdf_sets=args.pdf_sets,
            )

        make_pdf_plot(
            vals,
            errors,
            x_range,
            flavor=flavor,
            Q_scale=args.q_scale,
            labels=labels,
            colors=args.colors,
            outdir=outdir,
            args=args,
        )     

    if output_tools.is_eosuser_path(args.outpath):
        output_tools.copy_to_eos(outdir, args.outpath, "")


if __name__ == "__main__":
    main()
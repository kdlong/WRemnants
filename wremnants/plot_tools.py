import mplhep as hep
import matplotlib.pyplot as plt
from matplotlib import patches
from wremnants import boostHistHelpers as hh
from wremnants import histselections as sel
import math
import numpy as np
import re

hep.style.use(hep.style.ROOT)

def figureWithRatio(href, xlabel, ylabel, ylim, rlabel, rrange, xlim=None,
    grid_on_main_plot = False, grid_on_ratio_plot = False, plot_title = None
):
    hax = href.axes[0]
    width = math.ceil(hax.size/400)
    fig = plt.figure(figsize=(16*width,16))
    ax1 = fig.add_subplot(4, 1, (1, 3)) 
    ax2 = fig.add_subplot(4, 1, 4) 

    ax2.set_xlabel(xlabel, fontsize=60)
    ax1.set_xlabel(" ")
    ax1.set_ylabel(ylabel, fontsize=60)
    ax1.set_xticklabels([], fontsize=60)
    if not xlim:
        xlim = [href.axes[0].edges[0], href.axes[0].edges[href.axes[0].size-1]]
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)
    ax2.set_ylabel(rlabel, fontsize=60)
    ax2.set_ylim(rrange)
    if ylim:
        ax1.set_ylim(ylim)
    else:
        ax1.autoscale(axis='y')
    if grid_on_main_plot:  ax1.grid(which = "both")
    if grid_on_ratio_plot: ax2.grid(which = "both")
    if plot_title: ax1.set_title(plot_title)
    ax1.tick_params(axis='both', which='major', labelsize=50)
    ax2.tick_params(axis='both', which='major', labelsize=50)
    return fig,ax1,ax2

def addLegend(ax, ncols=2, extra_text=None):
    has_extra_text = extra_text is not None
    handles, labels = ax.get_legend_handles_labels()
    
    if has_extra_text:
        #handles.append(patches.Patch(color='none', label=extra_text))
        ax.plot([], [], ' ', ' ')

    shape = np.divide(*ax.get_figure().get_size_inches())
    #TODO: The goal is to leave the data in order, but it should be less hacky
    handles[:] = reversed(handles)
    labels[:] = reversed(labels)
    if len(handles) % 2 and ncols == 2:
        handles.insert(math.floor(len(handles)/2), patches.Patch(color='none', label = ' '))
        labels.insert(math.floor(len(labels)/2), ' ')
    #handles= reversed(handles)
    #labels= reversed(labels)
    ax.legend(handles=handles, labels=labels, prop={'size' : 40*(0.7 if shape == 1 else 1.3)}, ncol=ncols, loc='upper right')

def makeStackPlotWithRatio(
    histInfo, stackedProcs, histName="nominal", unstacked=None, 
    xlabel="", ylabel="Events/bin", rrange=[0.9, 1.1], ymax=None, xlim=None, nlegcols=2,
    binwnorm=None, select={},  action = (lambda x: x), extra_text=None, grid = False, plot_title = None
):
    stack = [action(histInfo[k][histName][select]) for k in stackedProcs if histInfo[k][histName]]
    colors = [histInfo[k]["color"] for k in stackedProcs if histInfo[k][histName]]
    labels = [histInfo[k]["label"] for k in stackedProcs if histInfo[k][histName]]
    fig, ax1, ax2 = figureWithRatio(stack[0], xlabel, ylabel, [0, ymax] if ymax else None, "Data/Pred.", rrange, xlim=xlim, grid_on_ratio_plot = grid, plot_title = plot_title)
    
    hep.histplot(
        stack,
        histtype="fill",
        color=colors,
        label=labels,
        stack=True,
        ax=ax1,
        binwnorm=binwnorm,
    )
    
    if unstacked:
        if type(unstacked) == str: unstacked = unstacked.split(",")
        for i, proc in enumerate(unstacked):
            unstack = action(histInfo[proc][histName][select])
            isdata = "data" in proc.lower()
            hep.histplot(
                unstack,
                yerr=True if isdata else False,
                histtype="errorbar" if isdata else "step",
                color=histInfo[proc]["color"],
                label=histInfo[proc]["label"],
                ax=ax1,
                zorder=1 if i == 0 else -1,
                binwnorm=binwnorm,
            )
            hep.histplot(
                hh.divideHists(unstack, sum(stack), cutoff=0.01),
                histtype="errorbar" if isdata else "step",
                color=histInfo[proc]["color"],
                label=histInfo[proc]["label"],
                yerr=True if isdata else False,
                ax=ax2
            )

    addLegend(ax1, nlegcols, extra_text)
    return fig

def makePlotWithRatioToRef(
    hists, labels, colors, xlabel="", ylabel="Events/bin", rlabel="x/nominal",
    rrange=[0.9, 1.1], ymax=None, xlim=None, nlegcols=2, binwnorm=None, alpha=1.,
    baseline=True, data=None, data_label=None, autorrange=None, grid = False, fill_between=False
):
    # nominal is always at first, data is always at last, if included
    ratio_hists = [hh.divideHists(h, hists[0], cutoff=0.00001) for h in hists[not baseline:]]
    fig, ax1, ax2 = figureWithRatio(hists[0], xlabel, ylabel, [0, ymax] if ymax else None, rlabel, rrange, xlim=xlim, grid_on_ratio_plot = grid)
    
    print("WE ARE HERE!")

    hep.histplot(
        hists,
        histtype="step",
        color=colors,
        label=labels,
        stack=False,
        ax=ax1,
        binwnorm=binwnorm,
        alpha=alpha,
    )
    
    if data:
        hep.histplot(
            data,
            histtype="errorbar",
            color="black",
            label=data_label,
            stack=False,
            ax=ax1,
            binwnorm=binwnorm,
        )
        hep.histplot(
            hh.divideHists(data, hists[0], cutoff=1.e-6),
            histtype="errorbar",
            color="black",
            label=data_label,
            xerr=False,
            yerr=False,
            stack=False,
            ax=ax2,
        )

    endidx = len(hists) if not fill_between else 1
    hep.histplot(
        ratio_hists[not baseline:endidx],
        histtype="step",
        color=colors[not baseline:endidx],
        label=labels[not baseline:endidx],
        yerr=False,
        stack=False,
        ax=ax2,
        alpha=alpha,
    )
    if fill_between and not len(hists) % 2:
        print("Fill")
        for h1,h2,color in zip(ratio_hists[1::2], ratio_hists[2::2], colors[1::2]):
            ax2.fill_between(h1.axes[0].edges[:-1], h1.values(), h2.values(), step='mid', color=color, alpha=alpha)
         
    addLegend(ax1, nlegcols)
    return fig

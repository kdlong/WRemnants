from scripts.plotting import makeDataMCStackPlot
from scripts.combine import setupCombine,pullsAndImpacts
import os

args = dict(
    hists=["ptll"],
    #axlim=[0, 20],
    baseName="nominal",
    rrange=[0.9, 1.1],
    eoscp=True,
    outpath='/eos/user/k/kelong/www/WMassAnalysis/2024May_TNPStudy/Z',
    selectAxis=["vars"], 
    fillBetween=2,
    doubleColors=False,
    yscale=1.3,
)


corrs = ['scetlib_dyturboN3p0LLCorr', 'scetlib_dyturboN2p1LLCorr','scetlib_dyturboN3p1LLCorr', 'scetlib_dyturboN4p0LLCorr',]
#corrs = ['scetlibN3p0LLCorr', 'scetlibN2p1LLCorr','scetlibN3p1LLCorr', 'scetlibN4p0LLCorr',]
#corrs = ["scetlib_dyturboN3p0LLCorr", "dyturboN3LLpCorr", "scetlib_dyturboN4LLMSHT20an3loCorr", "matrix_radishCorr", "dyturboMSHT20Corr"]

style_map = {
    'scetlib_dyturboN3p0LLCorr' : {
        "color" : "mediumpurple",
        "label" : 'N$^{3+0}$LL+NNLO', 
    },
    'scetlib_dyturboN2p1LLCorr' : {
        "label" : 'N$^{2+1}$LL+NNLO', 
        "color" : "darkred",
    },
    'scetlib_dyturboN3p1LLCorr' : {
        "color" : "purple",
        "label" : 'N$^{3+1}$LL+NNLO', 
    },
    'scetlib_dyturboN4p0LLCorr' : {
        "color" : "darkgreen",
        "label" : 'N$^{4+0}$LL+NNLO',
    },
    'scetlibN3p0LLCorr' : {
        "color" : "indigo",
        "label" : 'N$^{3+0}$LL', 
    },
    'scetlibN2p1LLCorr' : {
        "color" : "red",
        "label" : 'N$^{2+1}$LL', 
    },
    'scetlibN3p1LLCorr' : {
        "color" : "darkviolet",
        "label" : 'N$^{3+1}$LL', 
    },
    'scetlibN4p0LLCorr' : {
        "color" : "forestgreen",
        "label" : 'N$^{4+0}$LL',
    },
    "dyturboN3LLpCorr" : {
        "color" : "darkgreen",
        "label" : "DYTurbo N$^{3}LL\prime$+NNLO",
    },
    "scetlib_dyturboN4LLMSHT20an3loCorr" : {
        "color" : "darkorange",
        "label" : "SCETlib+DYTurbo N$^{4}LL\prime$+NNLO (MSHT20an3lo)",
    },
    "matrix_radishCorr" : {
        "color" : "lightgreen",
        "label" : "MATRIX+RadISH N$^{3}LL$+NNLO",
    },
    "dyturboMSHT20Corr" : {
        "color" : "darkred",
        "label" : "DYTurbo N$^{3}LL\prime$+NNLO (MSHT20nnlo)",
    },
}

def other_variations(arr, pnom):
    return arr[:pnom]+arr[pnom+1:]

def select_variations(arr, pnom):
    return arr[pnom:pnom+1]*2+other_variations(arr, pnom)

for i,corr in enumerate(corrs):
    #histfile = f'/scratch/submit/cms/kdlong/Analysis/mz_dilepton_{corr}_altCorrs.hdf5'
    histfile = f'/scratch/submit/cms/kdlong/Analysis/mz_dilepton_{corr}.hdf5'
    print("-----> Nominal corr", corr, "histfile:", histfile)
    unc = f"resumTNPXp{'0' if 'p0LL' in corr else '1'}"

    varis = select_variations(corrs, i)
    labels = [style_map[v]["label"] for v in varis]
    colors = [style_map[v]["color"] for v in varis]
    labels[0] += " TNP $\pm 1\sigma$"
    labels[1] = ''
    #postfix=f"{corr}_AltPred"
    postfix=f"{corr}_TNPx2"
    postfix += "_0To20" if "xlim" in args else "fullrange"
    args.update(
        dict(
            infile=histfile,
            varName=varis,
            varLabel=labels,
            colors=colors,
            selectEntries=[unc+"Up", unc+"Down"]+["0"]*(len(corrs)-1),
            postfix=postfix,
            linestyle=['solid']*3+['dashed', 'dotted'],
            outfolder=corr,
            scaleleg=0.7,
            scaleTNP=2,
        )
    )
    #postfix = f"{corr}_pdTNPs_fullrange"
    #postfix=f"{corr}_Data"
    plotoutf = "/scratch/submit/cms/kdlong/CombineStudies/TheoryBiasStudies/"
    makeDataMCStackPlot.run_script(args)
    setupCombine.run_script(dict(
        inputFile=[histfile],
        outfolder=plotoutf,
        postfix=postfix,
        pseudoData=other_variations(corrs, i),
        fitvar=["ptll"],
        hdf5=True,
        )
    )

    fitfile = f"{plotoutf}/ZMassDilepton_ptll_{postfix}/fitresults_123456789.hdf5"
    if os.path.isfile(fitfile):
        pullsAndImpacts.run_script(dict(
            mode='ungrouped',
            oneSidedImpacts=True,
            filters=['.*scetlibNP|resumTNP.*'],
            inputFile=fitfile,
            output_mode='output',
            outFolder=os.path.join(args['outpath'], corr), 
            outputFile='impacts_ptll0To20.html',
            eoscp=True,
            otherExtensions=['pdf', 'png'],
            num=None,
            noPulls=False,
            )
        )

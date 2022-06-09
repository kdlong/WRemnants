from wremnants import boostHistHelpers as hh
from wremnants import histselections as sel
from wremnants.datasets import datasets2016
import logging
import lz4.frame
import pickle
import narf
#import uproot
import ROOT
import re

class datagroups(object):
    def __init__(self, infile, combine=False):
        self.combine = combine
        self.lumi = 1.
        if ".root" not in infile[-5:]:
            with lz4.frame.open(infile) as f:
                self.results = pickle.load(f)
            self.rtfile = None
        else:
            self.rtfile = ROOT.TFile.Open(infile)
            self.results = None

        self.lumi = None
        if self.datasets and self.results:
            self.data = [x for x in self.datasets.values() if x.is_data]
            if self.data and not combine:
                self.lumi = sum([self.results[x.name]["lumi"] for x in self.data if x.name in self.results])
        self.groups = {}

        self.nominalName = "nominal"

    def setLumi(self, lumi):
        self.lumi = lumi

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]

    def setHists(self, baseName, syst, procsToRead=None, label=None, nominalIfMissing=True, 
            selectSignal=True, forceNonzero=True):
        if not label:
            label = syst if syst else baseName
        if not procsToRead:
            procsToRead = self.groups.keys()

        foundExact = False
        for procName in procsToRead:
            group = self.groups[procName]
            group[label] = None

            for member in group["members"]:
                scale = group["scale"] if "scale" in group else None
                try:
                    h = self.readHist(baseName, member, syst, scaleOp=scale, forceNonzero=forceNonzero)
                    foundExact = True
                except ValueError as e:
                    if nominalIfMissing:
                        logging.warning(f"{str(e)}. Using nominal hist {self.nominalName} instead")
                        h = self.readHist(self.nominalName, member, "", scaleOp=scale, forceNonzero=forceNonzero)
                    else:
                        logging.warning(str(e))
                        continue
                group[label] = h if not group[label] else hh.addHists(h, group[label])
            if selectSignal and group[label] and "signalOp" in group and group["signalOp"]:
                group[label] = group["signalOp"](group[label])
        # Avoid situation where the nominal is read for all processes for this syst
        if not foundExact:
            raise ValueError(f"Did not find systematic {syst} for any processes!")

    #TODO: Better organize to avoid duplicated code
    def setHistsCombine(self, baseName, syst, channel, procsToRead=None, excluded_procs=[], label=None, nominalIfMissing=True):
        if type(excluded_procs) == str: excluded_procs = excluded_procs.split(",")
        #TODO Set axis names properly
        if baseName == "x":
            axisNames=["eta", "pt"]

        if not label:
            label = syst
        if not procsToRead:
            procsToRead = list(filter(lambda x: x not in excluded_procs, self.groups.keys()))

        for procName in procsToRead:
            group = self.groups[procName]
            group[label] = None
            h = self.readHistCombine(baseName, procName, syst, channel, axisNames, nominalIfMissing)
            group[label] = h 

    def histNameCombine(self, procName, baseName, syst, channel):
        name = f"{baseName}_{procName}"
        if syst not in ["", "nominal"]:
            name += "_"+syst
        if channel:
            name += "_"+channel
        if re.search("^pdf.*_sum", procName): # for pseudodata from alternative pdfset
            return("_".join([procName, channel])) 
        return name

    def loadHistsForDatagroups(self, baseName, syst, procsToRead=None, excluded_procs=None, channel="", label="", nominalIfMissing=True,
            selectSignal=True, forceNonzero=True, pseudodata=False):
        if self.rtfile and self.combine:
            self.setHistsCombine(baseName, syst, channel=channel, procsToRead=procsToRead, excluded_procs=excluded_procs, 
                    label=label, nominalIfMissing=nominalIfMissing)
        else:
            self.setHists(baseName, syst, procsToRead, label, nominalIfMissing, selectSignal, forceNonzero)

    def getDatagroups(self, excluded_procs=[]):
        if type(excluded_procs) == str: excluded_procs = list(excluded_procs)
        return dict(filter(lambda x: x[0] not in excluded_procs, self.groups.items()))

    def getDatagroupsForHist(self, histName):
        filled = {}
        for k, v in self.groups.items():
            if histName in v:
                filled[k] = v
        return filled

    def resultsDict(self):
        return self.results

    def processes(self):
        return self.groups.keys()

    def addSummedProc(self, refname, name, label, histname="", color="red", procsToRead=None, channel="", nominalIfMissing=True):
        if not histname:
            histname = name
        self.loadHistsForDatagroups(refname, syst=name, label=name, channel=channel, procsToRead=procsToRead, nominalIfMissing=nominalIfMissing)
        self.groups[name] = dict(
            label=label,
            color=color,
            members=[],
        )
        self.groups[name][histname] = sum([self.groups[x][name] for x in procsToRead])

    def copyWithAction(self, action, name, refproc, refname, label=None, color=None):
        self.groups[name] = dict(
            label=label if label else self.groups[refproc]["label"],
            color=color if color else self.groups[refproc]["color"],
            members=[],
        )
        self.groups[name][refname] = action(self.groups[refproc][refname])

class datagroups2016(datagroups):
    def __init__(self, infile, combine=False, wlike=False, pseudodata_pdfset = None):
        self.datasets = {x.name : x for x in datasets2016.getDatasets()}
        super().__init__(infile, combine)
        self.hists = {} # container storing temporary histograms
        self.groups =  {
            "Data" : dict(
                members = [self.datasets["dataPostVFP"]],
                color = "black",
                label = "Data",
                signalOp = sel.signalHistWmass if not wlike else None,
            ),
            "Zmumu" : dict(
                members = [self.datasets["ZmumuPostVFP"]],
                label = r"Z$\to\mu\mu$",
                color = "lightblue",
                signalOp = sel.signalHistWmass if not wlike else None,
            ),   
            "Ztautau" : dict(
                members = [self.datasets["ZtautauPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
                signalOp = sel.signalHistWmass if not wlike else None,
            ),            
        }
        if pseudodata_pdfset:
            self.groups[f"pdf{pseudodata_pdfset.upper()}_sum"] = dict(
                label = f"Pseudodata pdf{pseudodata_pdfset.upper()}",
                color = "dimgray"
            )
        if not wlike:
            self.groups.update({
                "Fake" : dict(
                    members = list(self.datasets.values()),
                    scale = lambda x: 1. if x.is_data else -1,
                    label = "Nonprompt",
                    color = "grey",
                    signalOp = sel.fakeHistABCD,
                ),
                "Wtau" : dict(
                    members = [self.datasets["WminustaunuPostVFP"], self.datasets["WplustaunuPostVFP"]],
                    label = r"W$^{\pm}\to\tau\nu$",
                    color = "orange",
                    signalOp = sel.signalHistWmass,
                ),
                "Wmunu" : dict(
                    members = [self.datasets["WminusmunuPostVFP"], self.datasets["WplusmunuPostVFP"]],
                    label = r"W$^{\pm}\to\mu\nu$",
                    color = "darkred",
                    signalOp = sel.signalHistWmass,
                ),
                "Top" : dict(
                    members = list(filter(lambda y: y.group == "Top", self.datasets.values())),
                    label = "Top",
                    color = "green",
                    signalOp = sel.signalHistWmass,
                ), 
                "Diboson" : dict(
                    members = list(filter(lambda y: y.group == "Diboson", self.datasets.values())),
                    label = "Diboson",
                    color = "pink",
                    signalOp = sel.signalHistWmass,
                ), 
            })
        else:
            self.groups["Other"] = dict(
                members = [x for x in self.datasets.values() if not x.is_data and x.name not in ["ZmumuPostVFP", "ZtautauPostVFP"]],
                label = "Other",
                color = "grey",
            )

    def histName(self, baseName, procName, syst):
        # This is kind of hacky to deal with the different naming from combine
        if baseName != "x" and (syst == "" or syst == self.nominalName):
            return baseName
        if (baseName == "" or baseName == "x") and syst:
            return syst
        return "_".join([baseName,syst])
    
    def readHist(self, baseName, proc, syst, scaleOp=None, forceNonzero=True):
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")
        h = output[histname]
        if forceNonzero:
            h = hh.clipNegativeVals(h)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale

    def readHistCombine(self, baseName, procName, syst, channel, axisNames, nominalIfMissing=True):
        name = self.histNameCombine(procName, baseName, syst, channel)
        rthist = self.rtfile.Get(name)
        if not rthist and nominalIfMissing:
            logging.warning(f"Failed to find hist {name} in ROOT file. Using nominal hist {self.nominalName} instead")
            name = self.histNameCombine(procName, baseName, "nominal", channel)
            rthist = self.rtfile.Get(name)
        if not rthist:
            raise RuntimeError(f"Failed to load hist {name} from file")
        return narf.root_to_hist(rthist, axis_names=axisNames)

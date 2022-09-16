import ROOT
import pathlib
import hist
import narf
from utilities import rdf_tools


ROOT.gInterpreter.Declare('#include "muon_calibration.h"')
ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def make_muon_calibration_helpers(filename = data_dir + "/calibration/correctionResults_v207_gen.root", era = None):

    helper = ROOT.wrem.CVHCorrector(filename)

    uncertainty_hist = get_dummy_uncertainties()
    uncertainty_hist_cpp = narf.hist_to_pyroot_boost(uncertainty_hist, tensor_rank = 2)
    # min gen pt = 9 GeV to avoid threshold effects
    # max weight = 10 to protect against outliers
    uncertainty_helper = ROOT.wrem.calibration_uncertainty_helper[type(uncertainty_hist_cpp)](ROOT.std.move(uncertainty_hist_cpp), 9., 10.)

    down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
    uncertainty_helper.tensor_axes = (uncertainty_hist.axes["calvar"], down_up_axis)

    return helper, uncertainty_helper

'''
def get_dummy_uncertainties():
    axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
    axis_calvar = hist.axis.Integer(0, 6, underflow=False, overflow=False, name = "calvar")
    axis_calparm = hist.axis.Integer(0, 4, underflow=False, overflow=False, name = "calparm")

    h = hist.Hist(axis_eta, axis_calvar, axis_calparm)

    # forward b-field-like
    h.values()[:axis_eta.index(-1.2), 0, 0] = 1e-4

    # forward alignment-like
    h.values()[:axis_eta.index(-1.2), 1, 2] = 1e-4/40.

    # forward resolution
    h.values()[:axis_eta.index(-1.2), 2, 3] = 0.05

    # central b-field-like
    h.values()[axis_eta.index(-1.2):axis_eta.index(0.), 3, 0] = 1e-4

    # central alignment-like
    h.values()[axis_eta.index(-1.2):axis_eta.index(0.), 4, 2] = 1e-4/40.

    # central resolution
    h.values()[axis_eta.index(-1.2):axis_eta.index(0.), 5, 3] = 0.05

    #mirror to positive eta
    h.values()[axis_eta.index(0.):, ...] = h.values()[:axis_eta.index(0.), ...]

    # set underflow and overflow to match boundaries
    h.values(flow=True)[0, ...] = h.values(flow=True)[1, ...]
    h.values(flow=True)[-1, ...] = h.values(flow=True)[-2, ...]

    return h
'''

def get_dummy_uncertainties():
    axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
    axis_calvar = hist.axis.Integer(0, 1, underflow=False, overflow=False, name = "calvar")
    axis_calparm = hist.axis.Integer(0, 4, underflow=False, overflow=False, name = "calparm")

    h = hist.Hist(axis_eta, axis_calvar, axis_calparm)

    # forward b-field-like
    h.values()[..., 0, 0] = 1e-4

    # set underflow and overflow to match boundaries
    h.values(flow=True)[0, ...] = h.values(flow=True)[1, ...]
    h.values(flow=True)[-1, ...] = h.values(flow=True)[-2, ...]

    return h

def define_corrected_muons(df, helper):
    # split the nested vectors
    df = df.Define("Muon_cvhmergedGlobalIdxs", "wrem::splitNestedRVec(Muon_cvhmergedGlobalIdxs_Vals, Muon_cvhmergedGlobalIdxs_Counts)")
    df = df.Define("Muon_cvhbsJacRef", "wrem::splitNestedRVec(Muon_cvhbsJacRef_Vals, Muon_cvhbsJacRef_Counts)")

    df = df.Define("Muon_correctedMom4Charge", helper, ["Muon_cvhbsPt", "Muon_cvhbsEta", "Muon_cvhbsPhi", "Muon_cvhbsCharge", "Muon_cvhmergedGlobalIdxs", "Muon_cvhbsJacRef"])

    # split into individual vectors
    df = df.Define("Muon_correctedPt", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Pt(); } ); return res;")
    df = df.Define("Muon_correctedEta", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Eta(); } ); return res;")
    df = df.Define("Muon_correctedPhi", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Phi(); } ); return res;")
    df = df.Define("Muon_correctedCharge", "ROOT::VecOps::RVec<int> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.second; }); return res;")
    return df

def define_gen_muons_from_reco_match(df, reco_subset = "goodMuons"):
    df = df.Define("goodMuons_idx", f"Muon_genPartIdx[{reco_subset}]")
    df = df.Define("GenMuon_pt", f"ROOT::VecOps::Take(GenPart_pt, goodMuons_idx)")
    df = df.Define("GenMuon_eta",f"ROOT::VecOps::Take(GenPart_eta, goodMuons_idx)")
    df = df.Define("GenMuon_phi", f"ROOT::VecOps::Take(GenPart_phi, goodMuons_idx)")
    df = df.Define("GenMuon_pdgId", f"ROOT::VecOps::Take(GenPart_pdgId, goodMuons_idx)")
    df = df.Define(
        "GenMuon_charge", (
            "ROOT::VecOps::RVec<int> res(GenMuon_pdgId.size());"
            "for (int i = 0; i < GenMuon_pdgId.size(); i++) {"
            "res[i] = GenMuon_pdgId[i] > 0 ? -1 : 1;}"
            "return res;"
        )
    )
    return df

def define_gen_muons_by_MC_truth(df, sort_key = None, sorted_cols_suffix = "_sorted"):
    df = df.Define(
        "isPrompt", 
        (
            "ROOT::VecOps::RVec<bool> res(GenPart_statusFlags.size());"
            "for (int i = 0; i < GenPart_statusFlags.size(); i++) {"
            "    res[i] = (GenPart_statusFlags[i] & 0x01);" # true when first bit is 1
            "}"
            "return res"
        )
    )
    df = df.Define("GenMuon_pt", "GenPart_pt[abs(GenPart_pdgId) == 13 && GenPart_status == 1 && isPrompt]")
    df = df.Define("GenMuon_eta", "GenPart_eta[abs(GenPart_pdgId) == 13 && GenPart_status == 1 && isPrompt]")
    df = df.Define("GenMuon_phi", "GenPart_phi[abs(GenPart_pdgId) == 13 && GenPart_status == 1 && isPrompt]")
    df = df.Define("GenMuon_pdgId", "GenPart_pdgId[abs(GenPart_pdgId) == 13 && GenPart_status == 1  && isPrompt]")
    df = df.Define("GenMuon_charge", (
        "ROOT::VecOps::RVec<int> res(GenMuon_pdgId.size());"
        "for (int i = 0; i < GenMuon_pdgId.size(); i++) {"
        "res[i] = GenMuon_pdgId[i] > 0 ? -1 : 1;}"
        "return res;"
    ))
    if sort_key:
        df = rdf_tools.sort_rdf_cols( # this adds cols_to_sort + "_sorted" to df
            df, sort_key = sort_key, cols_to_sort = [
                "GenMuon_pt",
                "GenMuon_eta",
                "GenMuon_phi",
                "GenMuon_charge"
            ],
            sorted_cols_suffix = sorted_cols_suffix
        )
    return df

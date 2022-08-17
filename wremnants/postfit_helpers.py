import numpy as np
import hist
import uproot
from wremnants import common, boostHistHelpers as hh
import re

def postfit_hist_to_2D(h):
    newh = hist.Hist(common.axis_eta, common.axis_pt, common.axis_charge, storage=h._storage_type())
    newh[...] = np.moveaxis(np.reshape(h.view(), (common.axis_charge.size, common.axis_eta.size, common.axis_pt.size)), 0, -1)
    return newh

class PostfitHelper(object):
    def __init__(self, fitfile, combinefile):
        self.fitfile = uproot.open(fitfile)
        self.combinefile = uproot.open(combinefile)
        self.cov = None
        self.nuisances = []
        self.axis_names = [x.name for x in common.nominal_axes[:2]]
        self.nominal_hist = None
        self.hist_name = "x_Wmunu"
        self.channels = ["plus", "minus"]
        self.nuisance_scales = {}

    def read_combine_hist_by_chan(self, syst, channel):
        hname = f"{self.hist_name}_{syst}_{channel}" if syst else "_".join([self.hist_name, channel])
        h = self.combinefile[hname].to_hist()
        if len(self.axis_names) != len(h.axes):
            raise ValueError("Mismatch between number of axes and names")

        for ax,name in zip(h.axes, self.axis_names):
            ax.__dict__['name'] = name

        h = self.add_charge_axis(h, channel)
        return h

    def read_combine_hist(self, syst):
        h = self.read_combine_hist_by_chan(syst, self.channels[0])
        for chan in self.channels[1:]:
            h += self.read_combine_hist_by_chan(syst, chan)
        return h

    def add_charge_axis(self, h, channel):
        hcharge = hist.Hist(*h.axes, common.axis_charge, storage=h._storage_type())
        charge = 0 if channel not in ["plus", "minus"] else (1. if "plus" else -1)

        # TODO: this wouldn't work for the Z at the moment
        hcharge[..., common.axis_charge.index(charge)] = h.view(flow=True)
        return hcharge

    def read_covariance_matrix(self):
        self.cov = self.fitfile["covariance_matrix_channelnois"].to_hist()

    def set_nuisance_subset(self, nuisance_names, order_from_cov=True):
        self.nuisances = nuisance_names 
        if order_from_cov:
            self.order_nuisances_from_cov()

    def set_nuisance_subset_by_regex(self, match_expr):
        self.nuisances = self.match_nuisances_regex(match_expr)
        self.order_nuisances_from_cov()

    def match_nuisances_regex(self, match_expr):
        ax = self.cov.axes["xaxis"]
        return [ax.value(i) for i in range(ax.size) if re.match(match_expr, ax.value(i)) ]

    def set_nuisance_scaling(self, nuisance_expr, scale):
        for nuisance in self.match_nuisances_regex(nuisance_expr):
            self.nuisance_scales[nuisance] = scale
         
    def order_nuisances_from_cov(self):
        if not self.cov:
            self.cov = read_covariance_matrix()
        idxs = self.cov.axes['xaxis'].index(self.nuisances)
        minidx = np.min(idxs)
        maxidx = np.max(idxs)

        if not all(np.sort(idxs) == np.arange(minidx, maxidx+1)):
            raise ValueError("Must select contiguous uncertainties")

        self.nuisances = [self.cov.axes['xaxis'][i] for i in range(minidx, maxidx+1)]
        self.cov_subset = self.cov.values()[minidx:maxidx+1, minidx:maxidx+1]

    def pulls(self):
        pulls_arrays = self.fitfile["fitresults"].arrays(self.nuisances)
        return np.array([pulls_arrays[x][0] for x in self.nuisances])

    def syst_hist(self, isUp):
        if not self.nominal_hist:
            raise ValueError("Read nominal hist before trying to load systematic hists")
        
        var_hists = []
        for x in self.nuisances:
            var_name = x+"Up" if isUp else x+"Down"
            #var_name = x+"Up"
            if x not in self.nuisance_scales:
                var_hists.append(self.read_combine_hist(var_name))
            else:
                scale = self.nuisance_scales[x]
                print("Scaling", x, "by", scale)
                var_hist = self.read_combine_hist(var_name)
                hratio = hh.divideHists(var_hist, self.nominal_hist, cutoff=0.0001)
                var_hist.values(flow=True)[...] = self.nominal_hist.values(flow=True)*np.exp(scale*np.log(hratio.values(flow=True)))
                var_hists.append(var_hist)
                
        syst_axis = hist.axis.Integer(0, len(var_hists)+1, name="systIdx", flow=False)
        var_hist = hist.Hist(*self.nominal_hist.axes, syst_axis)
        var_hist[...,0] = self.nominal_hist.values()
        var_hist[...,1:] = np.stack([x.values() for x in var_hists], axis=-1)
        print("Sum is", var_hist.sum())
        return var_hist

    def set_hist_name(self, hist_name):
        self.hist_name = hist_name

    def set_channels(self, channels):
        self.channels = channels

    def postfit_variation_hists(self):

        if self.cov_subset is None:
            self.cov_subset = self.cov.values()

        w,v = np.linalg.eigh(self.cov_subset)

        pulls = self.pulls()
        postfit_up = self.postfit_variations(self.syst_hist(isUp=True), w, v, pulls, isUp=True)
        postfit_down = self.postfit_variations(self.syst_hist(isUp=False), w, v, pulls, isUp=False)

        return postfit_up, postfit_down

    def read_nominal_hist(self):
        self.nominal_hist = self.read_combine_hist("")

    def prefit_tot_uncertainty_hists(self):
        up = self.syst_hist(isUp=True)
        down = self.syst_hist(isUp=False)
        up_var = hh.addHists(up[{"systIdx" : 0}], hh.rssFromVariationHist(up, "systIdx"))
        down_var = hh.addHists(down[{"systIdx" : 0}], -1*hh.rssFromVariationHist(down, "systIdx"))
        return (up_var, down_var)

    def postfit_tot_uncertainty_hists(self):
        up, down = self.postfit_variation_hists()
        up_var = hh.addHists(up[{"systIdx" : 0}], hh.rssFromVariationHist(up, "systIdx"))
        down_var = hh.addHists(down[{"systIdx" : 0}], -1*hh.rssFromVariationHist(down, "systIdx"))
        return (up_var, down_var)

    def postfit_variations(self, var_hist, w, v, pulls, isUp):
        var_values = var_hist.values()[...,1:]
        out = np.ones_like(var_values)
        kappa = np.divide(var_values, var_hist.values()[...,:1], out=out, where=var_values > 0.0001)

        fac = 1 if isUp else -1
        theta = np.insert(pulls + fac*np.sqrt(w)*v, 0, values=pulls, axis=0)

        pull_and_constraint = np.prod(np.power(kappa[...,np.newaxis].T, theta[...,np.newaxis,np.newaxis,np.newaxis]), axis=1).T
        new_var_hist = var_hist.copy() 
        new_var_hist[...] = var_hist[...,:1].values()*pull_and_constraint
        return new_var_hist



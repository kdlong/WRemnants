import uproot
import argparse
from utilities.io_tools import combinetf_input
import re

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--ungroup", action='store_true', help="Use ungrouped nuisances")
    parser.add_argument("-n", "--nuisance", type=str, help="Only print value for specific nuiance")
    parser.add_argument("-s", "--sort", action='store_true', help="Sort nuisances by impact")
    parser.add_argument("--pulls", action='store_true', help="Sort nuisances by impact")
    parser.add_argument("inputFile", 
        default="fitresults_123456789.root", 
        help="fitresults output ROOT file from combinetf")
    return parser.parse_args()


def printImpacts(args,fitresult,poi='Wmass'):
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    unit = 'MeV' if poi and poi.startswith('mass') else 'n.u. %'
    
    filtered_labels = labels

    if args.nuisance:
        filtered_labels = list(filter(lambda x: re.match(args.nuisance, x), labels))
        if not filtered_labels:
            raise ValueError(f"Failed to match any nuisances to expression '{args.nuisance}'")
        for l in filtered_labels:
            print(f"Impact of nuisance {l} is {impacts[list(labels).index(l)]*100} {unit}")
    else:
        print(f"Impact of all systematics (in {unit})")
        print("\n".join([f"   {k}: {round(v*100, 2)}" for k,v in zip(labels, impacts)]))

    if args.pulls:
        if not args.ungroup:
            raise ValueError("Pulls can only be printed for ungrouped impacts")
        pulls,constraints,_ = combinetf_input.get_pulls_and_constraints(args.inputFile, filtered_labels)
        print("-"*80)
        print("Pulls")
        print("-"*80)
        for l,p,c in zip(filtered_labels,pulls,constraints):
            print(f"{l}: {p:0.3f} +/- {c:0.3f}")


if __name__ == '__main__':
    args = parseArgs()
    fitresult = combinetf_input.get_fitresult(args.inputFile)
    for poi in combinetf_input.get_poi_names(fitresult):
        printImpacts(args,fitresult,poi)

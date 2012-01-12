from optparse import OptionParser

########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel",
		help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-p", "--path",
		help="Path to data",
		metavar="PATH", default="./",dest="path")
parser.add_option("-f", "--files",
		help="data files in quotes!",
		metavar="FILES", default="", dest="dfiles")

(options, args) = parser.parse_args()


listOFfiles=options.dfiles.split(":")

print listOFfiles

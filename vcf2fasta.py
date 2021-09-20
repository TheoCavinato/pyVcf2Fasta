from pysam import VariantFile
import argparse, os, shutil

"""
This code take alleles from a vcf file, translate the combination of alleles with IUPAC and write a fasta file
Use:
python3 vcf2fasta.py -vcf <fasta_file.fasta> -s <windows_size>
"""
#########################
#	USER ARGUMENTS	#
#########################

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-vcf', '--vcf_file', required=True, help="input vcf file to convert")
parser.add_argument('-o', '--output_prefix', required=True, help="prefix added to each output file")
parser.add_argument('-w', '--windows_size', type=int, required=False, help="in SNPs")
parser.add_argument('-s', '--sample_id', required=False, help="sample names to include (comma separated)")
args = parser.parse_args()

prefix=args.output_prefix
out_dir=prefix+"_fasta_from_vcf"
if os.path.isdir(out_dir):
	shutil.rmtree(out_dir)
os.mkdir(out_dir)
vcf_file=VariantFile(args.vcf_file)

#	IUPAC code to convert base couples into one unique letter
IUPAC={"AT":'W', "CG":'S', "AC":'M', "GT":'K', "AG":'R', "CT":'Y', "AA":'A', "TT":'T', "GG":'G', "CC":'C', "TA":'W', "GC":'S', "CA":'M', "TG":'K', "GA":'R', "TC":'Y'}

#	Get Individuals
if not args.sample_id:
	samples=list(vcf_file.header.samples)

if args.sample_id:
	samples_id=args.sample_id
	samples=samples_id.split(",")
	print("---Following samples will be considered:"+str(samples))

#################################################
#	WRITE THE ALLELES TO IUPAC FORMAT	#
#################################################

prefix_list=[]
if not args.windows_size:
	file_writer=[open(out_dir+"/"+prefix+"_"+sample+".fasta", 'w') for sample in samples]
	prefix_list.append(prefix)
	for nsample, file in enumerate(file_writer):
				file.write('\n')
				file.write(">"+samples[nsample])
				file.write('\n')

out_dir_meta=prefix+"_meta_data"
meta_file=open(out_dir_meta, 'w')
meta_file.write("\t".join(["chr", "start", "end", "sites", "\n"]))
site_counter=0
for nrec, rec in enumerate(vcf_file):
	#	If we are in a new windows, create files for these windows
	if args.windows_size:
		site_counter+=1
		if not nrec%args.windows_size:
			if 'file_writer' in locals():
				for file in file_writer:
					file.close()

			new_prefix=prefix+'_w'+str(nrec)+'_'
			prefix_list.append(new_prefix)
			file_writer=[open(out_dir+"/"+new_prefix+sample+".fasta", 'w') for sample in samples]

			if 'last_pos' in locals():
				meta_file.write("\t")
				meta_file.write("\t".join([str(last_pos), str(site_counter)]))
				meta_file.write("\n")

			meta_file.write("\t".join([rec.chrom, str(rec.pos)]))

			for nsample, file in enumerate(file_writer):
				file.write('\n')
				file.write(">"+samples[nsample])
				file.write('\n')
			site_counter=0

		last_pos = rec.pos
	#	Convert alleles to IUPAC and then add them to the files
	for nsample, file in enumerate(file_writer):
		#dois-je prendre en compte les missing data??
		alleles=''.join(rec.samples[samples[nsample]].alleles)
		if alleles in IUPAC.keys():
			file.write(IUPAC[alleles])
		else:
			raise Exception("The pair of allele is not convertible into IUPAC code")

meta_file.write("\t")
meta_file.write("\t".join([str(last_pos), str(site_counter)]))
meta_file.close()

#close the files
for file in file_writer:
	file.close()


print("---Segment files written")

#########################################
#	CONCATENATE FILES TOGETHER	#
#########################################

path=out_dir+"/"
for prefix in prefix_list:
	files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and prefix in i]
	output_file=open(prefix+"final.fasta", 'w')
	for file in files:
		for line in open(path+file, 'r'):
			output_file.write(line)
	output_file.close()
	
print("---Segment files concatenated")

shutil.rmtree(out_dir+"/")
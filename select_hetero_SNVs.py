### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
				python select_hetero_SNVs.py
				--in <VCF(INPUT)>
				--out <VCF(OUTPUT)>
				"""

import sys
import matplotlib.pyplot as plt

# --- end of imports --- #

def main( arguments ):
	"""! @brief runs everything """
	
	input_vcf = arguments[ arguments.index('--in')+1 ]
	output_vcf = arguments[ arguments.index('--out')+1 ]
	
	counter = 0
	error_counter = 0
	cov_filtered = 0
	allele_freq_filtered = 0
	with open( output_vcf, "w" ) as out:
		with open(  input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if parts[6] == "PASS":
						try:
							x, y = map( float, parts[-1].split(':')[1].split(',') )
							if 0.1 < y / (x+y) < 0.9:
								if 50 < x+y < 500:
									out.write( line )
									counter += 1
								else:
									cov_filtered += 1
							else:
								allele_freq_filtered += 1
						except:
							error_counter += 1 
				else:
					out.write( line )
				line = f.readline()
	print "number of remaining variants: " + str( counter )
	print "number of failed variants: " + str( error_counter )
	print "cov filtered: " + str( cov_filtered )
	print "allele freq filtered: " + str( allele_freq_filtered )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

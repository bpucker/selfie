### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python select_candidates.py
					--anno <FUNCTIONAL_ANNOTATION_FILE>
					--gff <GFF3_FILE>
					--in <HIGH_IMPACT_SNPEFF_FILE>
					--out <OUTPUT_FILE>
					--start <START_POSITION>
					--end <END_POSITION>
					"""

import re, sys

# --- end of imports --- #

def load_genes_in_region( gff_file, start_cutoff, end_cutoff ):
	"""! @brief load all genes in region specified by start and end """
	
	genes = []
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[ 2 ] == "gene":
				start, end = map( int, parts[ 3:5 ] )
				if end > start:
					if start < end:
						genes.append( { 'start': start, 'end': end, 'id': re.findall( "g\d+", parts[-1] )[0] } )
			line = f.readline()
	return genes


def load_functional_anno( functional_annotation_file ):
	""""! @brief load functional annotation """
	
	anno = {}
	with open( functional_annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno.update( { parts[0].split('.')[0]: parts[-1] } )
			line = f.readline()
	return anno


def main( arguments ):
	"""! @brief run everything """
	
	functional_annotation_file = arguments[ arguments.index('--anno')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	input_variant_file = arguments[ arguments.index('--in')+1 ]
	output_variant_file = arguments[ arguments.index('--out')+1 ]
	
	start_cutoff = int( arguments[ arguments.index('--start')+1 ] )
	end_cutoff = int( arguments[ arguments.index('--end')+1 ] )


	#genes = load_genes_in_region( gff_file, start_cutoff, end_cutoff )
	annotation = load_functional_anno( functional_annotation_file )
	
	with open( output_variant_file, "w" ) as out:
		with open( input_variant_file, "r" ) as f:
			out.write( f.readline() )
			line = f.readline()
			while line:
				parts = line.strip().split('\t')[:-1]
				pos = int( parts[1] )
				if start_cutoff < pos < end_cutoff:	#variant is in region of interest
					try:
						parts.append( annotation[ parts[4] ] )
					except KeyError:
						parts.append( "n/a" )
					out.write( "\t".join( parts ) + '\n' )
				line = f.readline()


if '--anno' in sys.argv and '--gff' in sys.argv and '--in' in sys.argv and '--out' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

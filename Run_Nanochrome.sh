#!/bin/bash

CHECK=1
if hash python3 2>/dev/null; then 
	echo $'\n'"Found python3" 
else
	echo $'\n'"python3 doesn't appear to be in PATH!" 1>&2;
	CHECK=0
fi

if hash minimap2 2>/dev/null; then 
	echo $'\n'"Found minimap2" 
else
	echo $'\n'"minimap2 doesn't appear to be in PATH!" 1>&2;
	CHECK=0
fi

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

usage() { echo "
	Usage: $0
	-g <file> Genome to scaffold in fasta format
	-r <file(,file)> Chromium paired reads (if one file will assume interleaved)
	-n <file> Nanopore reads, can be fasta, or fastq
	-f <integer> 10x Library mean fragment length
	-p <string> Prefix for output/processing files
	
	OPTIONAL:
	-l <integer> 3 - 20, Leniency: error/contiguity trade-off. Higher = More Error [5]
	-c <flag> (clean up) Add flag if you wish to clean up non-essential processing files
	" 1>&2; exit 1; }

l=5
c=0

while getopts ":g:r:n:f:p:l:c" o; do
    case "${o}" in
        f)
            f=${OPTARG}
            if [[ -n ${f//[0-9]/} ]]; then
			echo $'\n'"-f fragment length must be an integer (probably 20000-70000)"
			usage
	    fi 
            ;;
		l)
            l=${OPTARG}
            if [[ -n ${l//[0-9]/} ]]; then
				echo $'\n'"-l leniency must be an integer between 3 and 20 (inclusive)"
				usage
			fi
			if $l < 3 || $l > 20; then
				echo $'\n'"-l leniency must be an integer between 3 and 20 (inclusive)"
				usage
			fi
            ;;
        g)
            g=${OPTARG}
            ;;
		r)
            rd=${OPTARG}
			rd=${rd//,/ }
            ;;
		n)
            n=${OPTARG}
            ;;
		p)
            p=${OPTARG}
            ;;
		c)
			c=1
			;;
        *)
	    echo $'\n'"Unrecognised entry: ${OPTARG}"
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${g}" ] || [ -z "${rd}" ] || [ -z "${n}" ] || [ -z "${f}" ] || [ -z "${p}" ]; then
    usage
fi

if [ "$CHECK" -eq "0" ]; then
	exit 1
fi

echo $'\n'"Genome is: 		${g}"
echo "10x Reads are: 		${rd}"
echo "Nanopores are: 		${n}"
echo "Output prefix is: 	${p}"
echo "Fragment length is: 	${f}"
echo "Leniency is set to:   	${l}"

echo $'\n'"Running minimap2 - mapping 10x reads to genome.."
echo $'\n'"CMD: minimap2 -x sr ${g} ${rd} > ${p}_short.paf"

minimap2 -x sr ${g} ${rd} > ${p}_short.paf

echo $'\n'"Running chrome_candidates.py to generate candidate physical links"
echo $'\n'"CMD: python3 $DIR/chrome_candidates.py ${g} ${p}_short.paf ${f} ${p}"

python3 $DIR/chrome_candidates.py ${g} ${p}_short.paf ${f} ${p}

echo $'\n'"Running minimap2 - mapping nanopore reads to candidate scaffolds"
echo $'\n'"CMD: minimap2 -x ${p}_candidates.fa ${n} > ${p}_aln.paf"

minimap2 -x map-ont ${p}_candidates.fa ${n} > ${p}_aln.paf

echo $'\n'"Running nano_confirms.py to put everything together!"
echo $'\n'"CMD: python3 $DIR/nano_confirms.py ${g} ${p}_table.tsv ${p}_aln.paf ${f} ${p}"

python3 $DIR/nano_confirms.py ${g} ${p}_table.tsv ${p}_aln.paf ${f} ${l} ${p}

if [ "$c" -eq "1" ]; then
	rm -f ${p}_aln.sam ${p}_aln.paf
	rm -f ${p}_table.tsv
	rm -f ${p}_candidates.fa
fi


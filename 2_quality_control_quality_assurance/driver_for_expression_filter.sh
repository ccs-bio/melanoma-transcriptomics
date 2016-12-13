#!/bin/sh

BASE_DIR='/path/to/base/analysis/directory'

SCRIPT_PATH=$BASE_DIR'/analysis_scripts_103'

# QUANT_TYPES=( 'counts' 'fpkm' )
# QUANT_TYPES=( 'counts_exonic' )
QUANT_TYPES=( 'fpkm_cuffdiff' )

ALN_TYPES=( 'unique' )


echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	"

# -------------------------------------------- Contaminating Biotype Removal ------------------------------------------

CONTAMINANT_BIOTYPE_LIST='/path/to/gene_ids/that/will/be/removed/1_gene_ids_for_Contaminating_Biotypes.txt'




# -------------------------------------------------- Expression Filter ------------------------------------------------

for QUANT_TYPE in "${QUANT_TYPES[@]}"
{
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	"
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Processing: "$QUANT_TYPE

	for ALN_TYPE in "${ALN_TYPES[@]}"
	{
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	-"$QUANT_TYPE', '$ALN_TYPE



		RESULTS_BASE_DIR=$BASE_DIR'/results-103/'$QUANT_TYPE'/'$ALN_TYPE

		INPUT_DIR=$RESULTS_BASE_DIR'/All'

		OUTPUT_DIR=$RESULTS_BASE_DIR'/Filtered'

		mkdir -p $OUTPUT_DIR


		# --------------------------- Polluted File -------------------------------

		RESULTS_FILE_NAME='aggregated_expression_profile_for_'$QUANT_TYPE'_AllGenes-Polluted.txt'
		RESULTS_FILE=$INPUT_DIR'/'$RESULTS_FILE_NAME



		# --------------------------- Biotype Removal -----------------------------

		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	Removing Contaminants Biotypes..."
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

		RESULTS_FILE_BIOTYPE_FILTERED=$INPUT_DIR'/aggregated_expression_profile_for_'$QUANT_TYPE'_AllGenes.txt'

		# fgrep -v -F -f $CONTAMINANT_BIOTYPE_LIST $RESULTS_FILE > $RESULTS_FILE_BIOTYPE_FILTERED



		# -------------------------------- Params ---------------------------------

		EXP_INTERPRETATION="NONE,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,C,C"

		TYPE_FLAG='--'$QUANT_TYPE
		
		if [ "$QUANT_TYPE" == "counts_exonic" ]; then
			TYPE_FLAG='--counts'
		fi
		
		if [ "$QUANT_TYPE" == "fpkm_cuffdiff" ]; then
			TYPE_FLAG='--fpkm'
		fi
		
		# --------------------------- Filter Stringency ---------------------------
		
		PERCENT_PRESENCE="10"

		MIN_EXPRESSION_VALUE="0"	# Minimum expression value for a gene to pass

		if [ "$QUANT_TYPE" == "counts" ] || [ "$QUANT_TYPE" == "counts_exonic" ]; then
			MIN_EXPRESSION_VALUE="800"
			PERCENT_PRESENCE="33"
		fi

		# The FPKM is considerable lower as the average FPKM for the Pseudogene biotype is very low, 0.097.
		if [ "$QUANT_TYPE" == "fpkm" ] || [ "$QUANT_TYPE" == "fpkm_cuffdiff" ]; then
			MIN_EXPRESSION_VALUE="1.0"
			PERCENT_PRESENCE="10"
		fi

		# ---------------------------------- Run ----------------------------------

		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	Filtering..."
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

		$SCRIPT_PATH/script_expression_filter.py 	-i $RESULTS_FILE_BIOTYPE_FILTERED \
													-a $QUANT_TYPE \
													-o $OUTPUT_DIR \
													--interpretation $EXP_INTERPRETATION \
													--min_exp_val $MIN_EXPRESSION_VALUE \
													--presence $PERCENT_PRESENCE \
													$TYPE_FLAG
	}
}

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]	"



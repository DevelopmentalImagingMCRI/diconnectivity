#!/bin/bash

. `which DIConnEnv.sh`

if [ -z "$1" ]
then
	echo "Usage: $0 <subject list>"
else
	SUBJECTLIST=$1

	if [ ! -f "$SUBJECTLIST" ]
	then
		echo "The subject list given does not exist"
		exit 1;
	else
#rm -fr $DWIREGTOSTRUCTQC
		mkdir -p $DWIREGTOSTRUCTQC
		for i in `cat $SUBJECTLIST`
		do
			echo $i
#slices_output $DWIREGTOSTRUCTQC/${i}_linear.png ANTSSkullStrip2mm/${i}BrainExtractionBrain $DWIREGTOSTRUCT/$i/bzero_brain_linear_reg
#slices_output $DWIREGTOSTRUCTQC/${i}_nonlinear.png ANTSSkullStrip2mm/${i}BrainExtractionBrain $DWIREGTOSTRUCT/$i/bzero_brain_ants_reg
			slices_output $DWIREGTOSTRUCTQC/${i}_linear.png $T1SKULLSTRIPPED/${i} $DWIREGTOSTRUCT/$i/bzero_brain_linear_reg &
			slices_output $DWIREGTOSTRUCTQC/${i}_nonlinear.png $T1SKULLSTRIPPED/${i} $DWIREGTOSTRUCT/$i/bzero_brain_ants_reg &
#slices_output $DWIREGTOSTRUCTQC/${i}_nonlinear_wm.png $T1SKULLSTRIPPED/${i} $DWIREGTOSTRUCT/$i/bzero_brain_ants_reg_wm &
#slices_output $DWIREGTOSTRUCTQC/${i}_nonlinear_fa.png $T1SKULLSTRIPPED/${i} $DWIREGTOSTRUCT/$i/bzero_brain_ants_reg_fa &
			wait;
#slices_output $DWIREGTOSTRUCTQC/${i}_nonlinear_combo.png $T1SKULLSTRIPPED/${i} $DWIREGTOSTRUCT/$i/bzero_brain_ants_reg_combo

		done
	fi
fi	

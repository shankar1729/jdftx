#!/bin/bash

echo 13 #expected lines of output

#H2 geometry:
awk 'BEGIN { nIons=0 }
$1=="ion" && $2=="H" {
	nIons++;
	if(nIons==1) { x=$3; y=$4; z=$5 }
	if(nIons==2)
	{	dx = $3-x;
		dy = $4-y;
		dz = $5-z;
		r = 0.5291772 * sqrt(dx*dx + dy*dy + dz*dz);
		print r, "0.7583 0.010 H2 bond length [A]"
	}
}
' H2.ionpos

#H2 vibrations:
awk 'BEGIN { nFreq=0; nIR=0; }
	/Frequency:/ { nFreq++; printf("%f %f 100 H2 vib frequency %d [inv-cm]\n", $5, nFreq==6 ? 4376 : 0, nFreq); }
	/IR intensity:/ { nIR++; printf("%f 0 10 H2 IR intensity %d [km/mol]\n", $6, nIR); }
' H2_vibrations.out

#!/bin/bash

xcHeader="$1"

echo "#ifndef JDFTX_XCMAP_H"
echo "#define JDFTX_XCMAP_H"

awk '$1=="#define" && $2~/XC_/ {
		id = $2;
		#Construct names:
		name = tolower(id);
		gsub("xc_", "", name);
		gsub("_", "-", name);
		names[id] = name;
		#Extract type:
		split(id, tokens, "_");
		if(id ~ /HYB/)
			type = tokens[4];
		else
			type = tokens[3];
		types[id] = type;
	}
	
	function printType(type)
	{
		iDone = 0;
		printf("\nEnumStringMap<int> xcMap_%s(", type);
		for(id in names)
			if(types[id] == type)
			{	if(iDone) printf(",");
				printf("\n\t%s, \"%s\"", id, names[id]);
				iDone++;
			}
		printf("\n);\n");
	}
	
	END {
		printType("X");
		printType("C");
		printType("XC");
		printType("K");
	}' $xcHeader

echo "#endif // JDFTX_XCMAP_H"

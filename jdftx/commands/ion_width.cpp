/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <commands/command.h>
#include <electronic/Everything.h>

struct CommandIonWidth : public Command
{
	CommandIonWidth() : Command("ion-width")
	{
		format = "Ecut | fftbox | <width>";
		comments = "Manually specify width of gaussian representations of nuclear charge in bohr\n"
			"or set automatically based on either energy cut-off (Ecut) or grid spacing (fftbox)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	string key; pl.get(key, string(), "width");
		if(!key.length() || key=="Ecut") e.iInfo.ionWidthMethod = IonInfo::IonWidthEcut;
		else if(key=="fftbox") e.iInfo.ionWidthMethod = IonInfo::IonWidthFFTbox;
		else
		{	istringstream iss(key);
			iss >> e.iInfo.ionChargeWidth;
			if(iss.fail()) throw string("<width> must be Ecut, fftbox or a value in bohrs");
			e.iInfo.ionWidthMethod = IonInfo::IonWidthManual;
		}
	}

	void printStatus(Everything& e, int iRep)
	{	switch(e.iInfo.ionWidthMethod)
		{	case IonInfo::IonWidthFFTbox: logPrintf("fftbox"); break;
			case IonInfo::IonWidthEcut: logPrintf("Ecut"); break;
			case IonInfo::IonWidthManual: logPrintf("%lg", e.cntrl.Ecut); break;
		}
	}
}
commandIonWidth;

#include <electronic/JDFTxWrapper.h>
#include <electronic/Everything.h>
#include <electronic/IonicMinimizer.h>
#include <electronic/LatticeMinimizer.h>
#include <commands/parser.h>


JDFTxWrapper::JDFTxWrapper(std::vector<std::pair<string, string>> inputs, bool variableCell)
{	e = std::make_shared<Everything>();
	parse(inputs, *e, true);
	e->setup();
	e->dump(DumpFreq_Init, 0);
	Citations::print();
	logPrintf("Initialization completed successfully at t[s]: %9.2lf\n\n", clock_sec());
	logFlush();
	
	if(variableCell)
	{	LatticeGradient grad; grad.init(e->iInfo);
		lmin = std::make_shared<LatticeMinimizer>(*e);
		lmin->compute(&grad, NULL);
	}
	else
	{	IonicGradient grad; grad.init(e->iInfo);
		imin = std::make_shared<IonicMinimizer>(*e);
		imin->compute(&grad, NULL);
	}
}

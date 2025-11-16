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
		lmin = std::make_shared<LatticeMinimizer>(*e);
	else
		imin = std::make_shared<IonicMinimizer>(*e);
}

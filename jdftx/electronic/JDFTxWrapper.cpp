#include <electronic/JDFTxWrapper.h>
#include <electronic/Everything.h>
#include <electronic/IonicMinimizer.h>
#include <electronic/LatticeMinimizer.h>
#include <commands/parser.h>


JDFTxWrapper::JDFTxWrapper(std::vector<std::pair<string, string>> inputs, bool variableCell)
: variableCell(variableCell)
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
		lmin->report(-1);
	}
	else
	{	IonicGradient grad; grad.init(e->iInfo);
		imin = std::make_shared<IonicMinimizer>(*e);
		imin->compute(&grad, NULL);
		imin->report(-1);
	}
}

void JDFTxWrapper::minimize()
{	if(variableCell)
		lmin->minimize(e->latticeMinParams);
	else
		imin->minimize(e->ionicMinParams);
}

void convertArray(const NDarray& arr, IonicGradient& grad, std::string name)
{	if(arr.shape.size() != 2) throw std::invalid_argument(name + " must be 2D");
	if(arr.shape[1] != 3) throw std::invalid_argument(name + ".shape[1] must be 3");
	std::vector<size_t> index({0, 0});
	size_t& i_atom = index[0];
	size_t& i_dir = index[1];
	for(std::vector<vector3<>>& grad_sp: grad)
		for(vector3<>& grad_sp_a: grad_sp)
		{	for(i_dir=0; i_dir<3; i_dir++)
				grad_sp_a[i_dir] = arr[index];
			i_atom++;
		}
	if(arr.shape[0] != i_atom)
	{	std::ostringstream err; err << "Expected " << name << ".shape[1] == " << i_atom;
		throw std::invalid_argument(err.str());
	}
}

void convertArray(const NDarray& arr, matrix3<>& m, std::string name)
{	if(arr.shape.size() != 2) throw std::invalid_argument(name + " must be 2D");
	if(arr.shape[0] != 3) throw std::invalid_argument(name + ".shape[0] must be 3");
	if(arr.shape[1] != 3) throw std::invalid_argument(name + ".shape[1] must be 3");
	std::vector<size_t> index({0, 0});
	size_t& i_row = index[0];
	size_t& i_col = index[1];
	for(i_row=0; i_row<3; i_row++)
		for(i_col=0; i_col<3; i_col++)
			m(i_row, i_col) = arr[index];
}

void JDFTxWrapper::move(NDarray delta_positions, NDarray delta_R)
{	LatticeGradient delta, grad;
	delta.init(e->iInfo);
	grad.init(e->iInfo);
	convertArray(delta_positions, delta.ionic, "delta_positions");
	convertArray(delta_R, delta.lattice, "delta_R");
	delta.ionic = e->gInfo.R * delta.ionic; //input fractional, step expects Cartesian
	if(variableCell)
	{	delta.lattice = delta.lattice * e->gInfo.invR; //convert delta(R) to relative Cartesian strain
		lmin->step(delta, 1.0);
		lmin->compute(&grad, NULL);
		lmin->report(-1);
	}
	else
	{	if(nrm2(delta.lattice) > 1E-12)
			throw std::invalid_argument("Cell change not allowed (variableCell = false)");
		imin->step(delta.ionic, 1.0);
		imin->compute(&grad.ionic, NULL);
		imin->report(-1);
	}
}

#include <electronic/JDFTxWrapper.h>
#include <electronic/Everything.h>
#include <electronic/IonicMinimizer.h>
#include <electronic/LatticeMinimizer.h>
#include <commands/command.h>
#include <commands/parser.h>
#include <core/LatticeUtils.h>
#include <fluid/PCM.h>

//Wrap cavity function to provide ScalarField and ScalarFieldArray to NDarray conversions:
struct CavityFunctionWrapper
{	setCavity_t setCavity;
	CavityFunctionWrapper(setCavity_t setCavity) : setCavity(setCavity) {}
	void operator()(const ScalarField& n, ScalarFieldArray& shape); //implemented below
};
struct CavityFunctionGradWrapper
{	setCavity_t setCavityGrad;
	CavityFunctionGradWrapper(setCavity_t setCavityGrad) : setCavityGrad(setCavityGrad) {}
	void operator()(const ScalarFieldArray& A_shape, ScalarField& A_n); //implemented below
};


JDFTxWrapper::JDFTxWrapper(std::vector<std::pair<string, string>> inputs, bool variableCell, setCavity_t setCavity[2])
: variableCell(variableCell)
{	logPrintf("\n------------------------ Initializing JDFTxWrapper -----------------------\n");
	e = std::make_shared<Everything>();
	parse(inputs, *e, true);
	if(setCavity[0])
	{	//Set cavity functions if any:
		FluidSolverParams& fsp = e->eVars.fluidParams;
		fsp.cavityFunction = CavityFunctionWrapper(setCavity[0]);
		fsp.cavityFunctionGrad = CavityFunctionGradWrapper(setCavity[1]);
	}
	e->setup();
	if(variableCell) e->iInfo.computeStress = true;
	e->dump(DumpFreq_Init, 0);
	Citations::print();
	logPrintf("Initialization completed successfully at t[s]: %9.2lf\n\n", clock_sec());
	logFlush();
	
	if(variableCell)
		lmin = std::make_shared<LatticeMinimizer>(*e);
	else
		imin = std::make_shared<IonicMinimizer>(*e);
	compute();
}

void JDFTxWrapper::minimize()
{	if(variableCell)
		lmin->minimize(e->latticeMinParams);
	else
		imin->minimize(e->ionicMinParams);
}

void convertArray(NDarray arr, IonicGradient& grad, std::string name)
{	if(arr.shape.size() != 2) throw std::invalid_argument(name + " must be 2D");
	if(arr.shape[1] != 3) throw std::invalid_argument(name + ".shape[1] must be 3");
	assert(not arr.onGpu);
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

void convertArray(NDarray arr, matrix3<>& m, std::string name)
{	if(arr.shape.size() != 2) throw std::invalid_argument(name + " must be 2D");
	if(arr.shape[0] != 3) throw std::invalid_argument(name + ".shape[0] must be 3");
	if(arr.shape[1] != 3) throw std::invalid_argument(name + ".shape[1] must be 3");
	assert(not arr.onGpu);
	std::vector<size_t> index({0, 0});
	size_t& i_row = index[0];
	size_t& i_col = index[1];
	for(i_row=0; i_row<3; i_row++)
		for(i_col=0; i_col<3; i_col++)
			m(i_row, i_col) = arr[index];
}

void JDFTxWrapper::move(NDarray delta_positions, NDarray delta_R)
{	LatticeGradient delta; delta.init(e->iInfo);
	convertArray(delta_positions, delta.ionic, "delta_positions");
	convertArray(delta_R, delta.lattice, "delta_R");
	delta.ionic = e->gInfo.R * delta.ionic; //input fractional, step expects Cartesian
	delta.lattice = delta.lattice * e->gInfo.invR; //convert delta(R) to relative Cartesian strain
	
	//Check symmetries of the new lattice:
	matrix3<> Rnew = e->gInfo.R + delta.lattice * e->gInfo.R;
	matrix3<> metric = (~Rnew) * Rnew;
	for(const SpaceGroupOp& op: e->symm.getMatrices())
		if(nrm2(metric - (~op.rot) * metric * op.rot) > symmThreshold * nrm2(metric))
		{	logPrintf("\nNew lattice breaks previous symmetries\n");
			throw std::invalid_argument("Symmetry change");
		}
	
	//Check symmetries of the change in positions:
	IonicGradient delta_ionic_err = clone(delta.ionic);
	IonicGradient force_dir = e->gInfo.RT * delta.ionic; //convert to contravariant lattice coords (force-like)
	e->symm.symmetrize(force_dir); //... becasue symmetrize expects forces in lattice coords
	delta.ionic = e->gInfo.invRT * force_dir; //save symmetrized version back in Cartesian coords
	delta_ionic_err -= delta.ionic;
	if(dot(delta_ionic_err, delta_ionic_err) > symmThresholdSq)
	{	logPrintf("\nNew positions break previous symmetries\n");
		throw std::invalid_argument("Symmetry change");
	}
	//Step to new lattice/positions:
	if(variableCell)
		lmin->step(delta, 1.0);
	else
	{	if(nrm2(delta.lattice) > 1E-12)
		{	logPrintf("\nCell change not allowed (variableCell = false)\n");
			throw std::invalid_argument("Cell change not allowed (variableCell = false)");
		}
		imin->step(delta.ionic, 1.0);
	}
	compute();
}

void JDFTxWrapper::dumpEnd() const
{	e->dump(DumpFreq_End, 0);
}

double JDFTxWrapper::getEnergy() const
{	return relevantFreeEnergy(*e);
}

NDarray JDFTxWrapper::getForces() const
{	NDarray result;
	result.data = (double*)&forces[0][0];
	result.shape = std::vector<size_t>({forces.size(), 3});
	result.onGpu = false;
	return result;
}

NDarray JDFTxWrapper::getStress() const
{	NDarray result;
	result.data = (double*)&e->iInfo.stress(0, 0);
	result.shape = std::vector<size_t>({3, 3});
	result.onGpu = false;
	return result;
}

std::vector<string> JDFTxWrapper::getCommands()
{	std::vector<string> commands;
	for(const auto& mapEntry: getCommandMap())
		commands.push_back(mapEntry.first);
	return commands;
}

string JDFTxWrapper::getCommandDoc(string name)
{	const Command& c = *getCommandMap()[name];
	return c.name + ' ' + c.format + '\n' + c.comments + '\n';
}

void JDFTxWrapper::compute()
{	LatticeGradient grad; grad.init(e->iInfo);
	if(variableCell)
	{	lmin->compute(&grad, NULL);
		lmin->report(-1);
	}
	else
	{	imin->compute(&grad.ionic, NULL);
		imin->report(-1);
	}
	
	//Save forces in flat array:
	forces.clear();
	for(const std::vector<vector3<>>& grad_sp: grad.ionic)
		for(const vector3<>& grad_sp_a: grad_sp)
			forces.push_back(-grad_sp_a); //force is negative gradient (grad is already Cartesian)
	
	//Write periodicity (needed for recreating geometry in restart):
	if(e->dump.count(std::make_pair(DumpFreq_Ionic, DumpIonicPositions)) and mpiWorld->isHead())
	{	string fname = e->dump.getFilename("pbc");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		auto isTruncated = e->coulombParams.isTruncated();
		for(int i_dir=0; i_dir<3; i_dir++)
			fprintf(fp, "%d ", isTruncated[i_dir] ? 0 : 1);
		fclose(fp);
		logPrintf("done\n"); logFlush();
	}
}

void CavityFunctionWrapper::operator()(const ScalarField& n, ScalarFieldArray& shape)
{	vector3<size_t> S(n->gInfo.S);
	nullToZero(shape, n->gInfo);
	std::vector<size_t> size({S[0], S[1], S[2]}), strides; //strides empty => contiguous
	NDarray n_arr({n->dataPref(), size, strides, isGpuEnabled()});
	std::vector<NDarray> shape_arr;
	for(ScalarField& shape_i: shape)
		shape_arr.push_back(NDarray({shape_i->dataPref(), size, strides, isGpuEnabled()}));
	setCavity(n_arr, shape_arr);
}

void CavityFunctionGradWrapper::operator()(const ScalarFieldArray& A_shape, ScalarField& A_n)
{	vector3<size_t> S(A_shape[0]->gInfo.S);
	nullToZero(A_n, A_shape[0]->gInfo);
	std::vector<size_t> size({S[0], S[1], S[2]}), strides;  //strides empty => contiguous
	NDarray A_n_arr({A_n->dataPref(), size, strides, isGpuEnabled()});
	std::vector<NDarray> A_shape_arr;
	for(const ScalarField& A_shape_i: A_shape)
		A_shape_arr.push_back(NDarray({A_shape_i->dataPref(), size, strides, isGpuEnabled()}));
	setCavityGrad(A_n_arr, A_shape_arr);
}


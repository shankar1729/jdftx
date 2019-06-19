/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <core/ManagedMemory.h>
#include <core/GpuUtil.h>
#include <fftw3.h>
#include <mutex>
#include <map>
#include <set>

//-------- Memory usage profiler ---------

namespace MemUsageReport
{
	enum Mode { Add, Remove, Print };
	
	//Add, remove or retrieve memory report based on mode
	void manager(Mode mode, string category=string(), size_t nBytes=0)
	{	
		#ifdef ENABLE_PROFILING
		struct Usage
		{	size_t current, peak; //!< current and peak memory usage (in unit of complex numbers i.e. 16 bytes)
			Usage() : current(0), peak(0) {}
			
			Usage& operator+=(size_t n)
			{	current += n;
				if(current > peak)
					peak = current;
				return *this;
			}
			
			Usage& operator-=(size_t n)
			{	current -= n;
				return *this;
			}
		};
		static std::map<string, Usage> usageMap;
		static Usage usageTotal;
		static std::mutex usageLock;
		static const double bytesToGB = 1./pow(1024.,3);
		
		switch(mode)
		{	case Add:
			{	usageLock.lock();
				/*
				//DEBUG: Uncomment and tweak this block to trace execution point of highest memory of a given category
				//NOTE: Avoid commiting changes to this block during memory optimization;
				//restore it to this state and remember to comment it out before commiting!
				if(category=="ColumnBundle"
					&& (usageMap[category].current+nBytes > usageMap[category].peak) )
				{	printStack(true);
					logPrintf("MEMUSAGE: %30s %12.6lf GB\n", category.c_str(), (usageMap[category].current+nBytes)*bytesToGB);
				}
				*/
				usageMap[category] += nBytes;
				usageTotal += nBytes;
				usageLock.unlock();
				assert(category.length());
				break;
			}
			case Remove:
			{	usageLock.lock();
				usageMap[category] -= nBytes;
				usageTotal -= nBytes;
				usageLock.unlock();
				assert(category.length());
				break;
			}
			case Print:
			{	for(auto entry: usageMap)
					logPrintf("MEMUSAGE: %30s %12.6lf GB\n", entry.first.c_str(), entry.second.peak * bytesToGB);
				logPrintf("MEMUSAGE: %30s %12.6lf GB\n", "Total", usageTotal.peak * bytesToGB);
				break;
			}
		}
		#endif //ENABLE_PROFILING
	}
}


//-------- Memory pool to reduce system alloc/free calls ---------

namespace MemPool
{
	//Pool memory allocations in a memory space abstracted by MemSpace
	//MemSpace is a tag class with static functions:
	// void* alloc(size_t);  //returns 0 when out of memory
	// void free(void*);     //assumed to not fail
	// void outOfMemory();   //exit with appropriate out of memory error
	template<typename MemSpace> class MemPool
	{	uint8_t* pool; //pointer to entire pool of memory (allocated once)
		std::mutex lock; //for thread safety
		//Allocated memory
		std::map<size_t,size_t> used; //start -> stop
		//Available 'holes' in memory:
		std::map<size_t,size_t> holes; //start -> stop
		std::map<size_t,std::set<size_t>> holesBySize; //size -> set of starts
		//--- helper functions for managing holes
		typedef std::map<size_t,size_t>::iterator MapIter;
		typedef std::map<size_t,std::set<size_t>>::iterator MapSetIter;
		inline void printHoles()
		{	//List holes (for debugging only):
			logPrintf("Holes:");
			for(auto entry: holes)
				logPrintf(" (%lu,%lu)", entry.first, entry.second);
			logPrintf("\n");
		}
		inline void addHole(size_t start, size_t size)
		{	size_t startEff = start;
			size_t stopEff = start+size;
			MapIter ubound = holes.upper_bound(start); //iterator to hole after start
			MapIter lbound = ubound; if(lbound!=holes.begin()) lbound--; //iterator to hole before start
			//Check for contiguous hole before:
			if((lbound!=holes.end()) && (lbound->second==startEff))
			{	startEff = lbound->first; //absorb into new hole
				removeHole(lbound->first, &lbound); //remove old hole
			}
			//Check for contiguous hole after:
			if((ubound!=holes.end()) && (ubound->first==stopEff))
			{	stopEff = ubound->second; //absorb into new hole
				removeHole(ubound->first, &ubound); //remove old hole
			}
			//Add new hole:
			holes[startEff] = stopEff;
			holesBySize[stopEff-startEff].insert(startEff);
			//Uncomment following to debug:
			//logPrintf("Added (%lu,%lu) expanded to (%lu,%lu)\t",start,start+size, startEff,stopEff); printHoles();
		}
		inline void removeHole(size_t start, MapIter* holesPtr=0, MapSetIter* holesBySizePtr=0)
		{	//Remove from holes:
			MapIter holesIter = holesPtr ? *holesPtr : holes.find(start);
			assert((holesIter!=holes.end()) && (holesIter->first==start));
			size_t size = holesIter->second - start;
			holes.erase(holesIter);
			//Remove from holesBySize:
			MapSetIter holesBySizeIter = holesBySizePtr ? *holesBySizePtr : holesBySize.find(size);
			assert((holesBySizeIter!=holesBySize.end()) && (holesBySizeIter->first==size));
			holesBySizeIter->second.erase(start);
			if(!holesBySizeIter->second.size()) //no more holes of this size
				holesBySize.erase(holesBySizeIter);
			//Uncomment following to debug:
			//logPrintf("Deleted (%lu,%lu)\t", start,start+size); printHoles();
		}
	public:
		MemPool() : pool(0)
		{	if(mempoolSize)
			{	pool = (uint8_t*)MemSpace::alloc(mempoolSize);
				if(!pool) MemSpace::outOfMemory();
				addHole(0, mempoolSize);
			}
		}
		~MemPool()
		{	if(pool) MemSpace::free(pool);
		}
		void* alloc(size_t sizeRequested)
		{	if(!mempoolSize) return MemSpace::alloc(sizeRequested); //pool not in use
			lock.lock();
			//Find size adjusted to chunk size:
			const size_t chunkSize = 4096; //typical page size
			const size_t chunkMask = chunkSize - 1;
			size_t size = (sizeRequested + chunkMask) & (~chunkMask); //round up to multiple of chunkSize
			//Find hole just big enough to fit it:
			MapSetIter ubound = holesBySize.upper_bound(size);
			if(ubound == holesBySize.end())
			{	//No hole big enough left, so allocate externally:
				lock.unlock();
				void* ptr = MemSpace::alloc(sizeRequested);
				if(!ptr) MemSpace::outOfMemory();
				return ptr;
			}
			else
			{	//Hole found, so allocate from it:
				size_t start = *(ubound->second.begin());
				size_t holeSize = ubound->first;
				used[start] = start+size; //mark allocated range
				removeHole(start, 0, &ubound); //remove old hole
				if(holeSize > size) addHole(start+size, holeSize-size); //add hole left behind (if any)
				lock.unlock();
				return (void*)(pool+start);
			}
		}
		void free(void* ptr)
		{	if(!mempoolSize) return MemSpace::free(ptr); //pool not in use
			lock.lock();
			//Find in used map:
			size_t start = ((uint8_t*)ptr) - pool;
			MapIter usedIter = used.find(start);
			if(usedIter == used.end())
			{	//Not found in used => allocated externally
				MemSpace::free(ptr); //free externally
			}
			else
			{	//Found in used => allocated in pool
				size_t size = usedIter->second - start;
				used.erase(usedIter); //remove from used
				addHole(start, size); //add corresponding hole
			}
			lock.unlock();
		}
	};
	
	//---- MemSpace classes for each memory space ----
	#if defined(GPU_ENABLED) && defined(PINNED_HOST_MEMORY)
	struct MemSpaceCPU
	{	static void* alloc(size_t size)
		{	void* ptr;
			cudaError_t ret = cudaMallocHost(&ptr, size);
			return (ret==cudaSuccess) ? ptr : 0;
		}
		static void free(void* ptr) { cudaFreeHost(ptr); }
		static void outOfMemory() die_alone("Host memory allocation failed (out of pinned memory)\n");
	};
	#else
	struct MemSpaceCPU
	{	static void* alloc(size_t size) { return fftw_malloc(size); }
		static void free(void* ptr) { fftw_free(ptr); }
		static void outOfMemory() die_alone("Memory allocation failed (out of memory)\n");
	};
	#endif
	#ifdef GPU_ENABLED
	struct MemSpaceGPU
	{	static void* alloc(size_t size)
		{	assert(isGpuMine());
			void* ptr;
			cudaError_t ret = cudaMalloc(&ptr, size);
			return (ret==cudaSuccess) ? ptr : 0;
		}
		static void free(void* ptr)
		{	assert(isGpuMine());
			cudaFree(ptr);
		}
		static void outOfMemory() die_alone("GPU memory allocation failed (out of memory)\n");
	};
	#endif
	
	//Pool accessor functions (to avoid file-level static variables):
	MemPool<MemSpaceCPU>& CPU() { static MemPool<MemSpaceCPU> pool; return pool; }
	#ifdef GPU_ENABLED
	MemPool<MemSpaceGPU>& GPU() { static MemPool<MemSpaceGPU> pool; return pool; }
	#endif
}


//---------- class ManagedMemoryBase -----------

void ManagedMemoryBase::reportUsage()
{	MemUsageReport::manager(MemUsageReport::Print);
}

//Free memory
void ManagedMemoryBase::memFree()
{	if(!nBytes) return; //nothing to free
	if(onGpu)
	{
		#ifdef GPU_ENABLED
		MemPool::GPU().free(c);
		#else
		assert(!"onGpu=true without GPU_ENABLED"); //Should never get here!
		#endif
	}
	else MemPool::CPU().free(c);
	MemUsageReport::manager(MemUsageReport::Remove, category, nBytes);
	onGpu = false;
	c = 0;
	nBytes = 0;
	category.clear();
}

//Allocate memory
void ManagedMemoryBase::memInit(string category, size_t nBytes, bool onGpu)
{	if(category==this->category && nBytes==this->nBytes && onGpu==this->onGpu) return; //already in required state
	memFree();
	this->category = category;
	this->nBytes = nBytes;
	this->onGpu = onGpu;
	if(onGpu)
	{
		#ifdef GPU_ENABLED
		c = MemPool::GPU().alloc(nBytes);
		#else
		assert(!"onGpu=true without GPU_ENABLED");
		#endif
	}
	else c = MemPool::CPU().alloc(nBytes);
	MemUsageReport::manager(MemUsageReport::Add, category, nBytes);
}

void ManagedMemoryBase::memMove(ManagedMemoryBase&& mOther)
{	std::swap(category, mOther.category);
	std::swap(nBytes, mOther.nBytes);
	std::swap(onGpu, mOther.onGpu);
	std::swap(c, mOther.c);
	//Now mOther will be empty, while *this will have all its contents
}

//Move data to CPU
void ManagedMemoryBase::toCpu() const
{	if(!onGpu || !c) return; //already on cpu, or no data
#ifdef GPU_ENABLED
	assert(isGpuMine());
	ManagedMemoryBase& me = *((ManagedMemoryBase*)this);
	void* cCpu = MemPool::CPU().alloc(nBytes);
	cudaMemcpy(cCpu, me.c, nBytes, cudaMemcpyDeviceToHost);
	MemPool::GPU().free(me.c); //Free GPU mem
	me.c = cCpu; //Make c a cpu pointer
	me.onGpu = false;
#endif
}

// Move data to GPU
void ManagedMemoryBase::toGpu() const
{	if(onGpu || !c) return; //already on gpu, or no data
#ifdef GPU_ENABLED
	assert(isGpuMine());
	ManagedMemoryBase& me = *((ManagedMemoryBase*)this);
	void* cGpu = MemPool::GPU().alloc(nBytes);
	cudaMemcpy(cGpu, me.c, nBytes, cudaMemcpyHostToDevice);
	MemPool::CPU().free(me.c); //Free CPU mem
	me.c = cGpu; //Make c a gpu pointer
	me.onGpu = true;
#else
	assert(!"toGpu() called without GPU_ENABLED");
#endif
}

//--------- ManagedMemory<complex> and ManagedMemory<double> operators ------

void scale(double alpha, ManagedMemory<double>& y)
{	callPref(eblas_dscal)(y.nData(), alpha, y.dataPref(), 1);
}
void scale(double alpha, ManagedMemory<complex>& y)
{	callPref(eblas_zdscal)(y.nData(), alpha, y.dataPref(), 1);
}
void scale(complex alpha, ManagedMemory<complex>& y)
{	callPref(eblas_zscal)(y.nData(), alpha, y.dataPref(), 1);
}
void scale(const ManagedMemory<complex>& x, ManagedMemory<complex>& y)
{	assert(x.nData() == y.nData());
	callPref(eblas_zmul)(x.nData(), x.dataPref(), 1, y.dataPref(), 1);
}

void axpy(complex alpha, const ManagedMemory<complex>& x, ManagedMemory<complex>& y)
{	assert(x.nData() == y.nData());
	callPref(eblas_zaxpy)(x.nData(), alpha, x.dataPref(), 1, y.dataPref(), 1);
}

double nrm2(const ManagedMemory<complex>& a)
{	return callPref(eblas_dznrm2)(a.nData(), a.dataPref(), 1);
}

complex dotc(const ManagedMemory<complex>& a, const ManagedMemory<complex>& b)
{	assert(a.nData() == b.nData());
	return callPref(eblas_zdotc)(a.nData(), a.dataPref(), 1, b.dataPref(), 1);
}

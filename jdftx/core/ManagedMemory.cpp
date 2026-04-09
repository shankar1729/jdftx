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

#ifdef GPU_ENABLED
extern bool prefetchSupported; //defined in GpuUtil.cpp
extern int memLocDevice, memLocHost; //defined in GpuUtil.cpp
#endif

size_t cacheSize = 0; //target size of total allocations (unlimited if 0)
double cacheMargin = 0.5; //fraction of allocation size that may be wasted when reusing a cached pointer
bool cacheDebug = false;

void initMemCache()
{	if(not isGpuEnabled()) return; //Cache only used for GPU-related allocations
	
	//Cache size:
	const char* cacheSizeSrc = "JDFTX_CACHE_SIZE";
	const char* cacheSizeStr = getenv(cacheSizeSrc);
	string cacheTargetStr = "unlimited";
	if(not cacheSizeStr)
	{	cacheSizeSrc = "JDFTX_MEMPOOL_SIZE";
		cacheSizeStr = getenv(cacheSizeSrc);
		if(cacheSizeStr)
			logPrintf("WARNING: JDFTX_MEMPOOL_SIZE is deprecated; specify JDFTX_CACHE_SIZE instead.\n");
	}
	if(cacheSizeStr)
	{	int cacheSizeMB;
		if((sscanf(cacheSizeStr, "%d", &cacheSizeMB) == 1) and (cacheSizeMB >= 0))
		{	cacheSize = ((size_t)cacheSizeMB) << 20; //convert to bytes
			ostringstream oss; oss << cacheSizeMB << " MB (per process)";
			cacheTargetStr = oss.str();
		}
		else die("Invalid %s=\"%s\".\n", cacheSizeSrc, cacheSizeStr);
	}
	
	//Allocation margin:
	const char* cacheMarginStr = getenv("JDFTX_CACHE_MARGIN");
	if(cacheMarginStr)
		if(not ((sscanf(cacheMarginStr, "%lf", &cacheMargin) == 1) and (cacheMargin > 0)))
			die("Invalid JDFTX_CACHE_MARGIN=\"%s\".\n", cacheMarginStr);
	
	//Report cache parameters:
	logPrintf("Allocation cache: target size: %s, margin: %lg\n", 
		cacheTargetStr.c_str(), cacheMargin);

	//Allocation margin:
	const char* cacheDebugStr = getenv("JDFTX_CACHE_DEBUG");
	cacheDebug = cacheDebugStr and (string(cacheDebugStr) == "yes");
	if(cacheDebug)
	{	fprintf(stderr, "CACHEDBG: cache allocation debug messages enabled.\n");
		fflush(stderr);
	}
}


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
					logFlush();
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


//-------- Cache memory allocations to reduce system alloc/free calls ---------

namespace MemCache
{
	//Cache memory allocations in a memory space abstracted by MemSpace
	//MemSpace is a tag class with static members:
	// const bool shouldCache;    //whether to cache allocations for this type
	// void* alloc(size_t, bool); //returns 0 when out of memory (bool selects GPU/CPU for unified case)
	// void free(void*);          //assumed to not fail
	// void outOfMemory();        //exit with appropriate out of memory error
	template<typename MemSpace> class MemCache
	{	std::mutex lock; //for thread safety
		
		//Cache pointers:
		std::map<void*, size_t> allocated;
		std::multimap<size_t, void*> cache;
		typedef std::multimap<size_t, void*>::iterator CacheIter;
		size_t allocatedTot, cacheTot, overallTot; //total sizes per category
		size_t nAllocs, nAllocsRaw, nFrees, nFreesRaw; //statistics of total and raw alloc/free calls
	public:
		void free(void* ptr)
		{	if(not MemSpace::shouldCache) {	MemSpace::free(ptr); return; }
			nFrees++;
			//Move ptr from allocated to cache (for future allocs)
			lock.lock();
			auto iter = allocated.find(ptr);
			assert(iter != allocated.end());
			size_t size = iter->second;
			cache.insert(std::make_pair(size, ptr));
			cacheTot += size;
			allocated.erase(iter);
			allocatedTot -= size;
			lock.unlock();
		}
		
		void* alloc(size_t size, bool onGpu)
		{	if(not MemSpace::shouldCache) return MemSpace::alloc(size, onGpu);
			nAllocs++;
			lock.lock();
			//Search cache for suitable allocations first:
			{	auto iter = cache.lower_bound(size); //smallest allocation >= required
				if(iter != cache.end())
				{	size_t sizeCached = iter->first;
					void* ptrCached = iter->second;
					if(sizeCached < size + size*cacheMargin) //ignore if too big (avoid wasting memory)
					{	allocated.insert(std::make_pair(ptrCached, sizeCached));
						allocatedTot += sizeCached;
						cache.erase(iter);
						cacheTot -= sizeCached;
						lock.unlock();
						return ptrCached;
					}
				}
			}
			
			//Clear in advance to avoid overflowing target total memory usage:
			size_t expectedSize = overallTot + size;
			if(cacheSize and (expectedSize > cacheSize))
			{	const size_t expectedOverflow = expectedSize - cacheSize;
				report("Overflow anticipated", expectedOverflow);
				for(auto iter=cache.begin(); iter!=cache.end();)
				{	iter = freeRaw(iter);
					if(overallTot + size <= cacheSize) break;
				}
				report("Overflow addressed", expectedOverflow);
			}
			
			//Allocate from underlying memory space:
			void* ptr = allocRaw(size, onGpu);
			if(ptr)
			{	lock.unlock();
				return ptr;
			}
			report("Allocate failed", size);
			//Deallocate (unsuitable) cached entries till allocation succeeds:
			for(auto iter=cache.begin(); iter!=cache.end();)
			{	iter = freeRaw(iter);
				void* ptr = allocRaw(size, onGpu);
				if(ptr)
				{	lock.unlock();
					report("Allocate fixed", size);
					return ptr;
				}
			}
			//Allocation unsuccessful: out of memory even after emptying cache:
			lock.unlock();
			MemSpace::outOfMemory();
			return NULL;
		}
		
		MemCache()
		{	allocatedTot = 0;
			cacheTot = 0;
			overallTot = 0;
			nAllocs = 0;
			nAllocsRaw = 0;
			nFrees = 0;
			nFreesRaw = 0;
		}
		
		~MemCache()
		{	for(auto entry: allocated) MemSpace::free(entry.first);
			for(auto entry: cache) MemSpace::free(entry.second);
			allocated.clear();
			cache.clear();
		}
		
	private:
		//Raw allocation, along with book-keeping for cache
		inline void* allocRaw(size_t size, bool onGpu)
		{	void* ptr = MemSpace::alloc(size, onGpu);
			nAllocsRaw++; //counted regardless of success
			if(ptr)
			{	allocated.insert(std::make_pair(ptr, size));
				allocatedTot += size;
				overallTot += size;
			}
			return ptr;
		}
		
		//Free underlying memory of entry in cache, and advance iterator to next entry
		inline CacheIter freeRaw(CacheIter iter)
		{	MemSpace::free(iter->second);
			nFreesRaw++;
			cacheTot -= iter->first;
			overallTot -= iter->first;
			return cache.erase(iter);
		}
		
		inline void report(const char* context, size_t size)
		{	static const size_t halfMB = (1L << 19);
			if(cacheDebug)
			{	fprintf(stderr,
					"CACHEDBG: %s %lu MB; total: %lu MB, cache: %lu MB (%lu entries), raw allocs: %.1lf%%, frees: %.1lf%%\n",
					context, (size + halfMB) >> 20, (overallTot + halfMB) >> 20, (cacheTot + halfMB) >> 20,
					cache.size(), nAllocsRaw*100.0/nAllocs, nFreesRaw*100.0/nFrees);
				fflush(stderr);
			}
		}
	};
	
	//---- MemSpace classes for each memory space ----
	#if defined(GPU_ENABLED) && defined(CUDA_MANAGED_MEMORY)
	struct MemSpaceUnified
	{	static const bool shouldCache = true;
		static void* alloc(size_t size, bool onGpu)
		{	void* ptr;
			cudaError_t ret = cudaMallocManaged(&ptr, size);
			if(not onGpu) cudaDeviceSynchronize();
			if(prefetchSupported) cudaMemPrefetchAsync(ptr, size, onGpu ? memLocDevice : memLocHost);
			cudaGetLastError(); //clear error since handled by checking pointer
			return (ret==cudaSuccess) ? ptr : 0;
		}
		static void free(void* ptr) { cudaFree(ptr); }
		static void outOfMemory() die_alone("Managed memory allocation failed (out of memory)\n");
	};
	typedef MemSpaceUnified MemSpaceCPU;
	typedef MemSpaceUnified MemSpaceGPU;
	#elif defined(GPU_ENABLED) && defined(PINNED_HOST_MEMORY)
	struct MemSpaceCPU
	{	static const bool shouldCache = true;
		static void* alloc(size_t size, bool onGpu)
		{	void* ptr;
			cudaError_t ret = cudaMallocHost(&ptr, size);
			cudaGetLastError(); //clear error since handled by checking pointer
			return (ret==cudaSuccess) ? ptr : 0;
		}
		static void free(void* ptr) { cudaFreeHost(ptr); }
		static void outOfMemory() die_alone("Host memory allocation failed (out of pinned memory)\n");
	};
	#else
	struct MemSpaceCPU
	{	static const bool shouldCache = false;
		static void* alloc(size_t size, bool onGpu) { return fftw_malloc(size); }
		static void free(void* ptr) { fftw_free(ptr); }
		static void outOfMemory() die_alone("Memory allocation failed (out of memory)\n");
	};
	#endif
	#if defined(GPU_ENABLED) && (!defined(CUDA_MANAGED_MEMORY))
	struct MemSpaceGPU
	{	static const bool shouldCache = true;
		static void* alloc(size_t size, bool onGpu)
		{	assert(isGpuMine());
			void* ptr;
			cudaError_t ret = cudaMalloc(&ptr, size);
			cudaGetLastError(); //clear error since handled by checking pointer
			return (ret==cudaSuccess) ? ptr : 0;
		}
		static void free(void* ptr)
		{	assert(isGpuMine());
			cudaFree(ptr);
		}
		static void outOfMemory() die_alone("GPU memory allocation failed (out of memory)\n");
	};
	#endif
	
	//Cache accessor functions (to avoid file-level static variables):
	MemCache<MemSpaceCPU>& CPU() { static MemCache<MemSpaceCPU> cache; return cache; }
	#ifdef GPU_ENABLED
	#ifdef CUDA_MANAGED_MEMORY
	MemCache<MemSpaceGPU>& GPU() { return CPU(); } //use single combined cache
	#else
	MemCache<MemSpaceGPU>& GPU() { static MemCache<MemSpaceGPU> cache; return cache; }
	#endif
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
		      MemCache::GPU().free(c);
		#else
		assert(!"onGpu=true without GPU_ENABLED"); //Should never get here!
		#endif
	}
	else MemCache::CPU().free(c);
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
		c = MemCache::GPU().alloc(nBytes, true);
		#else
		assert(!"onGpu=true without GPU_ENABLED");
		#endif
	}
	else
	{	c = MemCache::CPU().alloc(nBytes, false);
		#if defined(GPU_ENABLED) && defined(CUDA_MANAGED_MEMORY)
		cudaDeviceSynchronize();
		#endif
	}
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
	#ifdef CUDA_MANAGED_MEMORY
	cudaDeviceSynchronize(); //finish pending computes, but no explicit data move
	if(prefetchSupported) cudaMemPrefetchAsync(c, nBytes, memLocHost);
	gpuErrorCheck();
	#else
	void* cCpu = MemCache::CPU().alloc(nBytes, false);
	cudaMemcpy(cCpu, me.c, nBytes, cudaMemcpyDeviceToHost);
	   MemCache::GPU().free(me.c); //Free GPU mem
	me.c = cCpu; //Make c a cpu pointer
	#endif
	me.onGpu = false;
#endif
}

// Move data to GPU
void ManagedMemoryBase::toGpu() const
{	if(onGpu || !c) return; //already on gpu, or no data
#ifdef GPU_ENABLED
	assert(isGpuMine());
	ManagedMemoryBase& me = *((ManagedMemoryBase*)this);
	#ifdef CUDA_MANAGED_MEMORY
	if(prefetchSupported) cudaMemPrefetchAsync(c, nBytes, memLocDevice);
	gpuErrorCheck();
	#else
	void* cGpu = MemCache::GPU().alloc(nBytes, true);
	cudaMemcpy(cGpu, me.c, nBytes, cudaMemcpyHostToDevice);
	   MemCache::CPU().free(me.c); //Free CPU mem
	me.c = cGpu; //Make c a gpu pointer
	#endif
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

#ifndef CUDATEST_H
#define CUDATEST_H

namespace CUDATest {
	const int SURFACE_SIZE       = 5120;
	const int CONFIGURATION_SIZE = 20;
	const int CONFIGURATIONS     = 256;
	const int TIMESTEPS          = 240;

	float measure_time();
}

#endif
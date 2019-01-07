#include "cuda_surface.cuh"

#include <cstdio>
#include <cuda.h>
#include <cmath>

void CUDA_Surface::Points3D::allocate()
{
    cudaMalloc((void**) &x, size*sizeof(float));
    cudaMalloc((void**) &y, size*sizeof(float));
    cudaMalloc((void**) &z, size*sizeof(float));
}

void CUDA_Surface::Points3D::allocate(const int _size)
{
    size = _size;
    cudaMalloc((void**) &x, size*sizeof(float));
    cudaMalloc((void**) &y, size*sizeof(float));
    cudaMalloc((void**) &z, size*sizeof(float));
}

void CUDA_Surface::Points3D::free()
{
    cudaFree(x);
    cudaFree(y);
    cudaFree(z);
}

void CUDA_Surface::Points3D::load_from(const CUDA_Surface::Points3D &cpu) const
{
    cudaMemcpy(x, cpu.x, size*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(y, cpu.y, size*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(z, cpu.z, size*sizeof(float), cudaMemcpyHostToDevice);
}

void CUDA_Surface::Points3D::save_to(const CUDA_Surface::Points3D &cpu) const
{
    cudaMemcpy(cpu.x, x, size*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu.y, y, size*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu.z, z, size*sizeof(float), cudaMemcpyDeviceToHost);
}

void CUDA_Surface::Points1D::allocate()
{
    cudaMalloc((void**) &s, size*sizeof(float));
    cudaMemset(s, 0, size);
}

void CUDA_Surface::Points1D::allocate(const int _size)
{
    size = _size;
    cudaMalloc((void**) &s, size*sizeof(float));
    cudaMemset(s, 0, size);
}

void CUDA_Surface::Points1D::free()
{
    cudaFree(s);
}

void CUDA_Surface::Points1D::load_from(const CUDA_Surface::Points1D &cpu) const
{
    cudaMemcpy(s, cpu.s, size*sizeof(float), cudaMemcpyHostToDevice);
}

void CUDA_Surface::Points1D::save_to(const CUDA_Surface::Points1D &cpu) const
{
    cudaMemcpy(cpu.s, s, size*sizeof(float), cudaMemcpyDeviceToHost);
}


// __global__ void compute_surface_layer(float *surface, const float *cx, const float *cy, const float *cz, 
//                                                       const float *x, const float *y, const float *z, const float hor)
// {
//     const int SURFACE_SIZE       = 5120;
//     const int CONFIGURATION_SIZE = 20;
//     const int CONFIGURATIONS     = 100;
//     const int TIMESTEPS          = 240;

//     const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
//     const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point
//     const int global_time  = blockDim.z*blockIdx.z + threadIdx.z; // time

//     float surface_result = 0.0, sat_x = 0.0, sat_y = 0.0, sat_z = 0.0; 
//     for (int iorb = 0; iorb < CONFIGURATION_SIZE; iorb++)
//     {
//         sat_x = x[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb];
//         sat_y = y[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb];
//         sat_z = z[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb];
//         // sat_h = h[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb]    
//         surface_result |= (sat_x*cx[global_point] + sat_y*cy[global_point] + sat_z*cz[global_point] > hor);
//     } 
//     surface[(global_time*SURFACE_SIZE + global_point)*CONFIGURATIONS + global_conf] = surface_result;
// }

__global__ void compute_surface_layer(float *surface, 
    const float *cx, const float *cy, const float *cz, 
    const float *x,  const float *y,  const float *z, const float hor)
{
    const int SURFACE_SIZE       = 5120;
    const int CONFIGURATION_SIZE = 20;
    const int CONFIGURATIONS     = 128;

    const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
    const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point
    const int global_time  = blockDim.z*blockIdx.z + threadIdx.z; // time

    float surface_result = 0.0, sat_x = 0.0, sat_y = 0.0, sat_z = 0.0; 
    for (int iorb = 0; iorb < CONFIGURATION_SIZE; iorb++)
    {
        sat_x = x[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
        sat_y = y[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
        sat_z = z[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
        surface_result += (sat_x*cx[global_point] + sat_y*cy[global_point] + sat_z*cz[global_point] > hor);
    } 
    surface[(global_time*SURFACE_SIZE + global_point)*CONFIGURATIONS + global_conf] = surface_result;
}

float loc_horizon(const float H, const float alpha)
{
    const float delta = 10.0 * M_PI / 180.0; // требуемое возвышение спутника над горизонтом
    const float alpha_star = asin(cos(delta) / H);

    return (alpha < alpha_star) ? cos(asin(H*sin(alpha)) - alpha) : sin(delta + alpha_star);
    // return sin(delta + alpha_star);
}

float CUDA_Surface::compute_surface(Points1D &gpu_surface, const Points3D &gpu_centroids, const Points3D &gpu_pos)
{
    const int SURFACE_SIZE       = 5120;
    // const int CONFIGURATION_SIZE = 20;
    const int CONFIGURATIONS     = 128;
    const int TIMESTEPS = 240;

    dim3 dim_grid, dim_block;

    dim_block.x = CONFIGURATIONS;
    dim_block.y = 1;
    dim_block.z = 1;

    dim_grid.x = 1;
    dim_grid.y = SURFACE_SIZE;
    dim_grid.z = TIMESTEPS;

    const float H = (6371.0 + 1500.0) / 6371.0;
    const float alpha = (120.0 * M_PI / 180.0) / 2.0;
    const float hor = loc_horizon(H, alpha);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    compute_surface_layer<<<dim_grid, dim_block>>>(gpu_surface.s, gpu_centroids.x, gpu_centroids.y, gpu_centroids.z,
                                                                  gpu_pos.x,       gpu_pos.y,       gpu_pos.z, hor);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    float gpu_time = 0.0; 
    cudaEventElapsedTime(&gpu_time, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return gpu_time;          
} 
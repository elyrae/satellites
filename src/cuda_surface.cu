#include "cuda_surface.cuh"

#include <cstdio>
#include <cuda.h>
#include <cmath>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
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

void CUDA_Surface::Points1D::allocate(const int _size)
{
    size = _size;
    cudaMalloc((void**) &s, size*sizeof(int));
    cudaMemset(s, 0, size);
}

void CUDA_Surface::Points1D::free()
{
    cudaFree(s);
}

void CUDA_Surface::Points1D::load_from(const CUDA_Surface::Points1D &cpu) const
{
    cudaMemcpy(s, cpu.s, size*sizeof(int), cudaMemcpyHostToDevice);
}

void CUDA_Surface::Points1D::save_to(const CUDA_Surface::Points1D &cpu) const
{
    cudaMemcpy(cpu.s, s, size*sizeof(int), cudaMemcpyDeviceToHost);
}


// __global__ void compute_surface_layer(float *surface, const float4 *c, const float4 *p, const float hor)
// {
//     const int SURFACE_SIZE       = 5120;
//     const int CONFIGURATION_SIZE = 20;
//     const int CONFIGURATIONS     = 128;

//     const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
//     const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point
//     const int global_time  = blockDim.z*blockIdx.z + threadIdx.z; // time

//     float surface_result = 0.0, sat_x = 0.0, sat_y = 0.0, sat_z = 0.0;
//     const float cen_x = c[global_point].x;
//     const float cen_y = c[global_point].y;
//     const float cen_z = c[global_point].z;

//     for (int iorb = 0; iorb < CONFIGURATION_SIZE; iorb++)
//     {
//         sat_x = p[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf].x;
//         sat_y = p[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf].y;
//         sat_z = p[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf].z;
//         surface_result += (sat_x*cen_x + sat_y*cen_y + sat_z*cen_z > hor);
//     } 
//     surface[(global_time*SURFACE_SIZE + global_point)*CONFIGURATIONS + global_conf] = surface_result;
// }

__global__ void compute_surface_layers(float *surface, 
    const float *cx, const float *cy, const float *cz, 
    const float *x,  const float *y,  const float *z, const float hor)
{
    const int SURFACE_SIZE       = 5120;
    const int CONFIGURATION_SIZE = 20;
    const int CONFIGURATIONS     = 256;

    const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
    const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point
    const int global_time  = blockDim.z*blockIdx.z + threadIdx.z; // time

    const float cen_x = cx[global_point];
    const float cen_y = cy[global_point];
    const float cen_z = cz[global_point];

    float surface_result = 0.0; 
    float sat_x = 0.0, sat_y = 0.0, sat_z = 0.0; 
    for (int iorb = 0; iorb < CONFIGURATION_SIZE; iorb++)
    {
        sat_x = x[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
        sat_y = y[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
        sat_z = z[(global_time*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
        surface_result += (sat_x*cen_x + sat_y*cen_y + sat_z*cen_z > hor);
    } 
    surface[(global_time*SURFACE_SIZE + global_point)*CONFIGURATIONS + global_conf] = surface_result;
}

// __global__ void compute_surface_layers_minimal(float *max_time, 
//                                                const float *cx, const float *cy, const float *cz, 
//                                                const float *x,  const float *y,  const float *z, const float hor)
// {
//     // const int SURFACE_SIZE       = 5120;
//     const int CONFIGURATION_SIZE = 20;
//     const int CONFIGURATIONS     = 256;
//     const int TIMESTEPS          = 240;

//     const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
//     const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point

//     const float cen_x = cx[global_point];
//     const float cen_y = cy[global_point];
//     const float cen_z = cz[global_point];

//     float surf = 0.0, sat_x = 0.0, sat_y = 0.0, sat_z = 0.0, m = 0.0, time = 0.0;
//     for (int timestep = 0; timestep < TIMESTEPS; timestep++) {
//         for (int iorb = 0; iorb < CONFIGURATION_SIZE; iorb++) {
//             sat_x = x[(timestep*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
//             sat_y = y[(timestep*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
//             sat_z = z[(timestep*CONFIGURATION_SIZE + iorb)*CONFIGURATIONS + global_conf];
//             surf += (sat_x*cen_x + sat_y*cen_y + sat_z*cen_z > hor);
//         }
//         m    = m*(1.0 - surf) + fmaxf(time, m)*surf;
//         time = (time + 30.0)*(1.0 - surf); 
        

//         // m    = (surf == 0.0) ? fmaxf(m, time) :       m;
//         // time = (surf == 0.0) ? 0.0            : (time + 30.0);
//     }

//     max_time[global_point*CONFIGURATIONS + global_conf] = fmaxf(m, time);
// }

__global__ void reduce_time(float *surface, float *max_time)
{
    const int SURFACE_SIZE       = 5120;
    const int CONFIGURATIONS     = 256;
    const int TIMESTEPS          = 240;

    const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
    const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point

    float surf = 0.0; 
    float time = 0.0, m = 0.0; 
    for (int timestep = 1; timestep < TIMESTEPS; timestep++) {
        surf = surface[(timestep*SURFACE_SIZE + global_point)*CONFIGURATIONS + global_conf];
        m    = (surf == 0.0) ? fmaxf(m, time) :       m;
        time = (surf == 0.0) ? 0.0            : (time + 30.0);  
    }
    m = fmaxf(m, time);
    max_time[global_point*CONFIGURATIONS + global_conf] = m;
}

__global__ void reduce_points(float *max_time)
{
    const int SURFACE_SIZE   = 5120;
    const int CONFIGURATIONS = 256;

    const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 

    float max = 0.0; 
    for (int point = 0; point < SURFACE_SIZE; point++)
        max = fmaxf(max, max_time[point*CONFIGURATIONS + global_conf]);
    max_time[global_conf] = max;
}

// // ЛЕВАЯ ГРАНИЦА БЕЗ УГЛОВЫХ ТОЧЕК, ПОДСЧЕТ КОЛИЧЕСТВА ЖИВЫХ СОСЕДЕЙ
// u[0, j] :=   u[1, j]   + u[1, j+1] + u[1, j-1] // ВНУТРИ ОБЛАСТИ
//            + u[0, j+1] + u[0, j-1]           // НА ГРАНИЦЕ
//            + u[N, j]   + u[N, j-1] + u[N, j+1] // БЕРУТСЯ С ПРАВОЙ ГРАНИЦЫ

// // ЛЕВЫЙ ВЕРХНИЙ УГОЛ, ПОДСЧЕТ КОЛИЧЕСТВА ЖИВЫХ СОСЕДЕЙ
// u[0, N] :=   u[1, N-1]                       // ВНУТРИ ОБЛАСТИ
//            + u[1, N]   + u[0, N-1]           // НА ГРАНИЦЕ
//            + u[N, N] +  // БЕРУТСЯ С ПРАВОЙ ГРАНИЦЫ

// max_time = max_time*(1 - surf) + max(time, max_time)*surf;
// time     = (time + settings.deltaT)*(1 - surf);  

// for (size_t j = 0; j < surface.size(); ++j)
//     max_time[j] =                 max_time[j]*(1 - surface[j]) + max(time[j], max_time[j])*surface[j];
// for (size_t j = 0; j < surface.size(); ++j)
//     time[j]     = (time[j] + settings.deltaT)*(1 - surface[j]) +                         0*surface[j]; 

float loc_horizon(const float H, const float alpha)
{
    const float delta = 10.0 * M_PI / 180.0; // требуемое возвышение спутника над горизонтом
    const float alpha_star = asin(cos(delta) / H);

    return (alpha < alpha_star) ? cos(asin(H*sin(alpha)) - alpha) : sin(delta + alpha_star);
    // return sin(delta + alpha_star);
}

float CUDA_Surface::compute_surface(Points1D &gpu_surface, Points1D &gpu_max_time, 
                                    const Points3D &gpu_centroids, const Points3D &gpu_pos)
{
    const int SURFACE_SIZE       = 5120;
    const int CONFIGURATION_SIZE = 20;
    const int CONFIGURATIONS     = 256;
    const int TIMESTEPS = 240;

    dim3    dim_grid_surf(1, SURFACE_SIZE, TIMESTEPS),    dim_block_surf(CONFIGURATIONS, 1, 1);
    dim3 dim_grid_maxtime(1, SURFACE_SIZE, 1),         dim_block_maxtime(CONFIGURATIONS, 1, 1);
    dim3  dim_grid_points(1, 1,            1),          dim_block_points(CONFIGURATIONS, 1, 1);
    // dim3 dim_grid_surf_min(1, SURFACE_SIZE, 1), dim_block_surf_min(CONFIGURATIONS, 1, 1);
    // dim3 dim_grid(CONFIGURATIONS, 1, 1), dim_block_surf(1, SURFACE_SIZE, TIMESTEPS);

    // dim_block.x = CONFIGURATIONS;
    // dim_block.y = 1;
    // dim_block.z = 1;

    // dim_grid.x = 1;
    // dim_grid.y = SURFACE_SIZE;
    // dim_grid.z = TIMESTEPS;

    const float H = (6371.0 + 1500.0) / 6371.0;
    const float alpha = (120.0 * M_PI / 180.0) / 2.0;
    const float hor = loc_horizon(H, alpha);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    compute_surface_layers<<<dim_grid_surf, dim_block_surf>>>(gpu_surface.s, 
                                                              gpu_centroids.x, gpu_centroids.y, gpu_centroids.z,
                                                              gpu_pos.x,       gpu_pos.y,       gpu_pos.z, hor);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    float gpu_time = 0.0; 
    cudaEventElapsedTime(&gpu_time, start, stop);
    printf("surf = %f ms\n", gpu_time);

    cudaEventRecord(start, 0);
    reduce_time<<<dim_grid_maxtime, dim_block_maxtime>>>(gpu_surface.s, gpu_max_time.s);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    cudaEventElapsedTime(&gpu_time, start, stop);
    printf("time reduce = %f ms\n", gpu_time);
    
    cudaEventRecord(start, 0);
    reduce_points<<<dim_grid_points, dim_block_points>>>(gpu_max_time.s);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );    

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    cudaEventElapsedTime(&gpu_time, start, stop);
    printf("point reduce = %f ms\n", gpu_time);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return gpu_time;          
} 

    // compute_surface_layers<<<dim_grid_surf, dim_block_surf>>>(gpu_surface.s, gpu_centroids.x, gpu_pos.x, hor);
    
// compute_surface_layers_minimal<<<dim_grid_surf_min, dim_block_surf_min>>>(gpu_max_time.s, 
//                                                                           gpu_centroids.x, gpu_centroids.y, gpu_centroids.z,
//                                                                           gpu_pos.x,       gpu_pos.y,       gpu_pos.z, hor);

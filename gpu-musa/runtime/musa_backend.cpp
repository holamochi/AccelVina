// The upstream host issues a fixed sequence of OpenCL commands. Implement only
// the subset it uses, backed by one in-order MUSA stream and two AOT kernels.
#include "wrapcl.h"
#include "vina_musa_launch.h"
#include <musa_runtime_api.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

struct _cl_kernel {
    unsigned number;
    std::array<cl_mem, 8> buffers{};
    std::array<int, 4> scalars{};
};
struct _cl_mem { void* ptr; size_t bytes; };
struct _cl_event { musaEvent_t done{}, begin{}; unsigned number; };

namespace {
musaDeviceProp device_properties;
_cl_kernel grid_kernel{1}, docking_kernel{2};
size_t live_bytes, peak_bytes, uploaded, downloaded;
FILE* profile;
void mark(const char* stage) {
    if (!profile) return;
    using clock = std::chrono::steady_clock;
    static const auto start = clock::now();
    std::fprintf(profile, "host_ms %.3f %s\n",
                 std::chrono::duration<double, std::milli>(clock::now()-start).count(), stage);
    std::fflush(profile);
}
void checked(musaError_t code) {
    if (code != musaSuccess) {
        std::fprintf(stderr, "MUSA: %s\n", musaGetErrorString(code));
        std::exit(EXIT_FAILURE);
    }
}
void finish_profile() {
    if (profile) {
        std::fprintf(profile,"peak_buffer_bytes %zu\nh2d_bytes %zu\nd2h_bytes %zu\n", peak_bytes, uploaded, downloaded);
        std::fclose(profile);
    }
}
}

void checkErr(cl_int err) {
    if (err != CL_SUCCESS) {
        std::fprintf(stderr, "Invalid docking backend command (%d)\n", err);
        std::exit(EXIT_FAILURE);
    }
}
void SetupPlatform(cl_platform_id** p, cl_int* i) {
    static cl_platform_id selected;
    *p = &selected;
    *i = 0;
}
void SetupDevice(cl_platform_id*, cl_device_id** ids, cl_int) {
    int count;
    checked(musaGetDeviceCount(&count));
    if (count != 1) {
        std::fprintf(stderr, "Expose one device with MUSA_VISIBLE_DEVICES (visible: %d)\n", count);
        std::exit(EXIT_FAILURE);
    }
    checked(musaSetDevice(0));
    checked(musaGetDeviceProperties(&device_properties, 0));
    char bus[32];
    checked(musaDeviceGetPCIBusId(bus, sizeof(bus), 0));
    if (const char* expected = std::getenv("VINA_EXPECT_PCI")) {
        if (std::strcmp(expected, bus) != 0) {
            std::fprintf(stderr,"Device PCI mismatch: expected %s, found %s\n", expected, bus);
            std::exit(EXIT_FAILURE);
        }
    }
    std::printf("\nGPU Device: %s (%s)\n", device_properties.name, bus);
    static cl_device_id selected;
    *ids = &selected;
    if (const char* path = std::getenv("VINA_MUSA_PROFILE")) {
        profile = std::fopen(path, "w");
        if (!profile) { std::perror(path); std::exit(EXIT_FAILURE); }
        std::atexit(finish_profile);
        mark("device_ready");
    }
}
// Context, queue and program handles are only passed through in this AOT backend.
void SetupContext(cl_platform_id*, cl_device_id*, cl_context* c, cl_uint, cl_int) { *c = nullptr; }
void SetupQueue(cl_command_queue* q, cl_context, cl_device_id*) { *q = nullptr; }
cl_program SetupBuildProgramWithBinary(cl_context, cl_device_id*, const char*) { return nullptr; }
cl_int clUnloadPlatformCompiler(cl_platform_id) { return CL_SUCCESS; }
void SetupKernel(cl_kernel* k, cl_program*, size_t n, const char names[][50]) {
    if (n != 2 || std::strcmp(names[0], "kernel1") || std::strcmp(names[1], "kernel2"))
        checkErr(CL_INVALID_KERNEL_NAME);
    k[0] = &grid_kernel;
    k[1] = &docking_kernel;
}
cl_int clGetDeviceInfo(cl_device_id, cl_device_info key, size_t bytes, void* out, size_t*) {
    auto* values = static_cast<size_t*>(out);
    if (key == CL_DEVICE_MAX_WORK_GROUP_SIZE && bytes == sizeof(size_t)) {
        values[0] = device_properties.maxThreadsPerBlock;
    } else if (key == CL_DEVICE_MAX_WORK_ITEM_SIZES && bytes == 3*sizeof(size_t)) {
        for (int i = 0; i < 3; ++i) values[i] = device_properties.maxThreadsDim[i];
    } else return CL_INVALID_VALUE;
    return CL_SUCCESS;
}
void CreateDeviceBuffer(cl_mem* handle, cl_mem_flags, size_t size, cl_context) {
    void* ptr;
    checked(musaMalloc(&ptr, size));
    *handle = new _cl_mem{ptr, size};
    live_bytes += size;
    peak_bytes = std::max(peak_bytes, live_bytes);
}
cl_int clReleaseMemObject(cl_mem h) {
    checked(musaFree(h->ptr));
    live_bytes -= h->bytes;
    delete h;
    return CL_SUCCESS;
}
cl_int clEnqueueWriteBuffer(cl_command_queue, cl_mem h, cl_bool blocking, size_t offset, size_t bytes,
                            const void* src, cl_uint n, const cl_event*, cl_event* event) {
    if (n || event || offset > h->bytes || bytes > h->bytes-offset) return CL_INVALID_VALUE;
    checked(musaMemcpyAsync(static_cast<char*>(h->ptr)+offset, src, bytes, musaMemcpyHostToDevice, nullptr));
    if (blocking) checked(musaStreamSynchronize(nullptr));
    uploaded += bytes;
    return CL_SUCCESS;
}
cl_int clEnqueueReadBuffer(cl_command_queue, cl_mem h, cl_bool blocking, size_t offset, size_t bytes,
                           void* dst, cl_uint n, const cl_event*, cl_event* event) {
    if (n || event || offset > h->bytes || bytes > h->bytes-offset) return CL_INVALID_VALUE;
    checked(musaMemcpyAsync(dst, static_cast<char*>(h->ptr)+offset, bytes, musaMemcpyDeviceToHost, nullptr));
    if (blocking) checked(musaStreamSynchronize(nullptr));
    downloaded += bytes;
    return CL_SUCCESS;
}
cl_int clFinish(cl_command_queue) { checked(musaStreamSynchronize(nullptr)); return CL_SUCCESS; }
void SetKernelArg(cl_kernel k, cl_uint idx, size_t bytes, const void* src) {
    unsigned buffers = k->number == 1 ? 7 : 8;
    if (idx < buffers && bytes == sizeof(cl_mem)) {
        std::memcpy(&k->buffers[idx], src, bytes);
    } else if (idx >= buffers && idx < (k->number == 1 ? 9u : 12u) && bytes == sizeof(int)) {
        std::memcpy(&k->scalars[idx-buffers], src, bytes);
    } else checkErr(CL_INVALID_ARG_INDEX);
}
cl_int clEnqueueNDRangeKernel(cl_command_queue, cl_kernel k, cl_uint dims, const size_t* offset,
                             const size_t* global, const size_t* local, cl_uint n, const cl_event*, cl_event* out) {
    const size_t g1[]{128,128,128}, l1[]{4,4,4}, g2[]{512,32}, l2[]{16,2};
    bool grid = k->number == 1;
    if (!out || offset || n || dims != (grid ? 3u : 2u) ||
        !std::equal(global, global+dims, grid ? g1 : g2) || !std::equal(local, local+dims, grid ? l1 : l2))
        return CL_INVALID_WORK_GROUP_SIZE;
    mark(grid ? "submit_grid" : "submit_search");
    auto* ev = new _cl_event{};
    ev->number = k->number;
    if (profile) {
        checked(musaEventCreate(&ev->begin));
        checked(musaEventRecord(ev->begin, nullptr));
    }
    checked(musaEventCreateWithFlags(&ev->done, profile ? musaEventDefault : musaEventDisableTiming));
    auto b = [&](int i) { return k->buffers[i]->ptr; };
    const auto& s = k->scalars;
    int status = grid
        ? vina_musa_launch_kernel1(b(0),b(1),b(2),b(3),b(4),b(5),b(6),s[0],s[1],nullptr)
        : vina_musa_launch_kernel2(b(0),b(1),b(2),b(3),b(4),b(5),b(6),b(7),s[0],s[1],s[2],s[3],nullptr);
    checked(static_cast<musaError_t>(status));
    checked(musaEventRecord(ev->done, nullptr));
    mark(grid ? "grid_queued" : "search_queued");
    *out = ev;
    return CL_SUCCESS;
}
cl_int clWaitForEvents(cl_uint n, const cl_event* events) {
    for (cl_uint i=0; i<n; ++i) checked(musaEventSynchronize(events[i]->done));
    mark("kernel_completed");
    return CL_SUCCESS;
}
cl_int clReleaseEvent(cl_event ev) {
    if (profile) {
        float ms;
        checked(musaEventSynchronize(ev->done));
        checked(musaEventElapsedTime(&ms, ev->begin, ev->done));
        std::fprintf(profile,"kernel%u_ms %.6f\n", ev->number, ms);
        checked(musaEventDestroy(ev->begin));
    }
    checked(musaEventDestroy(ev->done));
    delete ev;
    return CL_SUCCESS;
}

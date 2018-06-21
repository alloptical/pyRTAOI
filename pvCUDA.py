# convert by gpu
import pycuda.autoinit
import pycuda.driver as cuda
from pycuda.compiler import SourceModule


mod = SourceModule("""

__global__ void sample_mean(short* matrix, int pixelsPerLine, 
				int linesPerFrame, int samplesPerPixel, int flipEvenRows, long* result)
{
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	int result_idx = 0;
	int col = 0;
	int num_sample = 0;
	int this_value = 0;
	
	if (idx<pixelsPerLine*linesPerFrame*samplesPerPixel - samplesPerPixel+1){
		if ((idx - (idx / samplesPerPixel)*samplesPerPixel) == 0){
			result_idx = idx / samplesPerPixel;
			col = result_idx - (result_idx / pixelsPerLine)*pixelsPerLine;
			if ((result_idx / pixelsPerLine) - ((result_idx / pixelsPerLine) / 2) * 2 != flipEvenRows){
				result_idx = result_idx + pixelsPerLine - 2 * col - 1;
			}

			for (int i = 0; i < samplesPerPixel; i++){
				if (matrix[idx + i]>8192){
					this_value += matrix[idx + i] - 8192;
					num_sample += 1;
				}
			}

			if (num_sample>0){ result[result_idx] = this_value / num_sample; }
		}
	}
}

""")

# get cuda func
sample_mean = mod.get_function('sample_mean')

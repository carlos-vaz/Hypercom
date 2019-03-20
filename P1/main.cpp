#include <stdlib>
#include <time.h>
#include "surface.h"

#define ERROR_THRESH  1e-12	// Max tolerable error (infinity norm)
#define UPDATE_PERIOD 1000	// Display error after this many grid updates

extern void VTK_out(const int N, const int M, const double *Xmin, const double *Xmax,
             const double *Ymin, const double *Ymax, const double *T,
             const int index);	// be sure to comment out `main` in vtk.c if using this line

int main(void) {
	int N = 21, M = 21;
	double X = 2, Y = 1;
	Surface S(X, Y, N, M);
	S.setExpBoundaryT();	// Set boundary temperatures
	//S.setConstBoundaryT(0);
	S.setExpS();		// Set source terms
	//S.setConstS(-0.2);
	S.setExpSolution();	// Set analytic solution (for error calculation)
	S.printS();		// Print source terms
	S.printT();		// print initial temperatures

	Surface::internal_iter iter = S.begin();
	double err = ERROR_THRESH+1; long long count = 0; long rounds = 0;
	// LOOP TIL CONVERGENCE
	time_t start = clock();
	while(err > ERROR_THRESH)
	{
		*iter = S.getnewT(iter);
		if(++count%(S.internal_size*UPDATE_PERIOD)==0)
		{
			err = S.queryErrorT_relative();
			std::cout << "Error this round: " << err << std::endl;
			rounds++;
		}
		iter++;
	}
	time_t elapsed = clock() - start;


	// ##########	HELPFUL INFORMATION NOW PRINTED   ##########
	std::cout << "\nCALCULATED SOLUTION: \n"; S.printT();
	if(S.solution_provided) {std::cout << "\nANALYTICAL SOLUTION: \n"; S.printSolution();}
	std::cout << "\nNumber of cell updates: " << count << " (" << rounds << " rounds)" << std::endl;
	
	// Error analysis
	if(S.solution_provided) {
		double error = S.queryErrorT_absolute();
		std::cout << "ERROR (inf norm): " << error << std::endl;
	}

	// .VTK output (for ParaView)
	double *T = (double *)malloc(M * N * sizeof(double));
	for(int i=0; i<M*N; i++)
		T[i]=S.getT(i%N, i/N);
	double Xmin=0, Xmax=X, Ymin=0, Ymax=Y;
	VTK_out(N, M, &Xmin, &Xmax, &Ymin, &Ymax, T, 1);

	// Output value of MIDDLE POINT of region
	std::cout << "MIDDLE POINT value: " << S.getT(M/2,N/2) << std::endl;
	// Output TIME ELAPSED
	std::cout << "TIME ELAPSED: " << (double)elapsed/CLOCKS_PER_SEC << " seconds\n";
	return 0;
}


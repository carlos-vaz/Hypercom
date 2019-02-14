#include <iostream>
#include <math.h>
#include <assert.h>

class Surface {
	int numx;
	int numy;
	double dimx; 
	double dimy;
	double *Tdata; 		// temperature field
	double *Sdata; 		// heat source/sink field
	double *Tsolution;	// if provided

   public:
	bool solution_provided = false;
	long internal_size;
	class internal_iter : std::iterator<
					std::forward_iterator_tag,
					double,		// value type
					int, 		// difference type
					double*, 	// pointer type
					double		// reference type
						>{
	   public:
		int index;
		explicit internal_iter(Surface& s) : MySurface(s) {index = MySurface.numx+1;}
		internal_iter& operator++() {
			index++;
			if(index<MySurface.numx)
				index=MySurface.numx+1;
			else if(index%MySurface.numx==MySurface.numx-1)
				index+=2;
			else if(index%MySurface.numx==0)
				index+=1;
			if(index>MySurface.numx*(MySurface.numy-1)-2)
				index=MySurface.numx+1;
			assert(index>MySurface.numx && index%MySurface.numx!=MySurface.numx-1 && index%MySurface.numx!=0 && index<MySurface.numx*(MySurface.numy-1)-1);
			return *this;
		}
		internal_iter operator++(int) {return ++(*this);}
		bool operator==(internal_iter other) const {return other.index==index;}
		bool operator!=(internal_iter other) const {return !(other.index==index);}
		double& operator*() const {return MySurface.Tdata[index];}

	   private: 
		Surface& MySurface;
		
	};

	internal_iter begin() {return internal_iter(*this);}
	internal_iter end() {internal_iter endIter(*this); endIter.index = numx*(numy-1)-2; return endIter;}
	

	Surface(double X, double Y, int x, int y) : numx(x), numy(y), dimx(X), dimy(Y)	{
						assert(numx>1 && numy>1);
						Tdata 	  = new double[numx*numy];
						Sdata 	  = new double[numx*numy];
						Tsolution = new double[numx*numy];
						internal_size = (numx-1)*(numy-1);	}
	void setT(double, int, int);
	void setS(double, int, int);
	double getT(int, int);
	double getS(internal_iter);
	double getSolution(internal_iter);

	void setExpBoundaryT(); // set boundaries accord. to xe^y
	void setConstBoundaryT(double);
	void setExpS();
	void setConstS(double);
	void setExpSolution();
	void setConstSolution(double);

	double laplace(int, int);
	double laplace(internal_iter);
	double getnewT(internal_iter);

	double queryErrorT_relative();
	double queryErrorT_absolute();
	
	void printT();
	void printS();
	void printSolution();


	// DEBUGGING ONLY
	int getnumx() {return numx;}

};

void Surface::setT(double val, int x, int y) {
	assert(x>=0 && x<numx && y>=0 && y<numy);
	Tdata[y*numx+x] = val;
}

void Surface::setS(double val, int x, int y) {
	assert(x>=0 && x<numx && y>=0 && y<numy);
	Sdata[y*numx+x] = val;
}

double Surface::getT(int x, int y) {
	return Tdata[y*numx+x];
}

double Surface::getS(internal_iter iter) {
	return Sdata[iter.index];
}

double Surface::getSolution(internal_iter iter) {
	return Tsolution[iter.index];
}

void Surface::printT() {
	std::cout << "\n";
	for(int i=0; i<numy; i++)
	{
		for(int j=0; j<numx; j++)
			std::cout << Tdata[i*numx+j] << ", ";	
		std::cout << "\n";
	}
}

void Surface::printS() {
	std::cout << "\n";
	for(int i=0; i<numy; i++)
	{
		for(int j=0; j<numx; j++)
			std::cout << Sdata[i*numx+j] << ", ";	
		std::cout << "\n";
	}
}

void Surface::printSolution() {
	std::cout << "\n";
	for(int i=0; i<numy; i++)
	{
		for(int j=0; j<numx; j++)
			std::cout << Tsolution[i*numx+j] << ", ";	
		std::cout << "\n";
	}
}


void Surface::setExpBoundaryT() {
	for(int i=0; i<numy; i++)
		for(int j=0; j<numx; j++)
		{
			if(i!=0 && i!=numy-1 && j!=0 && j!=numx-1)
				j+=std::max(numx-2, 0);
			double x = dimx*(double)j/(double)(numx-1);
			double y = dimy*(double)i/(double)(numy-1);
			Tdata[i*numx+j] = x*pow(2.71828, y);
		}
}

void Surface::setConstBoundaryT(double val) {
	for(int i=0; i<numy; i++)
		for(int j=0; j<numx; j++)
		{
			if(i!=0 && i!=numy-1 && j!=0 && j!=numx-1)
				j+=std::max(numx-2, 0);
			Tdata[i*numx+j] = val;
		}
}


void Surface::setExpS() {
	for(int i=0; i<numy; i++)
		for(int j=0; j<numx; j++)
		{
			double x = dimx*(double)j/(double)(numx-1);
			double y = dimy*(double)i/(double)(numy-1);
			Sdata[i*numx+j] = x*pow(2.71828, y);
		}
}

void Surface::setConstS(double val) {
	for(int i=0; i<numy; i++)
		for(int j=0; j<numx; j++)
			Sdata[i*numx+j] = val;
}

void Surface::setExpSolution() {
	for(int i=0; i<numy; i++)
		for(int j=0; j<numx; j++)
		{
			double x = dimx*(double)j/(double)(numx-1);
			double y = dimy*(double)i/(double)(numy-1);
			Tsolution[i*numx+j] = x*pow(2.71828, y);
		}
	solution_provided = true;
}

void Surface::setConstSolution(double val) {
	for(int i=0; i<numy; i++)
		for(int j=0; j<numx; j++)
		{
			Tsolution[i*numx+j] = val;
		}
	solution_provided = true;
}

double Surface::laplace(int x, int y) {
	assert(x>0 && x<numx-1 && y>0 && y<numy-1);
	double delx = dimx/(double)(numx-1);
	double dely = dimy/(double)(numy-1);
	return (getT(x-1,y)+getT(x+1,y)-2*getT(x,y))/pow(delx,2)+(getT(x,y-1)+getT(x,y+1)-2*getT(x,y))/pow(dely,2);
}

double Surface::laplace(internal_iter iter) {
	int x = iter.index%numx;
	int y = (iter.index-x)/numx;
	assert(x>0 && x<numx-1 && y>0 && y<numy-1);
	double delx = dimx/(double)(numx-1);
	double dely = dimy/(double)(numy-1);
	return (getT(x-1,y)+getT(x+1,y)-2*getT(x,y))/pow(delx,2)+(getT(x,y-1)+getT(x,y+1)-2*getT(x,y))/pow(dely,2);
}

double Surface::getnewT(internal_iter iter) {
	int x = iter.index%numx;
	int y = (iter.index-x)/numx;
	assert(x>0 && x<numx-1 && y>0 && y<numy-1);
	double delx = dimx/(double)(numx-1);
	double dely = dimy/(double)(numy-1);
	return (-1*getS(iter)*pow((delx*dely),2)+(getT(x-1,y)+getT(x+1,y))*pow(dely,2)+(getT(x,y-1)+getT(x,y+1))*pow(delx,2))/(2*pow(delx,2)+2*pow(dely,2));
}

double Surface::queryErrorT_relative() {
	internal_iter iter = begin();
	double maxErr = 0, thisErr;
	for(; true;)
	{
		thisErr = fabs(*iter - getnewT(iter));
		if(thisErr > maxErr)
			maxErr = thisErr;
		if(++iter==begin())
			break;
	}
	return maxErr;
}

double Surface::queryErrorT_absolute() {
	internal_iter iter = begin();
	double maxErr = 0, thisErr;
	for(; true;)
	{
		thisErr = fabs(*iter - getSolution(iter));
		if(thisErr > maxErr)
			maxErr = thisErr;
		if(++iter==begin())
			break;
	}
	return maxErr;
}


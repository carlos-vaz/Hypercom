#include "surface.h"

int main(void) {
	Surface S(2, 1, 11, 11);
	S.setExpBoundaryT();
	S.setExpS();
	S.printT();
	S.printS();

	Surface::internal_iter iter0(S);
	std::cout << "iter0 index: " << iter0.index << std::endl;
	auto iter1 = S.begin();
	std::cout << "iter1 index: " << iter1.index << std::endl;
	auto iter2 = S.end();
	std::cout << "iter2 index: " << iter2.index << std::endl;
	

	int i=8;
	int& r1 = i;
	r1 = 9;
	std::cout << "r1=" << r1 << ", i=" << i << std::endl;
	int &r2 = i;
	r2 = 10;
	std::cout << "r2=" << r2 << ", r1=" << r1 << ", i=" << i << std::endl;


	return 0;
}

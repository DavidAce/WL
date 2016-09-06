
#include "personality.hpp"


using namespace std;
using namespace EMC_constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, personality const &guy) {
	for (int i = 0; i < nGenes; i++) {
		os << guy.H << endl;
	}
	return os;
}



